#!/usr/bin/env python3
"""
13_windowed_embed.py
====================
Windowed ESM-2 re-embedding of the sequences that step 02 had to TRUNCATE.
(Reviewer 1 Comment 3, Option A: targeted windowed analysis of the truncated
anchors and of the long NLR / LRR-RK benchmark proteins.)

Step 02 truncates every sequence longer than `max_seq_len` (1,022 aa) to its
first 1,022 residues, which removes the C-terminal LRR region of the longest
immune receptors. This script re-embeds ONLY those long sequences, by tiling
each one into overlapping windows of `max_seq_len` residues and mean-pooling
each window separately. Short sequences are untouched: their step-02
embeddings are already exact.

The final window in every protein is anchored to the C-terminus
(start = L - window), so the C-terminal LRR region is always covered by
exactly one complete window regardless of the stride.

Nothing is collapsed here. All per-window embeddings are written out, so that
14_truncation_sensitivity.py can compare alternative pooling rules
(length-weighted mean over windows vs. max cosine similarity over windows)
without any further GPU time.

Run from the project root, once per species:
    python 13_windowed_embed.py --config config_arabidopsis.yaml
    python 13_windowed_embed.py --config config_rice.yaml
    python 13_windowed_embed.py --config config_grapevine.yaml

Outputs -> results/<species>/13_windowed_embed/
    windowed_proteome_embeddings.npz  - arrays: embeddings (MxD), protein_ids,
                                        window_index, win_start, win_end
    windowed_anchor_embeddings.npz    - same, for the truncated anchors
    windows_manifest.csv              - one row per long protein
    summary.yaml                      - counts, truncation fractions
"""

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent))

from shared import load_config, step_dir, get_logger

# Non-standard residues stripped in steps 00 and 02; repeated here so that
# window coordinates refer to the same cleaned sequence the pipeline scored.
STRIP_CHARS = "*XJBZUO"


# ── reuse the step-02 embedder verbatim (module name starts with a digit) ────

def _load_step02():
    path = Path(__file__).resolve().parent / "02_embed_proteins.py"
    if not path.exists():
        raise FileNotFoundError(f"Cannot find {path} (run from the project root)")
    spec = importlib.util.spec_from_file_location("_step02_embed", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ── sequence handling ───────────────────────────────────────────────────────

def load_fasta(path: Path) -> dict:
    """FASTA -> {id: cleaned sequence}. Identical to step 02's loader."""
    seqs = {}
    for rec in SeqIO.parse(str(path), "fasta"):
        seq = str(rec.seq)
        for ch in STRIP_CHARS:
            seq = seq.replace(ch, "")
        key = rec.id.split("|")[0] if "|" in rec.id else rec.id
        seqs[key] = seq
    return seqs


def window_bounds(length: int, window: int, stride: int):
    """
    Tile [0, length) into windows of `window` residues, advancing by `stride`.
    The last window is always flush with the C-terminus.
    Returns a list of (start, end) with end - start == window.
    """
    if length <= window:
        return [(0, length)]
    starts = list(range(0, length - window + 1, stride))
    if starts[-1] != length - window:
        starts.append(length - window)
    return [(s, s + window) for s in starts]


def build_window_dict(seqs: dict, window: int, stride: int):
    """
    Returns:
        win_seqs  {window_id: subsequence}   for proteins longer than `window`
        rows      list of manifest dicts     one per long protein
        coords    {window_id: (protein_id, k, start, end)}
    """
    win_seqs, rows, coords = {}, [], {}
    for pid, seq in seqs.items():
        L = len(seq)
        if L <= window:
            continue
        bounds = window_bounds(L, window, stride)
        for k, (s, e) in enumerate(bounds):
            wid = f"{pid}__w{k}"
            win_seqs[wid] = seq[s:e]
            coords[wid] = (pid, k, s, e)
        rows.append({
            "protein_id": pid,
            "length": L,
            "n_windows": len(bounds),
            "last_window_start": bounds[-1][0],
            "residues_lost_by_truncation": L - window,
        })
    return win_seqs, rows, coords


def save_windowed(out_path: Path, emb, win_ids, coords, logger):
    """Write per-window embeddings plus their coordinates."""
    prot = np.array([coords[w][0] for w in win_ids])
    kidx = np.array([coords[w][1] for w in win_ids], dtype=np.int32)
    start = np.array([coords[w][2] for w in win_ids], dtype=np.int32)
    end = np.array([coords[w][3] for w in win_ids], dtype=np.int32)
    np.savez_compressed(
        out_path,
        embeddings=emb.astype(np.float32),
        window_ids=np.array(win_ids),
        protein_ids=prot,
        window_index=kidx,
        win_start=start,
        win_end=end,
    )
    logger.info(f"  wrote {emb.shape[0]:,} window embeddings -> {out_path.name}")


# ── main ────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Windowed ESM-2 re-embedding of truncated (>max_seq_len) proteins"
    )
    ap.add_argument("--config", default="config.yaml",
                    help="pipeline config (default: config.yaml)")
    ap.add_argument("--stride", type=int, default=None,
                    help="window stride in residues (default: max_seq_len // 2)")
    ap.add_argument("--batch-size", type=int, default=None,
                    help="override config batch_size (windows are all full length)")
    ap.add_argument("--limit", type=int, default=None,
                    help="smoke test: embed only the first N long proteins")
    args = ap.parse_args()

    cfg = load_config(args.config)
    out = step_dir(cfg, "13_windowed_embed")
    logger = get_logger("13_windowed_embed", out)

    window = int(cfg["max_seq_len"])
    stride = args.stride if args.stride else window // 2
    batch_size = args.batch_size if args.batch_size else cfg["batch_size"]

    logger.info("=" * 65)
    logger.info("STEP 13  - Windowed re-embedding of truncated sequences")
    logger.info("=" * 65)
    logger.info(f"Config  : {args.config}  (species: {cfg.get('species')})")
    logger.info(f"Window  : {window} aa   Stride: {stride} aa")

    base = Path(cfg["base_output_dir"])
    proteome_fasta = base / "00_download_proteome" / "proteome.fasta"
    anchor_fasta = base / "01_fetch_anchors" / "anchors.fasta"
    for p in (proteome_fasta, anchor_fasta):
        if not p.exists():
            raise FileNotFoundError(f"Missing {p}. Run steps 00-02 first.")

    proteome_seqs = load_fasta(proteome_fasta)
    anchor_seqs = load_fasta(anchor_fasta)
    logger.info(f"Proteome: {len(proteome_seqs):,} proteins")
    logger.info(f"Anchors : {len(anchor_seqs)} proteins")

    # ── identify the truncated sequences ─────────────────────────────────
    p_win, p_rows, p_coords = build_window_dict(proteome_seqs, window, stride)
    a_win, a_rows, a_coords = build_window_dict(anchor_seqs, window, stride)

    n_long_p, n_long_a = len(p_rows), len(a_rows)
    frac_p = n_long_p / len(proteome_seqs) if proteome_seqs else 0.0
    logger.info(f"\nTruncated in step 02:")
    logger.info(f"  proteome: {n_long_p:,} / {len(proteome_seqs):,} ({frac_p:.2%}) "
                f"-> {len(p_win):,} windows")
    logger.info(f"  anchors : {n_long_a} / {len(anchor_seqs)} "
                f"-> {len(a_win):,} windows")
    if n_long_p == 0:
        logger.warning("No sequences exceed max_seq_len; nothing to do.")
        return

    if args.limit:
        keep = {r["protein_id"] for r in p_rows[: args.limit]}
        p_win = {w: s for w, s in p_win.items() if p_coords[w][0] in keep}
        p_rows = [r for r in p_rows if r["protein_id"] in keep]
        logger.warning(f"--limit {args.limit}: SMOKE TEST, "
                       f"{len(p_win):,} windows only. Do not use for the paper.")

    # ── manifest ─────────────────────────────────────────────────────────
    manifest = pd.DataFrame(
        [{**r, "source": "proteome"} for r in p_rows]
        + [{**r, "source": "anchor"} for r in a_rows]
    ).sort_values(["source", "length"], ascending=[True, False])
    manifest.to_csv(out / "windows_manifest.csv", index=False)

    # ── embed ────────────────────────────────────────────────────────────
    step02 = _load_step02()
    embedder = step02.ESM2Embedder(
        model_name=cfg["esm_model"],
        batch_size=batch_size,
        max_seq_len=window,     # windows are already <= window, so no truncation
        device=cfg["device"],
        logger=logger,
    )

    if a_win:
        logger.info("\nEmbedding truncated ANCHOR windows ...")
        a_emb, a_ids = embedder.embed(a_win, out / "_cache_anchor_windows.npz")
        save_windowed(out / "windowed_anchor_embeddings.npz",
                      a_emb, a_ids, a_coords, logger)

    logger.info("\nEmbedding truncated PROTEOME windows ...")
    p_emb, p_ids = embedder.embed(p_win, out / "_cache_proteome_windows.npz")
    save_windowed(out / "windowed_proteome_embeddings.npz",
                  p_emb, p_ids, p_coords, logger)

    # ── summary ──────────────────────────────────────────────────────────
    lengths = manifest.loc[manifest["source"] == "proteome", "length"]
    summary = {
        "species": cfg.get("species"),
        "window": window,
        "stride": stride,
        "esm_model": cfg["esm_model"],
        "n_proteome_total": int(len(proteome_seqs)),
        "n_proteome_truncated": int(n_long_p),
        "frac_proteome_truncated": round(float(frac_p), 4),
        "n_proteome_windows": int(len(p_win)),
        "n_anchors_total": int(len(anchor_seqs)),
        "n_anchors_truncated": int(n_long_a),
        "n_anchor_windows": int(len(a_win)),
        "truncated_length_max": int(lengths.max()) if len(lengths) else 0,
        "truncated_length_median": float(lengths.median()) if len(lengths) else 0.0,
        "smoke_test": bool(args.limit),
    }
    with open(out / "summary.yaml", "w") as fh:
        yaml.dump(summary, fh, default_flow_style=False, sort_keys=False)

    logger.info(f"\nOutput -> {out}")
    logger.info("The _cache_*_windows.npz files are resume caches and may be "
                "deleted once the windowed_*_embeddings.npz files exist.")
    logger.info("Done OK")


if __name__ == "__main__":
    main()
