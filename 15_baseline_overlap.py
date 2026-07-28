#!/usr/bin/env python3
"""
15_baseline_overlap.py
======================
Complete overlap breakdown of the curated defense proteins recovered by each
method.  (Reviewer 1 Comment 6.)

12_benchmark_baselines.py already reports a partial breakdown, but it collapses
the two non-embedding methods into a single "standard" class (domain OR
alignment), so only the ESM-2-only count (the 70 quoted in the manuscript) can
be read off it.  The Reviewer asks for the full picture: proteins recovered by
all methods, by each method alone, the pairwise overlaps, and the proteins
missed by every method.

Three methods, defined EXACTLY as in 12_benchmark_baselines.py so that the
numbers are consistent with the metrics already reported:

  ESM-2       - the moderate-tier candidate set (defense_moderate flag)
  Alignment   - best local-alignment score to any of the 33 anchors
                (NCBI blastp bitscore if BLAST+ is on PATH, otherwise the
                Biopython Smith-Waterman / BLOSUM62 fallback), candidate set
                = top-N by score with N matched to the ESM-2 moderate count
  Domain      - proteins carrying any diagnostic defense InterPro domain
                (NB-ARC / PR / LRR+kinase)

Ground truth is the same curated GO defense set used in 09_benchmark_curated.py
(GO:0006952 defense response + GO:0002376 immune system process, mapped to
RefSeq accessions).

The three methods give 2^3 = 8 mutually exclusive cells that sum to the curated
positive set, which is the compact table the Reviewer asks for.  Inclusive
pairwise overlaps are reported separately.

The alignment baseline is the slow part.  Its per-protein scores are cached, so
the first run costs the same as 12_benchmark_baselines.py and every later run
is instant.  Use --refresh to force recomputation.

Run from the project root:
    python 15_baseline_overlap.py

Outputs -> results/15_baseline_overlap/
    overlap_cells_<species>.csv        - the 8 mutually exclusive cells
    overlap_pairwise_<species>.csv     - totals and inclusive pairwise overlaps
    overlap_table_<species>.csv        - manuscript-ready compact table
    positive_membership_<species>.csv  - one row per curated positive, with the
                                         three method flags and the description
    summary.txt                        - human-readable summary (also printed)
"""

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

RESULTS = Path("results")
OUTDIR = RESULTS / "15_baseline_overlap"

# the ESM-2-only count currently stated in the manuscript; checked, not assumed
MANUSCRIPT_ESM_ONLY = 70


def load_step12():
    """Import 12_benchmark_baselines.py (module name starts with a digit)."""
    path = Path(__file__).resolve().parent / "12_benchmark_baselines.py"
    if not path.exists():
        raise FileNotFoundError(
            f"Cannot find {path} (run this from the project root)")
    spec = importlib.util.spec_from_file_location("_step12_baselines", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def find_fasta(folder: str, patterns) -> Path:
    """Same lookup logic as 12_benchmark_baselines.py."""
    for pat in patterns:
        hits = sorted(RESULTS.glob(pat))
        if hits:
            return hits[0]
    return Path("___missing___")


def homology_scores(step12, folder: str, cache_dir: Path, refresh: bool, log):
    """{versionless_protein_id: best score to any anchor}, cached to CSV."""
    import shutil

    cache_file = cache_dir / f"homology_scores_{folder}.csv"
    if cache_file.exists() and not refresh:
        log(f"  [ALIGN] using cached scores -> {cache_file.name}")
        cached = pd.read_csv(cache_file)
        return dict(zip(cached["pid"], cached["score"])), "anchor alignment (cached)"

    prot_fa = find_fasta(folder, [f"{folder}/00_download_proteome/proteome.fasta",
                                  f"{folder}/**/proteome.fasta",
                                  f"{folder}/**/*.fasta"])
    anchor_fa = find_fasta(folder, [f"{folder}/01_fetch_anchors/anchors.fasta",
                                    "01_fetch_anchors/anchors.fasta",
                                    "**/anchors.fasta"])
    log(f"  proteome FASTA: {prot_fa}")
    log(f"  anchor   FASTA: {anchor_fa}")

    if shutil.which("blastp") and shutil.which("makeblastdb"):
        best = step12.blast_best_bitscore(prot_fa, anchor_fa, log)
        label = "BLAST-to-anchor"
    else:
        best = step12.alignment_best_score(prot_fa, anchor_fa, log)
        label = "alignment-to-anchor (SW)"

    if best is None:
        return None, label

    pd.DataFrame({"pid": list(best.keys()), "score": list(best.values())}
                 ).to_csv(cache_file, index=False)
    log(f"  [ALIGN] cached scores -> {cache_file.name}")
    return best, label


def cell_table(pos_df):
    """8 mutually exclusive cells over the three method flags."""
    e = pos_df["ESM2"].values
    a = pos_df["alignment"].values
    d = pos_df["domain"].values
    n = len(pos_df)

    cells = [
        ("Recovered by all three methods",        e & a & d),
        ("ESM-2 and alignment only",              e & a & ~d),
        ("ESM-2 and domain only",                 e & ~a & d),
        ("Alignment and domain only",             ~e & a & d),
        ("ESM-2 only",                            e & ~a & ~d),
        ("Alignment only",                        ~e & a & ~d),
        ("Domain only",                           ~e & ~a & d),
        ("Missed by all three methods",           ~e & ~a & ~d),
    ]
    rows = [{"cell": name,
             "n": int(m.sum()),
             "pct_of_curated": round(100 * int(m.sum()) / n, 1) if n else 0.0}
            for name, m in cells]
    return pd.DataFrame(rows)


def pairwise_table(pos_df):
    """Inclusive totals and pairwise overlaps (these deliberately overlap)."""
    e = pos_df["ESM2"].values
    a = pos_df["alignment"].values
    d = pos_df["domain"].values
    n = len(pos_df)
    rows = [
        ("Curated defense proteins in proteome", n),
        ("Recovered by ESM-2 (total)", int(e.sum())),
        ("Recovered by alignment (total)", int(a.sum())),
        ("Recovered by domain retrieval (total)", int(d.sum())),
        ("ESM-2 and alignment (inclusive)", int((e & a).sum())),
        ("ESM-2 and domain (inclusive)", int((e & d).sum())),
        ("Alignment and domain (inclusive)", int((a & d).sum())),
        ("Recovered by at least one method", int((e | a | d).sum())),
    ]
    return pd.DataFrame(rows, columns=["quantity", "n"])


def manuscript_table(cells, pairs, n_pos):
    """The compact table for the manuscript: the 8 cells plus method totals."""
    out = cells.rename(columns={"cell": "Category", "n": "Curated proteins (n)",
                                "pct_of_curated": "% of curated set"})
    totals = pairs[pairs["quantity"].str.contains("total|at least one")].copy()
    totals["% of curated set"] = (100 * totals["n"] / n_pos).round(1)
    totals = totals.rename(columns={"quantity": "Category",
                                    "n": "Curated proteins (n)"})
    blank = pd.DataFrame([{"Category": "", "Curated proteins (n)": "",
                           "% of curated set": ""}])
    return pd.concat([out, blank, totals], ignore_index=True)


def main():
    ap = argparse.ArgumentParser(
        description="Full method-overlap breakdown among curated defense proteins")
    ap.add_argument("--refresh", action="store_true",
                    help="recompute the alignment scores instead of using the cache")
    ap.add_argument("--species", default=None,
                    help="restrict to one results folder, e.g. arabidopsis_thaliana")
    args = ap.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    cache = OUTDIR / "_cache"
    cache.mkdir(exist_ok=True)

    lines = []
    def log(s=""):
        print(s)
        lines.append(s)

    step12 = load_step12()

    log("=" * 72)
    log("Method overlap among curated defense proteins  (Reviewer 1, comment 6)")
    log("ESM-2 (moderate)  |  alignment-to-anchor (matched-N)  |  InterPro domain")
    log("=" * 72)

    species = step12.SPECIES
    if args.species:
        species = {k: v for k, v in species.items() if k == args.species}
        if not species:
            sys.exit(f"ERROR: unknown species folder '{args.species}'")

    for folder, (display, taxid) in species.items():
        val = RESULTS / folder / "04_validate_annotations" / "validated_results.csv"
        if not val.exists():
            log(f"\n[SKIP] {val} not found")
            continue

        log(f"\n### {display}  (taxid {taxid})")
        df = pd.read_csv(val)
        if "protein_id" not in df.columns:
            df = df.rename(columns={df.columns[0]: "protein_id"})
        df["pid"] = df["protein_id"].map(step12.strip_ver)

        curated = step12.fetch_curated_refseq(taxid, cache, log)
        y_true = df["pid"].isin(curated).values
        n_pos = int(y_true.sum())
        log(f"  curated positives mapped: {n_pos:,} / {len(df):,}")

        if n_pos < step12.MIN_CURATED:
            log(f"  [SKIP] fewer than {step12.MIN_CURATED} curated positives; "
                f"an overlap breakdown would not be interpretable here.")
            continue

        # --- method 1: ESM-2 moderate tier ---
        esm_mask = df["defense_moderate"].values.astype(bool)
        n_esm = int(esm_mask.sum())

        # --- method 2: alignment to anchors, top-N matched to the ESM-2 count ---
        best, hlabel = homology_scores(step12, folder, cache, args.refresh, log)
        if best is None:
            log("  [ABORT] no alignment baseline available; cannot build the "
                "three-way breakdown for this species.")
            continue
        score = df["pid"].map(lambda p: best.get(p, 0.0)).values
        order = np.argsort(-score, kind="stable")
        align_mask = np.zeros(len(df), bool)
        align_mask[order[:n_esm]] = True

        # --- method 3: InterPro defense domains ---
        domain_ids = step12.fetch_domain_refseq(taxid, cache, log)
        domain_mask = df["pid"].isin(domain_ids).values

        log(f"  method sets: ESM-2 N={n_esm:,}  {hlabel} N={int(align_mask.sum()):,}  "
            f"domain N={int(domain_mask.sum()):,}")

        # --- restrict to the curated positives and break down ---
        pos_df = pd.DataFrame({
            "protein_id": df.loc[y_true, "protein_id"].values,
            "pid": df.loc[y_true, "pid"].values,
            "ESM2": esm_mask[y_true],
            "alignment": align_mask[y_true],
            "domain": domain_mask[y_true],
        })
        if "z_max" in df.columns:
            pos_df["z_max"] = df.loc[y_true, "z_max"].values.round(3)
        if "best_category" in df.columns:
            pos_df["esm_category"] = df.loc[y_true, "best_category"].values
        if "has_defense_keyword" in df.columns:
            pos_df["refseq_keyword"] = df.loc[y_true, "has_defense_keyword"].values
        if "description" in df.columns:
            pos_df["description"] = df.loc[y_true, "description"].values

        cells = cell_table(pos_df)
        pairs = pairwise_table(pos_df)
        table = manuscript_table(cells, pairs, n_pos)

        cells.to_csv(OUTDIR / f"overlap_cells_{folder}.csv", index=False)
        pairs.to_csv(OUTDIR / f"overlap_pairwise_{folder}.csv", index=False)
        table.to_csv(OUTDIR / f"overlap_table_{folder}.csv", index=False)
        pos_df.to_csv(OUTDIR / f"positive_membership_{folder}.csv", index=False)

        log("")
        log(cells.to_string(index=False))
        log("")
        log(pairs.to_string(index=False))

        # --- consistency checks ---
        total = int(cells["n"].sum())
        log("")
        log(f"  check: cells sum to {total:,} (curated positives = {n_pos:,}) -> "
            f"{'OK' if total == n_pos else 'MISMATCH'}")

        esm_only = int(cells.loc[cells["cell"] == "ESM-2 only", "n"].iloc[0])
        if folder == "arabidopsis_thaliana":
            log(f"  check: ESM-2 only = {esm_only} vs {MANUSCRIPT_ESM_ONLY} "
                f"stated in the manuscript -> "
                f"{'OK' if esm_only == MANUSCRIPT_ESM_ONLY else 'DIFFERS - update the manuscript'}")

        old = RESULTS / "12_benchmark_baselines" / "recovered_breakdown.csv"
        if old.exists():
            prev = pd.read_csv(old)
            row = prev[prev["species"] == display]
            if not row.empty and "ESM2_only_vs_standard" in row.columns:
                prev_only = int(row["ESM2_only_vs_standard"].iloc[0])
                log(f"  check: ESM-2 only = {esm_only} vs {prev_only} from "
                    f"12_benchmark_baselines.py -> "
                    f"{'OK' if esm_only == prev_only else 'MISMATCH - methods drifted'}")

    (OUTDIR / "summary.txt").write_text("\n".join(lines), encoding="utf-8")
    log(f"\nOutputs -> {OUTDIR}")
    log("Done.")


if __name__ == "__main__":
    if not RESULTS.exists():
        sys.exit("ERROR: run from the project root (no 'results/' folder here).")
    main()
