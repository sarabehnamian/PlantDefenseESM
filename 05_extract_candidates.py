#!/usr/bin/env python3
"""
05_extract_candidates.py
========================
Produce clean, publication-ready candidate gene lists from the validated
results: one table per threshold level, plus a curated high-confidence
novel candidate list.

Outputs -> results/05_extract_candidates/
    candidates_strict.csv     - Z >= 2.5
    candidates_moderate.csv   - Z >= 2.0
    candidates_lenient.csv    - Z >= 1.5
    novel_candidates.csv      - moderate-threshold, NO defense keyword
    candidate_sequences.fasta - FASTA of moderate-threshold candidates
    summary.yaml
"""

from pathlib import Path

import pandas as pd
import yaml
from Bio import SeqIO

from shared import load_config, step_dir, get_logger


def main():
    cfg = load_config()
    out = step_dir(cfg, "05_extract_candidates")
    logger = get_logger("05_extract_candidates", out)

    logger.info("=" * 65)
    logger.info("STEP 05  - Extract candidate lists")
    logger.info("=" * 65)

    # Load validated results from step 04
    val_path = (
        Path(cfg["base_output_dir"])
        / "04_validate_annotations"
        / "validated_results.csv"
    )
    df = pd.read_csv(val_path, index_col="protein_id")

    # Columns to export
    keep_cols = [
        "z_max", "best_category", "has_defense_keyword",
        "keyword_hits", "novelty", "description",
    ]
    # Add per-category z-scores
    z_cols = [c for c in df.columns if c.startswith("z_") and c != "z_max"]
    export_cols = z_cols + keep_cols

    summary = {}

    for label in ("strict", "moderate", "lenient"):
        cands = df[df[f"defense_{label}"]].copy()
        cands = cands.sort_values("z_max", ascending=False)
        out_path = out / f"candidates_{label}.csv"
        cands[export_cols].to_csv(out_path)
        summary[label] = int(len(cands))
        logger.info(f"  {label:8s}: {len(cands):>5,} candidates -> {out_path.name}")

    # ── Novel candidates (moderate, no keyword) ──────────────────────────
    novel = df[(df["defense_moderate"]) & (~df["has_defense_keyword"])].copy()
    novel = novel.sort_values("z_max", ascending=False)
    novel[export_cols].to_csv(out / "novel_candidates.csv")
    summary["novel_moderate"] = int(len(novel))
    logger.info(f"\n  Novel candidates (moderate, no keyword): {len(novel):,}")

    # ── FASTA of moderate candidates ─────────────────────────────────────
    proteome_fasta = (
        Path(cfg["base_output_dir"]) / "00_download_proteome" / "proteome.fasta"
    )
    cand_ids = set(df[df["defense_moderate"]].index)
    n_written = 0
    with open(out / "candidate_sequences.fasta", "w") as fh:
        for rec in SeqIO.parse(str(proteome_fasta), "fasta"):
            if rec.id in cand_ids:
                cat = df.loc[rec.id, "best_category"]
                z = df.loc[rec.id, "z_max"]
                fh.write(f">{rec.id}|{cat}|z={z:.2f} {rec.description}\n")
                seq = str(rec.seq)
                for i in range(0, len(seq), 80):
                    fh.write(seq[i : i + 80] + "\n")
                n_written += 1
    logger.info(f"  Wrote {n_written:,} sequences -> candidate_sequences.fasta")

    # ── Summary ──────────────────────────────────────────────────────────
    with open(out / "summary.yaml", "w") as fh:
        yaml.dump(summary, fh, default_flow_style=False)

    logger.info(f"\nOutput -> {out}")
    logger.info("Done OK")


if __name__ == "__main__":
    main()
