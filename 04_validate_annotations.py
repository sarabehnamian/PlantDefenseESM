#!/usr/bin/env python3
"""
04_validate_annotations.py
==========================
Independent validation: check how well ESM-2 candidates overlap with
defense-related annotations from the genome project (RefSeq / InterPro).

This is NOT circular: annotations come from domain databases, while our
predictions come purely from sequence embeddings.

Outputs -> results/04_validate_annotations/
    validated_results.csv     - full table with keyword hits, novelty class
    enrichment_tests.csv      - Fisher's exact test at each threshold
    novelty_breakdown.csv     - known / novel / keyword-only counts
    summary.yaml
"""

from pathlib import Path

import pandas as pd
import yaml
from scipy.stats import fisher_exact

from shared import load_config, step_dir, get_logger, VALIDATION_KEYWORDS


def load_descriptions(proteome_stats: Path) -> dict:
    """protein_id -> description string."""
    df = pd.read_csv(proteome_stats)
    return dict(zip(df["protein_id"], df["description"].fillna("")))


def keyword_scan(description: str, keywords: list) -> list:
    """Return which keywords appear in a protein description."""
    desc_lower = description.lower()
    return [kw for kw in keywords if kw.lower() in desc_lower]


def main():
    cfg = load_config()
    out = step_dir(cfg, "04_validate_annotations")
    logger = get_logger("04_validate_annotations", out)

    logger.info("=" * 65)
    logger.info("STEP 04  - Validate against keyword annotations")
    logger.info("=" * 65)

    # Load Z-score results from step 03
    z_path = (
        Path(cfg["base_output_dir"]) / "03_classify_defense" / "zscore_matrix.csv"
    )
    z_df = pd.read_csv(z_path, index_col="protein_id")

    # Load protein descriptions from step 00
    desc_path = (
        Path(cfg["base_output_dir"]) / "00_download_proteome" / "proteome_stats.csv"
    )
    descs = load_descriptions(desc_path)

    # ── Keyword scan ─────────────────────────────────────────────────────
    logger.info(f"Scanning {len(z_df):,} proteins against "
                f"{len(VALIDATION_KEYWORDS)} defense keywords ...")

    z_df["description"] = z_df.index.map(lambda x: descs.get(x, ""))
    z_df["keyword_hits"] = z_df["description"].apply(
        lambda d: "|".join(keyword_scan(d, VALIDATION_KEYWORDS))
    )
    z_df["has_defense_keyword"] = z_df["keyword_hits"].str.len() > 0

    n_total = len(z_df)
    n_kw = int(z_df["has_defense_keyword"].sum())
    bg_rate = n_kw / n_total
    logger.info(f"Background: {n_kw:,} / {n_total:,} ({bg_rate:.1%}) have "
                f"defense keywords")

    # ── Enrichment test per threshold ────────────────────────────────────
    enrichment_rows = []
    for label in ("strict", "moderate", "lenient"):
        col = f"defense_{label}"
        cands = z_df[z_df[col]]
        n_cand = len(cands)
        n_cand_kw = int(cands["has_defense_keyword"].sum())
        cand_rate = n_cand_kw / n_cand if n_cand else 0
        enrichment = cand_rate / bg_rate if bg_rate > 0 else float("inf")

        # 2x2 contingency
        a = n_cand_kw
        b = n_cand - n_cand_kw
        c = n_kw - n_cand_kw
        d = n_total - n_cand - c
        odds, pval = fisher_exact([[a, b], [c, d]], alternative="greater")

        enrichment_rows.append({
            "threshold": label,
            "percentile_cutoff": cfg.get(f"percentile_{label}", "N/A"),
            "n_candidates": n_cand,
            "n_with_keyword": n_cand_kw,
            "rate": round(cand_rate, 4),
            "background_rate": round(bg_rate, 4),
            "fold_enrichment": round(enrichment, 2),
            "odds_ratio": round(odds, 2),
            "fisher_p": f"{pval:.2e}",
        })
        logger.info(
            f"\n  {label} (P >= {cfg.get(f'percentile_{label}', 'N/A')} | top-N):"
            f"\n    Candidates      : {n_cand:,}"
            f"\n    With keywords   : {n_cand_kw:,} ({cand_rate:.1%})"
            f"\n    Fold enrichment : {enrichment:.2f}x"
            f"\n    Fisher p-value  : {pval:.2e}"
        )

    enrich_df = pd.DataFrame(enrichment_rows)
    enrich_df.to_csv(out / "enrichment_tests.csv", index=False)

    # ── Novelty classification ───────────────────────────────────────────
    logger.info("\nNovelty classification (moderate threshold):")
    z_df["novelty"] = "non_candidate"
    mask_def = z_df["defense_moderate"]
    mask_kw = z_df["has_defense_keyword"]

    z_df.loc[mask_def & mask_kw, "novelty"] = "known_defense"
    z_df.loc[mask_def & ~mask_kw, "novelty"] = "novel_candidate"
    z_df.loc[~mask_def & mask_kw, "novelty"] = "keyword_only"

    novelty = z_df["novelty"].value_counts()
    for cat, n in novelty.items():
        logger.info(f"  {cat:20s}  {n:>6,}")

    novelty_df = novelty.reset_index()
    novelty_df.columns = ["class", "count"]
    novelty_df.to_csv(out / "novelty_breakdown.csv", index=False)

    # ── Save full validated table ────────────────────────────────────────
    z_df.to_csv(out / "validated_results.csv")

    # ── Summary ──────────────────────────────────────────────────────────
    summary = {
        "n_proteins": n_total,
        "n_with_defense_keyword": n_kw,
        "enrichment": enrichment_rows,
        "novelty_moderate": novelty.to_dict(),
    }
    with open(out / "summary.yaml", "w") as fh:
        yaml.dump(summary, fh, default_flow_style=False)

    logger.info(f"\nOutput -> {out}")
    logger.info("Done OK")


if __name__ == "__main__":
    main()
