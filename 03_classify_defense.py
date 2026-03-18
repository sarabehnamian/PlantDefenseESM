#!/usr/bin/env python3
"""
03_classify_defense.py
======================
For each protein in the proteome:
  1. Compute cosine similarity to the centroid of each defense category.
  2. Convert similarities to Z-scores (null = whole-proteome distribution).
  3. Compute per-category percentile ranks.
  4. Assign defense candidate status using TWO complementary methods:
     (a) Per-category percentile thresholds (ensures coverage across categories)
     (b) Per-category top-N ranking (guarantees minimum representation)

The union of both methods defines the final candidate set.

Outputs -> results/03_classify_defense/
    similarity_matrix.csv     - raw cosine similarities (N x 6 categories)
    zscore_matrix.csv         - Z-scores + percentile ranks + candidate flags
    category_centroids.npz    - centroid vectors for downstream use
    top_per_category.csv      - top candidates within each category
    summary.yaml              - counts at each threshold
"""

from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats as sp_stats

from shared import load_config, step_dir, get_logger


def load_embeddings(npz_path: Path):
    """Load .npz -> (embeddings, protein_ids)."""
    data = np.load(npz_path, allow_pickle=True)
    return data["embeddings"], list(data["protein_ids"])


def compute_centroids(anchor_emb, anchor_ids, cfg, logger):
    """
    Group anchor embeddings by category and compute L2-normalized centroid.
    """
    base = Path(cfg["base_output_dir"]) / "01_fetch_anchors"
    meta = pd.read_csv(base / "anchors_metadata.csv")
    id_to_cat = dict(zip(meta["uniprot_id"], meta["category"]))

    id_to_idx = {pid: i for i, pid in enumerate(anchor_ids)}

    centroids = {}
    for category in sorted(set(id_to_cat.values())):
        cat_ids = [uid for uid, cat in id_to_cat.items()
                   if cat == category and uid in id_to_idx]
        if not cat_ids:
            logger.warning(f"No anchors for {category}")
            continue
        idxs = [id_to_idx[uid] for uid in cat_ids]
        c = anchor_emb[idxs].mean(axis=0)
        c = c / np.linalg.norm(c)
        centroids[category] = c
        logger.info(f"  {category:20s}  {len(cat_ids)} anchors")

    return centroids


def cosine_similarity_matrix(proteome_emb, centroids, categories):
    """Return (N, C) matrix of cosine similarities."""
    norms = np.linalg.norm(proteome_emb, axis=1, keepdims=True)
    norms = np.where(norms == 0, 1.0, norms)
    normed = proteome_emb / norms
    C = np.stack([centroids[cat] for cat in categories])
    return normed @ C.T


def main():
    cfg = load_config()
    out = step_dir(cfg, "03_classify_defense")
    logger = get_logger("03_classify_defense", out)

    # Config parameters with defaults
    top_n_per_cat = cfg.get("top_n_per_category", 50)
    pctl_strict = cfg.get("percentile_strict", 99.5)
    pctl_moderate = cfg.get("percentile_moderate", 99.0)
    pctl_lenient = cfg.get("percentile_lenient", 97.0)

    logger.info("=" * 65)
    logger.info("STEP 03  - Classify defense candidates")
    logger.info("=" * 65)

    # Load embeddings from step 02
    emb_dir = Path(cfg["base_output_dir"]) / "02_embed_proteins"
    anchor_emb, anchor_ids = load_embeddings(emb_dir / "anchor_embeddings.npz")
    proteome_emb, proteome_ids = load_embeddings(
        emb_dir / "proteome_embeddings.npz"
    )
    logger.info(f"Anchors  : {anchor_emb.shape}")
    logger.info(f"Proteome : {proteome_emb.shape}")

    # -- Centroids ---------------------------------------------------------
    logger.info("\nCategory centroids:")
    centroids = compute_centroids(anchor_emb, anchor_ids, cfg, logger)
    categories = sorted(centroids.keys())

    np.savez(
        out / "category_centroids.npz",
        **{cat: centroids[cat] for cat in categories},
        categories=np.array(categories),
    )

    # -- Cosine similarity -------------------------------------------------
    logger.info("\nComputing cosine similarities ...")
    sim = cosine_similarity_matrix(proteome_emb, centroids, categories)

    sim_df = pd.DataFrame(
        sim, index=proteome_ids,
        columns=[f"sim_{cat}" for cat in categories],
    )
    sim_df.index.name = "protein_id"
    sim_df.to_csv(out / "similarity_matrix.csv")
    logger.info(f"Similarity matrix: {sim_df.shape}")

    # -- Z-scores (global) -------------------------------------------------
    logger.info("Computing Z-scores ...")
    z_df = pd.DataFrame(index=proteome_ids)
    z_df.index.name = "protein_id"

    for cat in categories:
        col = f"sim_{cat}"
        z_df[f"z_{cat}"] = sp_stats.zscore(sim_df[col].values)

    z_cols = [f"z_{cat}" for cat in categories]
    z_df["z_max"] = z_df[z_cols].max(axis=1)
    z_df["best_category"] = (
        z_df[z_cols].idxmax(axis=1).str.replace("z_", "", regex=False)
    )

    # -- Per-category percentile ranks -------------------------------------
    logger.info("Computing per-category percentile ranks ...")
    for cat in categories:
        sim_col = f"sim_{cat}"
        z_df[f"pctl_{cat}"] = sp_stats.rankdata(
            sim_df[sim_col].values
        ) / len(proteome_ids) * 100

    pctl_cols = [f"pctl_{cat}" for cat in categories]
    z_df["pctl_max"] = z_df[pctl_cols].max(axis=1)
    z_df["pctl_best_category"] = (
        z_df[pctl_cols].idxmax(axis=1).str.replace("pctl_", "", regex=False)
    )

    # -- METHOD A: Per-category percentile thresholds ----------------------
    for label, pctl in [
        ("strict", pctl_strict),
        ("moderate", pctl_moderate),
        ("lenient", pctl_lenient),
    ]:
        z_df[f"defense_pctl_{label}"] = z_df["pctl_max"] >= pctl

    # -- METHOD B: Top-N per category (rank-based) -------------------------
    logger.info(f"Selecting top {top_n_per_cat} per category ...")
    top_per_cat_rows = []
    top_ids_union = set()

    for cat in categories:
        sim_col = f"sim_{cat}"
        cat_sorted = sim_df[sim_col].nlargest(top_n_per_cat)
        for pid, score in cat_sorted.items():
            top_per_cat_rows.append({
                "protein_id": pid,
                "category": cat,
                "similarity": round(float(score), 6),
                "z_score": round(float(z_df.loc[pid, f"z_{cat}"]), 4),
                "percentile": round(float(z_df.loc[pid, f"pctl_{cat}"]), 2),
            })
            top_ids_union.add(pid)

    top_cat_df = pd.DataFrame(top_per_cat_rows)
    top_cat_df.to_csv(out / "top_per_category.csv", index=False)

    z_df["defense_topN"] = z_df.index.isin(top_ids_union)

    # -- COMBINED: union of percentile + topN = final candidate set --------
    for label in ("strict", "moderate", "lenient"):
        z_df[f"defense_{label}"] = (
            z_df[f"defense_pctl_{label}"] | z_df["defense_topN"]
        )

    z_df.to_csv(out / "zscore_matrix.csv")

    # -- Summary -----------------------------------------------------------
    n_total = len(proteome_ids)
    summary = {"n_proteins": n_total, "top_n_per_category": top_n_per_cat}

    logger.info("\n--- Per-category percentile thresholds ---")
    summary["pctl_thresholds"] = {}
    for label, pctl in [("strict", pctl_strict), ("moderate", pctl_moderate),
                         ("lenient", pctl_lenient)]:
        n = int(z_df[f"defense_pctl_{label}"].sum())
        summary["pctl_thresholds"][label] = {
            "percentile": float(pctl),
            "n_candidates": n,
            "pct": round(100 * n / n_total, 2),
        }
        logger.info(f"  P >= {pctl:5.1f} ({label:8s}): {n:>5,} ({100*n/n_total:.2f}%)")

    logger.info(f"\n--- Top-{top_n_per_cat} per category ---")
    n_topN = int(z_df["defense_topN"].sum())
    summary["topN_candidates"] = n_topN
    logger.info(f"  Unique proteins in top-{top_n_per_cat}: {n_topN:,}")

    logger.info("\n--- Combined candidate set (percentile | top-N) ---")
    summary["combined"] = {}
    for label in ("strict", "moderate", "lenient"):
        n = int(z_df[f"defense_{label}"].sum())
        summary["combined"][label] = {
            "n_candidates": n,
            "pct": round(100 * n / n_total, 2),
        }
        logger.info(f"  {label:8s}: {n:>5,} ({100*n/n_total:.2f}%)")

    # Category breakdown for combined moderate
    mod = z_df[z_df["defense_moderate"]]
    cat_breakdown = mod["best_category"].value_counts().to_dict()
    summary["category_breakdown_moderate"] = {
        k: int(v) for k, v in cat_breakdown.items()
    }
    logger.info("\nCandidates by best category (combined moderate):")
    for cat, n in sorted(cat_breakdown.items(), key=lambda x: -x[1]):
        logger.info(f"  {cat:20s}  {n:>5,}")

    # Per-category top-N breakdown
    logger.info(f"\nTop-{top_n_per_cat} per category breakdown:")
    for cat in categories:
        cat_data = top_cat_df[top_cat_df["category"] == cat]
        n_cat = len(cat_data)
        min_z = cat_data["z_score"].min()
        max_z = cat_data["z_score"].max()
        logger.info(f"  {cat:20s}  {n_cat:>3} proteins  Z=[{min_z:.2f}, {max_z:.2f}]")

    with open(out / "summary.yaml", "w") as fh:
        yaml.dump(summary, fh, default_flow_style=False)

    logger.info(f"\nOutput -> {out}")
    logger.info("Done OK")


if __name__ == "__main__":
    main()
