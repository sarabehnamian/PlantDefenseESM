#!/usr/bin/env python3
"""
r2_param_sensitivity.py
=======================
Threshold / parameter sensitivity analysis for the PlantDefenseESM revision
(Reviewer 2: parameter robustness; Reviewer 3: threshold calibration).

Runs ENTIRELY from cached step-03 outputs. No GPU, no re-embedding.

For each species it:
  1. Sweeps the percentile cutoff (the "strict/moderate/lenient" knob).
  2. Sweeps the top-N-per-category rank rule.
  3. Reports, for every setting: candidate count, % of proteome, and
     Jaccard overlap vs the default moderate set (pctl 99.0 + top-50).
  4. Reports how stable the per-category composition is.

Usage
-----
    python3 r2_param_sensitivity.py
    # or point it at custom dirs:
    python3 r2_param_sensitivity.py results/arabidopsis_thaliana results/oryza_sativa

Outputs (written next to this script):
    r2_param_sensitivity_<species>.csv   - full sweep table per species
    r2_param_sensitivity_summary.csv     - one combined tidy table
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Default species result directories (edit if yours differ)
DEFAULT_DIRS = [
    "results/arabidopsis_thaliana",
    "results/oryza_sativa",
    "results/vitis_vinifera",
]

# Sweep grids
PERCENTILES = [95.0, 96.0, 97.0, 98.0, 99.0, 99.5, 99.9]
TOP_NS = [0, 10, 25, 50, 100]

# Default operating point (the published "moderate" tier)
DEF_PCTL = 99.0
DEF_TOPN = 50


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


def topn_union(sim_df: pd.DataFrame, sim_cols, n: int) -> set:
    """Union of the top-N proteins per category by raw similarity."""
    if n <= 0:
        return set()
    ids = set()
    for col in sim_cols:
        ids.update(sim_df[col].nlargest(n).index.tolist())
    return ids


def candidate_set(z_df, sim_df, sim_cols, pctl, topn) -> set:
    """Combined set = (pctl_max >= pctl) UNION (top-N per category)."""
    pctl_ids = set(z_df.index[z_df["pctl_max"] >= pctl].tolist())
    rank_ids = topn_union(sim_df, sim_cols, topn)
    return pctl_ids | rank_ids


def analyse_species(result_dir: Path) -> pd.DataFrame:
    species = result_dir.name
    cls = result_dir / "03_classify_defense"
    z_df = pd.read_csv(cls / "zscore_matrix.csv", index_col="protein_id")
    sim_df = pd.read_csv(cls / "similarity_matrix.csv", index_col="protein_id")

    sim_cols = [c for c in sim_df.columns if c.startswith("sim_")]
    cats = [c.replace("sim_", "") for c in sim_cols]
    n_total = len(z_df)

    # best-category label per protein (recompute from sims so it is self-contained)
    best_cat = sim_df[sim_cols].idxmax(axis=1).str.replace("sim_", "", regex=False)

    # reference / default set
    ref_set = candidate_set(z_df, sim_df, sim_cols, DEF_PCTL, DEF_TOPN)
    ref_comp = best_cat.loc[list(ref_set)].value_counts()
    ref_comp = ref_comp.reindex(cats, fill_value=0)

    rows = []
    for pctl in PERCENTILES:
        for topn in TOP_NS:
            s = candidate_set(z_df, sim_df, sim_cols, pctl, topn)
            comp = best_cat.loc[list(s)].value_counts().reindex(cats, fill_value=0)
            # composition drift = L1 distance of category fractions vs reference
            ref_frac = ref_comp / max(ref_comp.sum(), 1)
            cur_frac = comp / max(comp.sum(), 1)
            comp_l1 = float(np.abs(cur_frac - ref_frac).sum())
            row = {
                "species": species,
                "percentile": pctl,
                "top_n": topn,
                "n_candidates": len(s),
                "pct_proteome": round(100 * len(s) / n_total, 3),
                "jaccard_vs_default": round(jaccard(s, ref_set), 4),
                "frac_default_retained": round(
                    len(s & ref_set) / max(len(ref_set), 1), 4),
                "category_comp_L1_drift": round(comp_l1, 4),
                "is_default": (pctl == DEF_PCTL and topn == DEF_TOPN),
            }
            rows.append(row)

    df = pd.DataFrame(rows)
    out = Path(__file__).resolve().parent / f"r2_param_sensitivity_{species}.csv"
    df.to_csv(out, index=False)

    print(f"\n===== {species}  (proteome N={n_total:,}) =====")
    print(f"Default set (pctl {DEF_PCTL}, top-{DEF_TOPN}): "
          f"{len(ref_set):,} candidates")
    # Compact view: vary percentile at the default top-N
    print(f"\n  Percentile sweep at top-{DEF_TOPN}:")
    sub = df[df.top_n == DEF_TOPN][
        ["percentile", "n_candidates", "pct_proteome", "jaccard_vs_default"]]
    print(sub.to_string(index=False))
    # Compact view: vary top-N at the default percentile
    print(f"\n  Top-N sweep at percentile {DEF_PCTL}:")
    sub2 = df[df.percentile == DEF_PCTL][
        ["top_n", "n_candidates", "pct_proteome", "jaccard_vs_default",
         "category_comp_L1_drift"]]
    print(sub2.to_string(index=False))
    return df


def main():
    dirs = sys.argv[1:] if len(sys.argv) > 1 else DEFAULT_DIRS
    all_dfs = []
    for d in dirs:
        p = Path(d)
        if not (p / "03_classify_defense" / "zscore_matrix.csv").exists():
            print(f"[skip] {d}: no 03_classify_defense/zscore_matrix.csv")
            continue
        all_dfs.append(analyse_species(p))

    if all_dfs:
        combined = pd.concat(all_dfs, ignore_index=True)
        out = Path(__file__).resolve().parent / "r2_param_sensitivity_summary.csv"
        combined.to_csv(out, index=False)
        print(f"\nWrote combined summary -> {out}")
        print("Per-species tables -> r2_param_sensitivity_<species>.csv")


if __name__ == "__main__":
    main()
