#!/usr/bin/env python3
"""
r3_solve_keyword_split.py
=========================
Solve for the exact 8 "broad/generic" keywords. The calculation logic is
already confirmed: with NO terms removed it reproduces the manuscript's
all-keyword folds (4.22 / 3.35 / 3.87) exactly. So the only unknown is which
8 keywords were placed in the broad group. This script finds the subset that
reproduces all three published defense-specific folds simultaneously:

        A. thaliana 6.89   V. vinifera 5.97   O. sativa 3.77

It runs from cached results only (no GPU). It prunes to the keywords that
actually move the enrichment, then brute-forces 8-term subsets of that pool.

Run from the project root:
    python r3_solve_keyword_split.py
"""

from pathlib import Path
import itertools
import sys
import numpy as np
import pandas as pd

from shared import VALIDATION_KEYWORDS as K

SPECIES = {
    "arabidopsis_thaliana": "A. thaliana",
    "vitis_vinifera":       "V. vinifera",
    "oryza_sativa":         "O. sativa",
}
TARGET = {"A. thaliana": 6.89, "V. vinifera": 5.97, "O. sativa": 3.77}
RESULTS = Path("results")
TOL = 0.06

# ---- load each species: keyword-hit matrix, candidate mask, total hits ------
data = {}
for folder, disp in SPECIES.items():
    val = RESULTS / folder / "04_validate_annotations" / "validated_results.csv"
    if not val.exists():
        sys.exit(f"Missing {val}")
    df = pd.read_csv(val)
    desc = pd.Series(df["description"].fillna("").astype(str).str.lower().values)
    M = np.zeros((len(df), len(K)), dtype=bool)
    for j, kw in enumerate(K):
        M[:, j] = desc.str.contains(kw.lower(), regex=False).values
    cand = df["defense_moderate"].astype(bool).values
    data[disp] = (M, cand, M.sum(1))


def fold(disp, broad_idx):
    """Moderate-tier fold using only the SPECIFIC keywords (= all minus broad)."""
    M, cand, total = data[disp]
    if broad_idx:
        spec = (total - M[:, list(broad_idx)].sum(1)) >= 1
    else:
        spec = total >= 1
    bg = spec.mean()
    return spec[cand].mean() / bg if bg else float("nan")


def main():
    print("Sanity (broad = none) vs manuscript all-keyword folds 4.22/3.35/3.87:")
    for disp in SPECIES.values():
        print(f"  {disp:12s} {fold(disp, []):.2f}")

    # pool = keywords whose single removal moves the fold in ANY species
    pool = []
    for j in range(len(K)):
        moved = any(abs(fold(d, [j]) - fold(d, [])) >= 0.008 for d in SPECIES.values())
        if moved:
            pool.append(j)
    print(f"\nActive pool ({len(pool)}): {[K[j] for j in pool]}")

    # brute-force 8-term broad subsets of the pool
    exact = []
    scored = []
    for combo in itertools.combinations(pool, 8):
        folds = {d: fold(d, combo) for d in SPECIES.values()}
        err = sum(abs(folds[d] - TARGET[d]) for d in SPECIES.values())
        scored.append((err, combo, folds))
        if all(abs(folds[d] - TARGET[d]) <= TOL for d in SPECIES.values()):
            exact.append((combo, folds))

    print(f"\nExact matches (all 3 within {TOL}): {len(exact)}")
    for combo, folds in exact[:25]:
        print("  BROAD =", [K[j] for j in combo],
              {d: round(v, 2) for d, v in folds.items()})

    if not exact:
        print("\nNo exact 8-term match in pool. Closest 12 subsets:")
        for err, combo, folds in sorted(scored)[:12]:
            print(f"  err={err:.2f}  BROAD =", [K[j] for j in combo],
                  {d: round(v, 2) for d, v in folds.items()})

    # write the table for whichever broad set is chosen (best match)
    best = (exact[0][0] if exact else sorted(scored)[0][1])
    broad = [K[j] for j in best]
    specific = [k for k in K if k not in broad]
    rows = ([{"keyword": k, "group": "defense-specific"} for k in specific] +
            [{"keyword": k, "group": "broad/generic"} for k in broad])
    pd.DataFrame(rows).to_csv("keyword_groups.csv", index=False)
    print("\nWrote keyword_groups.csv for the best-matching split.")


if __name__ == "__main__":
    if not RESULTS.exists():
        sys.exit("ERROR: run from the project root (no 'results/' folder here).")
    main()
