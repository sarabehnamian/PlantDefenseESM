#!/usr/bin/env python3
"""
r3_word_boundary_check.py
=========================
Reviewer 3, minor comment: "clarify whether keyword matching used word
boundaries. Substring matching could introduce avoidable false positives."

The pipeline (04_validate_annotations.py) uses case-insensitive SUBSTRING
matching: `kw.lower() in description.lower()`. This script quantifies whether
that inflates the enrichment, by recomputing the moderate-tier enrichment with
WHOLE-WORD (word-boundary) matching and comparing.

It imports VALIDATION_KEYWORDS from shared.py, so the keyword list is identical
to the pipeline. As a faithfulness check, the substring fold it computes should
reproduce the published all-47 moderate folds (4.22 / 3.35 / 3.87).

Run from the project root (the folder containing results/ and shared.py):
    python r3_word_boundary_check.py
"""

from pathlib import Path
import re
import sys
import pandas as pd

from shared import VALIDATION_KEYWORDS as K

SPECIES = {
    "arabidopsis_thaliana": "A. thaliana",
    "vitis_vinifera":       "V. vinifera",
    "oryza_sativa":         "O. sativa",
}
RESULTS = Path("results")

# word-boundary matcher: keyword not flanked by another letter/digit
# (handles hyphenated tokens like "NB-ARC", "PR-", "R-gene" correctly, and
#  prevents short tokens like "STS"/"PAL"/"TIR" matching inside longer words)
PATS = [re.compile(r'(?<![a-z0-9])' + re.escape(k.lower()) + r'(?![a-z0-9])') for k in K]
KLOW = [k.lower() for k in K]


def has_substring(desc: str) -> bool:
    return any(k in desc for k in KLOW)


def has_wordboundary(desc: str) -> bool:
    return any(p.search(desc) for p in PATS)


def main():
    rows = []
    print(f"{'species':12s} {'sub fold':>9s} {'wb fold':>9s} "
          f"{'bg sub':>7s} {'bg wb':>7s} {'changed':>9s}")
    for folder, disp in SPECIES.items():
        val = RESULTS / folder / "04_validate_annotations" / "validated_results.csv"
        if not val.exists():
            print(f"  [SKIP] {val} not found")
            continue
        df = pd.read_csv(val)
        desc = df["description"].fillna("").astype(str).str.lower()
        cand = df["defense_moderate"].astype(bool).values

        hs = desc.apply(has_substring)
        hw = desc.apply(has_wordboundary)
        bg_s, bg_w = hs.mean(), hw.mean()
        fold_s = hs[cand].mean() / bg_s
        fold_w = hw[cand].mean() / bg_w
        changed = int((hs != hw).sum())
        pct = 100 * changed / len(df)

        print(f"{disp:12s} {fold_s:>9.2f} {fold_w:>9.2f} "
              f"{bg_s:>7.3f} {bg_w:>7.3f} {changed:>6d} ({pct:.2f}%)")
        rows.append({
            "species": disp,
            "n_proteome": len(df),
            "fold_substring_moderate": round(fold_s, 2),
            "fold_wordboundary_moderate": round(fold_w, 2),
            "background_substring": round(bg_s, 4),
            "background_wordboundary": round(bg_w, 4),
            "n_proteins_status_changed": changed,
            "pct_proteins_changed": round(pct, 3),
        })

    if rows:
        pd.DataFrame(rows).to_csv("r3_word_boundary_check.csv", index=False)
        print("\nWrote r3_word_boundary_check.csv")
        print("Sanity check: 'sub fold' should read 4.22 / 3.35 / 3.87 "
              "(matches the published all-47 moderate folds).")


if __name__ == "__main__":
    if not RESULTS.exists():
        sys.exit("ERROR: run from the project root (no 'results/' folder here).")
    main()
