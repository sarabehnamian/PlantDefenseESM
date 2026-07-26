#!/usr/bin/env python3
"""
08_make_supplementary_tables.py
===============================
Builds clean, machine-readable supplementary tables of predicted defense
candidates, one per species, in response to Reviewer 1 Comment 2.

It does NOT run any new analysis: it simply reformats the existing
`04_validate_annotations/validated_results.csv` of each species into a
reader-friendly, well-documented table of the MODERATE-tier candidates
(the set underlying the abstract, Table 2, and all figures), containing,
for every predicted candidate:
    - protein identifier
    - assigned defense category
    - confidence (max Z-score, max percentile, stringency tier)
    - per-category Z-scores
    - keyword-validation status and the matched keyword evidence
    - annotation_support class
    - RefSeq functional description

Run from the project root (the folder that contains `results/`):
    python 08_make_supplementary_tables.py

Output -> results/08_supplementary_tables/
    PlantDefenseESM_Supplementary_Candidates.xlsx   (one sheet per species + dictionary + summary)
"""

from pathlib import Path
import sys
import pandas as pd

# species folder name  ->  display name used in the manuscript
SPECIES = {
    "arabidopsis_thaliana": "Arabidopsis thaliana",
    "oryza_sativa":         "Oryza sativa",
    "vitis_vinifera":       "Vitis vinifera",
}

RESULTS = Path("results")
OUTDIR = RESULTS / "08_supplementary_tables"
OUTDIR.mkdir(parents=True, exist_ok=True)


def highest_tier(row):
    """Strictest tier a protein qualifies for (strict > moderate > lenient)."""
    if row.get("defense_strict", False):
        return "strict"
    if row.get("defense_moderate", False):
        return "moderate"
    if row.get("defense_lenient", False):
        return "lenient"
    return "none"


def build_species_table(species_dir: str) -> pd.DataFrame:
    val_path = RESULTS / species_dir / "04_validate_annotations" / "validated_results.csv"
    if not val_path.exists():
        print(f"  [SKIP] {val_path} not found")
        return None

    df = pd.read_csv(val_path)
    # protein_id may be the index column
    if "protein_id" not in df.columns:
        df = df.rename(columns={df.columns[0]: "protein_id"})

    # candidates = moderate-tier set (matches abstract / Table 2 / figures)
    if "defense_moderate" in df.columns:
        cand = df[df["defense_moderate"]].copy()
    else:  # fall back: union of any defense_* flag
        flags = [c for c in df.columns if c.startswith("defense_") and df[c].dtype == bool]
        cand = df[df[flags].any(axis=1)].copy()

    cand["tier"] = cand.apply(highest_tier, axis=1)

    # per-category z-score columns, in a stable order
    z_cat_cols = [c for c in cand.columns if c.startswith("z_") and c != "z_max"]

    out = pd.DataFrame()
    out["protein_id"] = cand["protein_id"]
    out["assigned_category"] = cand.get("best_category")
    out["tier"] = cand["tier"]
    out["confidence_zmax"] = cand.get("z_max").round(3)
    if "pctl_max" in cand.columns:
        out["confidence_percentile"] = cand["pctl_max"].round(2)
    out["validation_status"] = cand.get("has_defense_keyword").map(
        {True: "annotation-supported", False: "no defense keyword"}
    )
    out["matched_keywords"] = cand.get("keyword_hits").fillna("")
    out["annotation_support_class"] = cand.get("annotation_support")
    out["refseq_description"] = cand.get("description").fillna("")
    # append per-category z-scores
    for c in z_cat_cols:
        out[c] = cand[c].round(3)

    out = out.sort_values("confidence_zmax", ascending=False).reset_index(drop=True)
    return out


DICTIONARY = [
    ("protein_id", "RefSeq protein accession (proteome identifier)."),
    ("assigned_category", "Defense category with the highest Z-score (one of six)."),
    ("tier", "Stringency tier passed: strict (>=99.5th pct) or moderate (>=99.0th). Table is restricted to moderate-tier candidates."),
    ("confidence_zmax", "Maximum per-category Z-score; primary confidence score."),
    ("confidence_percentile", "Maximum per-category percentile rank (0-100)."),
    ("validation_status", "'annotation-supported' if the RefSeq description contains >=1 defense keyword, else 'no defense keyword'."),
    ("matched_keywords", "Defense keyword(s) found in the RefSeq description (evidence; '|'-separated)."),
    ("annotation_support_class", "annotation_supported_candidate (candidate WITH a defense keyword) or candidate_without_defense_keyword (candidate WITHOUT any defense keyword)."),
    ("refseq_description", "Full RefSeq functional annotation string."),
    ("z_<category>", "Z-score of this protein against each of the six defense category centroids."),
]


def main():
    tables = {}
    summary_rows = []

    print("=" * 60)
    print("Building supplementary candidate tables")
    print("=" * 60)

    for folder, display in SPECIES.items():
        print(f"\n{display}  ({folder})")
        tbl = build_species_table(folder)
        if tbl is None:
            continue
        tables[display] = tbl

        n = len(tbl)
        by_tier = tbl["tier"].value_counts().to_dict()
        n_novel = int((tbl["annotation_support_class"] == "candidate_without_defense_keyword").sum())
        print(f"  rows: {n:,}")
        print(f"  by tier: {by_tier}")
        print(f"  candidate_without_defense_keyword: {n_novel:,}")
        print(f"  columns: {list(tbl.columns)}")
        print("  preview (top 3 rows):")
        with pd.option_context("display.max_colwidth", 45, "display.width", 200):
            print(tbl.head(3).to_string(index=False))

        summary_rows.append({
            "Species": display,
            "Moderate-tier candidates": n,
            "Strict": by_tier.get("strict", 0),
            "Moderate only": by_tier.get("moderate", 0),
            "Candidates without a defense keyword": n_novel,
            "Annotation-supported": int((tbl["validation_status"] == "annotation-supported").sum()),
        })

    # combined Excel workbook
    xlsx_path = OUTDIR / "PlantDefenseESM_Supplementary_Candidates.xlsx"
    try:
        with pd.ExcelWriter(xlsx_path, engine="openpyxl") as xl:
            pd.DataFrame(summary_rows).to_excel(xl, sheet_name="Summary", index=False)
            pd.DataFrame(DICTIONARY, columns=["Column", "Description"]).to_excel(
                xl, sheet_name="Data dictionary", index=False)
            for display, tbl in tables.items():
                sheet = display.replace(" ", "_")[:31]
                tbl.to_excel(xl, sheet_name=sheet, index=False)
        print(f"\nExcel workbook -> {xlsx_path}")
    except Exception as e:
        print(f"\n[ERROR] Could not write Excel ({e}). "
              f"Install openpyxl with: pip install openpyxl, then re-run.")

    print("\nDone.")


if __name__ == "__main__":
    if not RESULTS.exists():
        sys.exit("ERROR: run this from the project root (no 'results/' folder here).")
    main()
