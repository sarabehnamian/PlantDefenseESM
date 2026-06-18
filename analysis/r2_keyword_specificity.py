#!/usr/bin/env python3
"""
r2_keyword_specificity.py
=========================
Keyword-specificity / annotation-bias sensitivity analysis for the
PlantDefenseESM revision (Reviewer 2: keyword validation is heuristic and
broad terms can inflate enrichment; Reviewer 3: split the validation
keywords into broad vs defense-specific and report them as a table).

Runs ENTIRELY from cached step-04 outputs. No GPU, no re-embedding.

The 47 validation keywords are partitioned into two groups:
  - BROAD   : generic terms that also occur widely outside innate immunity
              (e.g. "kinase", "peroxidase", "oxidase", "receptor-like",
              "lignin", "ethylene").
  - SPECIFIC: defense-specific terms unlikely to be hit by housekeeping
              proteins (e.g. "NB-ARC", "pathogenesis-related", "defensin",
              "powdery mildew").

For each species it re-derives the keyword-validation flag under three
keyword sets (ALL 47, SPECIFIC-only, BROAD-only) and recomputes, at each
stringency tier, the candidate keyword rate, background rate, fold
enrichment and one-sided Fisher's exact p-value.  If the enrichment holds
(or strengthens) under SPECIFIC-only keywords, the signal is not an
artefact of broad terms.

It also recomputes the moderate-tier novelty counts under SPECIFIC-only
keywords, since "novel" = candidate without ANY defense keyword and the
broad terms make that definition more conservative.

Run from the project root (the folder that contains `results/`):
    python3 r2_keyword_specificity.py
    # or point it at custom dirs:
    python3 r2_keyword_specificity.py results/arabidopsis_thaliana results/oryza_sativa

Outputs (written next to this script):
    r2_keyword_groups.csv                 - the broad/specific keyword table
    r2_keyword_enrichment_<species>.csv   - enrichment per keyword set x tier
    r2_keyword_specificity_summary.csv    - one combined tidy table
    r2_keyword_novelty_<species>.csv      - novelty counts: ALL vs SPECIFIC
"""

import sys
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact

# Default species result directories (edit if yours differ)
DEFAULT_DIRS = [
    "results/arabidopsis_thaliana",
    "results/oryza_sativa",
    "results/vitis_vinifera",
]

# Partition of the 47 VALIDATION_KEYWORDS (shared.py) into broad vs specific.
# BROAD = generic terms that frequently occur outside innate immunity and
# could inflate enrichment; SPECIFIC = defense-specific terms.
BROAD_KEYWORDS = [
    "kinase", "receptor-like",
    "peroxidase", "oxidase",
    "lignin", "ethylene",
    "resistance", "disease",
]

SPECIFIC_KEYWORDS = [
    "defense", "defence", "pathogen",
    "NBS", "NB-ARC", "LRR", "TIR", "CC-NBS",
    "PR-", "pathogenesis-related", "chitinase", "glucanase",
    "thaumatin", "defensin", "osmotin",
    "stilbene synthase", "STS", "resveratrol",
    "phenylalanine ammonia-lyase", "PAL",
    "WRKY", "NPR1", "EDS1", "PAD4",
    "callose",
    "jasmonic", "salicylic",
    "hypersensitive", "programmed cell death",
    "R gene", "R-gene", "immune", "immunity",
    "downy mildew", "powdery mildew", "Botrytis", "Plasmopara",
    "Erysiphe", "Phytophthora",
]

KEYWORD_SETS = {
    "all": BROAD_KEYWORDS + SPECIFIC_KEYWORDS,
    "specific": SPECIFIC_KEYWORDS,
    "broad": BROAD_KEYWORDS,
}


def keyword_scan(description: str, keywords: list) -> bool:
    """True if any keyword (case-insensitive substring) appears in description."""
    desc_lower = str(description).lower()
    return any(kw.lower() in desc_lower for kw in keywords)


def enrichment_rows(df: pd.DataFrame, species: str, set_name: str,
                    keywords: list) -> list:
    """Fisher enrichment at each tier for one keyword set."""
    has_kw = df["description"].apply(lambda d: keyword_scan(d, keywords))
    n_total = len(df)
    n_kw = int(has_kw.sum())
    bg_rate = n_kw / n_total if n_total else 0.0

    rows = []
    for tier in ("strict", "moderate", "lenient"):
        col = f"defense_{tier}"
        if col not in df.columns:
            continue
        mask = df[col].astype(bool)
        n_cand = int(mask.sum())
        n_cand_kw = int((mask & has_kw).sum())
        cand_rate = n_cand_kw / n_cand if n_cand else 0.0
        enrichment = cand_rate / bg_rate if bg_rate > 0 else float("inf")

        a = n_cand_kw
        b = n_cand - n_cand_kw
        c = n_kw - n_cand_kw
        d = n_total - n_cand - c
        odds, pval = fisher_exact([[a, b], [c, d]], alternative="greater")

        rows.append({
            "species": species,
            "keyword_set": set_name,
            "n_keywords": len(keywords),
            "tier": tier,
            "n_candidates": n_cand,
            "n_with_keyword": n_cand_kw,
            "candidate_rate": round(cand_rate, 4),
            "background_rate": round(bg_rate, 4),
            "fold_enrichment": round(enrichment, 2),
            "odds_ratio": round(float(odds), 2),
            "fisher_p": f"{pval:.2e}",
        })
    return rows


def novelty_rows(df: pd.DataFrame, species: str) -> list:
    """Moderate-tier novelty counts under ALL vs SPECIFIC-only keywords."""
    out = []
    mod = df["defense_moderate"].astype(bool) if "defense_moderate" in df.columns \
        else pd.Series(False, index=df.index)
    n_mod = int(mod.sum())
    for set_name in ("all", "specific"):
        has_kw = df["description"].apply(
            lambda d: keyword_scan(d, KEYWORD_SETS[set_name]))
        known = int((mod & has_kw).sum())
        novel = int((mod & ~has_kw).sum())
        out.append({
            "species": species,
            "keyword_set": set_name,
            "moderate_candidates": n_mod,
            "known_defense": known,
            "novel_candidate": novel,
            "novel_fraction": round(novel / n_mod, 4) if n_mod else 0.0,
        })
    return out


def analyse_species(result_dir: Path) -> tuple:
    species = result_dir.name
    val_path = result_dir / "04_validate_annotations" / "validated_results.csv"
    df = pd.read_csv(val_path)
    if "protein_id" not in df.columns:
        df = df.rename(columns={df.columns[0]: "protein_id"})
    if "description" not in df.columns:
        df["description"] = ""
    df["description"] = df["description"].fillna("")

    enr = []
    for set_name, kws in KEYWORD_SETS.items():
        enr.extend(enrichment_rows(df, species, set_name, kws))
    enr_df = pd.DataFrame(enr)
    out = Path(__file__).resolve().parent / f"r2_keyword_enrichment_{species}.csv"
    enr_df.to_csv(out, index=False)

    nov_df = pd.DataFrame(novelty_rows(df, species))
    nov_out = Path(__file__).resolve().parent / f"r2_keyword_novelty_{species}.csv"
    nov_df.to_csv(nov_out, index=False)

    # Compact console view
    print(f"\n===== {species}  (proteome N={len(df):,}) =====")
    print("  Fold enrichment by keyword set (moderate tier):")
    sub = enr_df[enr_df.tier == "moderate"][
        ["keyword_set", "n_keywords", "background_rate",
         "candidate_rate", "fold_enrichment", "fisher_p"]]
    print(sub.to_string(index=False))
    print("\n  Moderate-tier novelty (ALL vs SPECIFIC keywords):")
    print(nov_df[["keyword_set", "known_defense", "novel_candidate",
                  "novel_fraction"]].to_string(index=False))

    return enr_df, nov_df


def main():
    # Write the broad/specific keyword table once
    here = Path(__file__).resolve().parent
    groups = (
        [{"keyword": k, "group": "specific"} for k in SPECIFIC_KEYWORDS]
        + [{"keyword": k, "group": "broad"} for k in BROAD_KEYWORDS]
    )
    pd.DataFrame(groups).to_csv(here / "r2_keyword_groups.csv", index=False)
    print(f"Wrote keyword groups table -> r2_keyword_groups.csv "
          f"({len(SPECIFIC_KEYWORDS)} specific, {len(BROAD_KEYWORDS)} broad)")

    dirs = sys.argv[1:] if len(sys.argv) > 1 else DEFAULT_DIRS
    all_enr = []
    for d in dirs:
        p = Path(d)
        if not (p / "04_validate_annotations" / "validated_results.csv").exists():
            print(f"[skip] {d}: no 04_validate_annotations/validated_results.csv")
            continue
        enr_df, _ = analyse_species(p)
        all_enr.append(enr_df)

    if all_enr:
        combined = pd.concat(all_enr, ignore_index=True)
        out = here / "r2_keyword_specificity_summary.csv"
        combined.to_csv(out, index=False)
        print(f"\nWrote combined summary -> {out}")
        print("Per-species tables -> r2_keyword_enrichment_<species>.csv, "
              "r2_keyword_novelty_<species>.csv")


if __name__ == "__main__":
    if not Path("results").exists():
        sys.exit("ERROR: run from the project root (no 'results/' folder here).")
    main()
