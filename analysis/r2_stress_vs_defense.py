#!/usr/bin/env python3
"""
r2_stress_vs_defense.py
=======================
Distinguish "defense-related" from "stress-related / developmental / metabolic"
among the moderate-tier candidates (Reviewer 3: keyword/category overlap).

Runs ENTIRELY from cached step-04 outputs. No GPU, no re-embedding.

For each species it scans the RefSeq description of every moderate-tier
candidate and tags it with three (non-exclusive) annotation flags:

  * defense_specific : carries a strict biotic-defense term
                       (the 39-term defense-specific subset)
  * broad_stress     : carries a generic/broad term that overlaps abiotic
                       stress, development, or general metabolism
                       (drought, cold, auxin, growth, kinase, oxidase, ...)
  * unannotated      : carries no validation keyword at all

It then reports, per species:
  - how many moderate candidates are defense-specific only,
  - how many are flagged ONLY by a broad/overlapping term (the proteins the
    reviewer worries are stress/development, not defense),
  - how many carry an explicit abiotic-stress / developmental term,
  - the same for the keyword-negative ("novel") subset.

Outputs (written next to this script):
    r2_stress_vs_defense_<species>.csv   - per-protein flags (candidates only)
    r2_stress_vs_defense_summary.csv     - one combined tidy summary table
"""

import sys
from pathlib import Path

import pandas as pd

DEFAULT_DIRS = [
    "results/arabidopsis_thaliana",
    "results/oryza_sativa",
    "results/vitis_vinifera",
]

# --- term sets -------------------------------------------------------------
# 39 defense-specific terms (the strict biotic-defense subset). These are the
# validation keywords MINUS the 8 broad/generic terms identified in the revision.
DEFENSE_SPECIFIC = [
    "disease", "defense", "defence", "pathogen",
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

# 8 broad/generic terms already separated out in the revision. These overlap
# heavily with abiotic stress, development, and general metabolism.
BROAD_TERMS = [
    "kinase", "receptor-like", "peroxidase", "oxidase",
    "lignin", "ethylene", "resistance", "disease",
]
# NOTE: "disease"/"resistance" appear in both BMC lists historically; for the
# overlap analysis we treat the 8-term broad set exactly as published.
BROAD_TERMS = ["kinase", "receptor-like", "peroxidase", "oxidase",
               "lignin", "ethylene", "resistance", "disease"]

# Explicit abiotic-stress / developmental / general-metabolism vocabulary.
# A candidate flagged here, but NOT by any defense-specific term, is the kind
# of overlap the reviewer asks us to distinguish.
ABIOTIC_DEV_META = [
    # abiotic stress
    "drought", "cold", "heat", "salt", "salinity", "osmotic", "dehydration",
    "wound", "wounding", "abscisic", "ABA", "freezing", "chilling", "hypoxia",
    "oxidative stress", "heavy metal", "UV",
    # development / growth / hormone-growth
    "auxin", "cytokinin", "gibberellin", "brassinosteroid",
    "growth", "development", "developmental", "flowering", "meristem",
    "root", "pollen", "seed", "embryo", "senescence",
    # general metabolism
    "metabolic", "metabolism", "biosynthesis", "transferase",
    "hydrolase", "dehydrogenase", "transporter", "transport",
    "photosynth", "ribosom",
]


def scan(desc: str, terms) -> bool:
    d = (desc or "").lower()
    return any(t.lower() in d for t in terms)


def analyse_species(result_dir: Path) -> pd.DataFrame:
    species = result_dir.name
    vr = pd.read_csv(
        result_dir / "04_validate_annotations" / "validated_results.csv"
    )

    # restrict to moderate-tier candidates
    cand = vr[vr["defense_moderate"]].copy()
    desc = cand["description"].fillna("")

    cand["is_defense_specific"] = desc.apply(lambda d: scan(d, DEFENSE_SPECIFIC))
    cand["is_broad"] = desc.apply(lambda d: scan(d, BROAD_TERMS))
    cand["is_abiotic_dev_meta"] = desc.apply(lambda d: scan(d, ABIOTIC_DEV_META))
    cand["has_keyword"] = cand["has_defense_keyword"]

    # mutually informative buckets
    cand["bucket"] = "unannotated"
    cand.loc[cand["is_broad"] & ~cand["is_defense_specific"], "bucket"] = "broad_only"
    cand.loc[cand["is_defense_specific"], "bucket"] = "defense_specific"

    out = Path(__file__).resolve().parent / f"r2_stress_vs_defense_{species}.csv"
    cols = ["protein_id", "description", "keyword_hits",
            "is_defense_specific", "is_broad", "is_abiotic_dev_meta", "bucket"]
    cols = [c for c in cols if c in cand.columns]
    cand[cols].to_csv(out, index=False)

    n = len(cand)
    n_def = int(cand["is_defense_specific"].sum())
    n_broad_only = int((cand["is_broad"] & ~cand["is_defense_specific"]).sum())
    n_abiotic = int(cand["is_abiotic_dev_meta"].sum())
    n_abiotic_only = int(
        (cand["is_abiotic_dev_meta"] & ~cand["is_defense_specific"]).sum())
    n_unann = int((~cand["has_keyword"]).sum())

    # keyword-negative ("novel") subset: how many actually look abiotic/dev/meta?
    kwneg = cand[~cand["has_keyword"]]
    n_kwneg = len(kwneg)
    n_kwneg_abiotic = int(kwneg["is_abiotic_dev_meta"].sum())

    print(f"\n===== {species}  (moderate-tier candidates N={n:,}) =====")
    print(f"  defense-specific term present : {n_def:,} ({n_def/n:.1%})")
    print(f"  broad term ONLY (no specific) : {n_broad_only:,} ({n_broad_only/n:.1%})")
    print(f"  abiotic/dev/metabolism term   : {n_abiotic:,} ({n_abiotic/n:.1%})")
    print(f"    ... of which NOT defense-spec: {n_abiotic_only:,} ({n_abiotic_only/n:.1%})")
    print(f"  no validation keyword at all  : {n_unann:,} ({n_unann/n:.1%})")
    print(f"  keyword-negative subset       : N={n_kwneg:,}; "
          f"{n_kwneg_abiotic:,} ({(n_kwneg_abiotic/n_kwneg if n_kwneg else 0):.1%}) "
          f"carry an abiotic/dev/metabolism term")

    return pd.DataFrame([{
        "species": species,
        "n_moderate_candidates": n,
        "n_defense_specific": n_def,
        "pct_defense_specific": round(100 * n_def / n, 2),
        "n_broad_only": n_broad_only,
        "pct_broad_only": round(100 * n_broad_only / n, 2),
        "n_abiotic_dev_meta": n_abiotic,
        "pct_abiotic_dev_meta": round(100 * n_abiotic / n, 2),
        "n_abiotic_only_not_defense": n_abiotic_only,
        "pct_abiotic_only_not_defense": round(100 * n_abiotic_only / n, 2),
        "n_unannotated": n_unann,
        "pct_unannotated": round(100 * n_unann / n, 2),
        "n_keyword_negative": n_kwneg,
        "n_keyword_negative_abiotic": n_kwneg_abiotic,
        "pct_keyword_negative_abiotic": round(
            100 * n_kwneg_abiotic / n_kwneg, 2) if n_kwneg else 0.0,
    }])


def main():
    dirs = sys.argv[1:] if len(sys.argv) > 1 else DEFAULT_DIRS
    all_dfs = []
    for d in dirs:
        p = Path(d)
        vr = p / "04_validate_annotations" / "validated_results.csv"
        if not vr.exists():
            print(f"[skip] {vr} not found", file=sys.stderr)
            continue
        all_dfs.append(analyse_species(p))
    if all_dfs:
        summary = pd.concat(all_dfs, ignore_index=True)
        out = Path(__file__).resolve().parent / "r2_stress_vs_defense_summary.csv"
        summary.to_csv(out, index=False)
        print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
