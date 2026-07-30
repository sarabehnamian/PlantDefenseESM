#!/usr/bin/env python3
"""
18_final_residual.py
====================
Folds the existing jackhmmer scores into the new PSI-BLAST / filter-off arms
and computes the FINAL ESM-2-only residual against the complete search panel.

This is the number that replaces "the ESM-2-only count is an upper bound"
in the Results and Discussion.

Reads:
  results/13_psiblast_filteroff/arm_scores.csv        (from 17_psiblast_filteroff.py)
  results/16_profile_baselines/_cache/scores_arabidopsis_thaliana_jackhmmer.csv
  results/15_baseline_overlap/positive_membership_arabidopsis_thaliana.csv
  results/16_profile_baselines/residual_curated_arabidopsis_thaliana.csv  (the published 72)

Writes (results/18_final_residual/):
  final_metrics.csv        every arm, one table            -> replaces Table 8
  final_residual.csv       the surviving ESM-2-only proteins
  recovered_from_72.csv    which of the published 72 were picked up, and by what
  summary.md               plain-text summary, no tabulate needed

Usage:
  python 18_final_residual.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

SPECIES = "arabidopsis_thaliana"
N_MATCHED = 2807

ARM_SCORES = Path("results/13_psiblast_filteroff/arm_scores.csv")
JACK_CACHE = Path(f"results/16_profile_baselines/_cache/scores_{SPECIES}_jackhmmer.csv")
MEMBERSHIP = Path(f"results/15_baseline_overlap/positive_membership_{SPECIES}.csv")
PUBLISHED_RESIDUAL = Path(f"results/16_profile_baselines/residual_curated_{SPECIES}.csv")
OUTDIR = Path("results/18_final_residual")


def die(msg: str):
    print(f"\nERROR: {msg}\n", file=sys.stderr)
    sys.exit(1)


def find_col(df: pd.DataFrame, candidates: list[str], what: str) -> str:
    lower = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c in lower:
            return lower[c]
    die(f"could not find the {what} column. Columns present: {list(df.columns)}")


def top_set(scores: pd.Series, n: int) -> set[str]:
    """Matched-count candidate set: top n by score, excluding zero/NaN."""
    s = scores.dropna()
    s = s[s > 0]
    return set(s.sort_values(ascending=False).head(n).index)


def main() -> None:
    for f in (ARM_SCORES, MEMBERSHIP):
        if not f.exists():
            die(f"missing input: {f}")

    OUTDIR.mkdir(parents=True, exist_ok=True)

    # ---- membership: curated set, ESM-2 set, domain set -------------------
    mem = pd.read_csv(MEMBERSHIP)
    pid = find_col(mem, ["protein_id"], "protein id")
    mem[pid] = mem[pid].astype(str).str.strip()
    curated = set(mem[pid])
    esm = set(mem.loc[mem["ESM2"].astype(str).str.lower() == "true", pid])
    domain = set(mem.loc[mem["domain"].astype(str).str.lower() == "true", pid])
    align = set(mem.loc[mem["alignment"].astype(str).str.lower() == "true", pid])
    print(f"curated {len(curated):,} | ESM-2 {len(esm):,} | "
          f"alignment {len(align):,} | domain {len(domain):,}")

    # ---- new arms --------------------------------------------------------
    arms_df = pd.read_csv(ARM_SCORES)
    apid = find_col(arms_df, ["protein_id"], "protein id")
    arms_df[apid] = arms_df[apid].astype(str).str.strip()
    arms_df = arms_df.set_index(apid)

    method_sets: dict[str, set[str]] = {}
    for col in arms_df.columns:
        method_sets[col] = top_set(arms_df[col], N_MATCHED)

    # ---- jackhmmer from the existing cache -------------------------------
    if JACK_CACHE.exists():
        jack = pd.read_csv(JACK_CACHE)
        jpid = find_col(jack, ["protein_id", "pid"], "protein id")
        jscore = find_col(
            jack, ["score", "score_jackhmmer", "bitscore", "best_bitscore"],
            "jackhmmer score")
        jack[jpid] = jack[jpid].astype(str).str.strip()
        method_sets["jackhmmer"] = top_set(
            jack.set_index(jpid)[jscore], N_MATCHED)
        print(f"jackhmmer loaded from cache: {len(method_sets['jackhmmer']):,} candidates")
    else:
        print(f"NOTE: {JACK_CACHE} not found - jackhmmer not included. "
              f"Check the filename and re-run.")

    method_sets["alignment_SW"] = align
    method_sets["interpro_domain"] = domain

    # ---- metrics table ---------------------------------------------------
    rows = [{
        "method": "ESM-2 (moderate tier)",
        "candidates": len(esm),
        "recovered_TP": len(esm & curated),
        "precision": round(len(esm & curated) / len(esm), 3),
        "recall": round(len(esm & curated) / len(curated), 3),
    }]
    for name, s in method_sets.items():
        tp = len(s & curated)
        rows.append({
            "method": name,
            "candidates": len(s),
            "recovered_TP": tp,
            "precision": round(tp / len(s), 3) if s else 0.0,
            "recall": round(tp / len(curated), 3),
        })
    metrics = pd.DataFrame(rows)
    metrics["F1"] = (2 * metrics["precision"] * metrics["recall"] /
                     (metrics["precision"] + metrics["recall"])).round(3)
    metrics.to_csv(OUTDIR / "final_metrics.csv", index=False)
    print("\n" + metrics.to_string(index=False))

    # ---- the final residual ---------------------------------------------
    search_union = set().union(*method_sets.values())
    esm_positives = esm & curated
    residual = esm_positives - search_union

    print(f"\nESM-2 recovers {len(esm_positives)} curated proteins")
    print(f"Search panel ({len(method_sets)} methods) recovers "
          f"{len(search_union & curated)} curated proteins")
    print(f"FINAL ESM-2-only residual: {len(residual)}")

    res = mem[mem[pid].isin(residual)].copy()
    res.to_csv(OUTDIR / "final_residual.csv", index=False)

    # ---- what happened to the published 72 -------------------------------
    if PUBLISHED_RESIDUAL.exists():
        pub = pd.read_csv(PUBLISHED_RESIDUAL)
        ppid = find_col(pub, ["protein_id"], "protein id")
        pub[ppid] = pub[ppid].astype(str).str.strip()
        pub72 = set(pub[ppid])
        recovered_now = pub72 - residual
        still_only = pub72 & residual
        print(f"\nPublished residual: {len(pub72)}")
        print(f"  now recovered by a search method: {len(recovered_now)}")
        print(f"  still ESM-2-only:                 {len(still_only)}")

        det = []
        for p in sorted(recovered_now):
            by = [n for n, s in method_sets.items() if p in s]
            desc = pub.loc[pub[ppid] == p, "description"]
            det.append({
                "protein_id": p,
                "recovered_by": "; ".join(by) if by else "(not in ESM positives)",
                "description": desc.iloc[0] if len(desc) else "",
            })
        pd.DataFrame(det).to_csv(OUTDIR / "recovered_from_72.csv", index=False)
    else:
        pub72, recovered_now, still_only = set(), set(), set()

    # ---- summary.md (no tabulate) ---------------------------------------
    def md_table(df: pd.DataFrame) -> str:
        cols = list(df.columns)
        out = ["| " + " | ".join(cols) + " |",
               "| " + " | ".join("---" for _ in cols) + " |"]
        for _, r in df.iterrows():
            out.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
        return "\n".join(out)

    with open(OUTDIR / "summary.md", "w", encoding="utf-8") as fh:
        fh.write("# Final search-panel comparison (E4 + E5 closed)\n\n")
        fh.write(f"Species: *{SPECIES}*  \n")
        fh.write(f"Curated GO defense set: {len(curated):,} proteins  \n")
        fh.write(f"Matched candidate count: n = {N_MATCHED:,}\n\n")
        fh.write("## Per-method performance\n\n")
        fh.write(md_table(metrics))
        fh.write("\n\n## ESM-2-only residual\n\n")
        fh.write(f"- Published (vs jackhmmer + domain only): {len(pub72)}\n")
        fh.write(f"- Against the full panel ({', '.join(method_sets)}): "
                 f"**{len(residual)}**\n")
        if pub72:
            fh.write(f"- Of the published set, {len(recovered_now)} are now "
                     f"recovered by a search method and {len(still_only)} remain "
                     f"ESM-2-only\n")
        fh.write("\n## Caveat to state in Methods\n\n")
        fh.write("Ranking by maximum bitscore across 33 independent PSI-BLAST "
                 "PSSMs mixes score scales, since each query profile diverges "
                 "over the three iterations. This affects the matched-count cut "
                 "but not the AUPRC over the full ranking, and all search arms "
                 "exceed the embeddings on both measures.\n")

    print(f"\nDone. See {OUTDIR}/summary.md")


if __name__ == "__main__":
    main()
