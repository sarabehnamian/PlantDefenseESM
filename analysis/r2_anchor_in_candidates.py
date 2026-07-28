#!/usr/bin/env python3
"""
r2_anchor_in_candidates.py
==========================
Reviewer-3 check: are the 33 anchor proteins (or near-identical homologs of
them) included in the candidate counts and enrichment analyses?

Anchors are UniProt entries; the scored proteomes are RefSeq, so an anchor can
never enter the candidate set by ID. The real question is whether RefSeq
ORTHOLOGS / near-identical homologs of the (mostly Arabidopsis) anchors land in
the moderate-tier candidate set, and whether removing them changes enrichment.

For each species this script:
  1. Loads moderate-tier candidates from validated_results.csv.
  2. Aligns every anchor against every candidate (global, BLOSUM62) and flags a
     candidate as an "anchor homolog" at >=98% and >=95% identity.
  3. Reports how many candidates are anchor homologs.
  4. Recomputes the moderate-tier keyword fold-enrichment with those homologs
     removed, to show the enrichment does not depend on anchor self-inclusion.

Run from project root. Pure pip (biopython only). No GPU.
Output -> r2_anchor_in_candidates.csv
"""

from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

from shared import CORE_ANCHORS, load_config

SPECIES_CONFIGS = {
    "arabidopsis_thaliana": "config_arabidopsis.yaml",
    "oryza_sativa": "config_rice.yaml",
    "vitis_vinifera": "config_grapevine.yaml",
}

ANCHOR_FASTA = Path("results/01_fetch_anchors/anchors.fasta")
ANCHOR_CACHE = Path("results/01_fetch_anchors/_uniprot_cache")
ID_THRESHOLDS = (0.98, 0.95)


def load_anchor_seqs():
    """Load anchor sequences from anchors.fasta (>UID|CAT|name).

    Anchors are identical across species, so use the first species folder
    that has the file."""
    candidates = [
        Path("results") / sp / "01_fetch_anchors" / "anchors.fasta"
        for sp in SPECIES_CONFIGS
    ]
    seqs = {}
    for fp in candidates:
        if fp.exists():
            for rec in SeqIO.parse(str(fp), "fasta"):
                uid = rec.id.split("|")[0]
                seqs[uid] = str(rec.seq).replace("*", "")
            return seqs
    return seqs


def make_aligner():
    a = PairwiseAligner()
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score = -11
    a.extend_gap_score = -1
    a.mode = "global"
    return a


def pct_identity(aligner, s1, s2):
    """Fast global %identity over the shorter sequence length."""
    aln = aligner.align(s1, s2)[0]
    a, b = aln[0], aln[1]
    matches = sum(1 for x, y in zip(a, b) if x == y and x != "-")
    return matches / min(len(s1), len(s2))


def main():
    aligner = make_aligner()
    anchor_seqs = load_anchor_seqs()
    print(f"Loaded {len(anchor_seqs)} anchor sequences\n")
    if len(anchor_seqs) == 0:
        raise SystemExit(
            "ERROR: 0 anchors loaded. Check that "
            "results/01_fetch_anchors/anchors.fasta exists."
        )

    rows = []
    for species, cfg_name in SPECIES_CONFIGS.items():
        cfg = load_config(cfg_name)
        base = Path(cfg["base_output_dir"])

        val = pd.read_csv(
            base / "04_validate_annotations" / "validated_results.csv",
            index_col="protein_id",
        )
        cands = val[val["defense_moderate"]].copy()

        # background keyword rate (whole proteome)
        bg_rate = val["has_defense_keyword"].mean()

        # candidate sequences
        prot_fasta = base / "00_download_proteome" / "proteome.fasta"
        cand_ids = set(cands.index.astype(str))
        cand_seqs = {}
        for rec in SeqIO.parse(str(prot_fasta), "fasta"):
            rid = rec.id
            if rid in cand_ids:
                cand_seqs[rid] = str(rec.seq).replace("*", "")

        # align anchors vs candidates; flag homologs
        best_id = {cid: 0.0 for cid in cand_seqs}
        for uid, aseq in anchor_seqs.items():
            for cid, cseq in cand_seqs.items():
                # length pre-filter: skip if lengths differ by >25%
                lo, hi = sorted((len(aseq), len(cseq)))
                if lo / hi < 0.75:
                    continue
                pid = pct_identity(aligner, aseq, cseq)
                if pid > best_id[cid]:
                    best_id[cid] = pid

        n_cand = len(cands)
        for thr in ID_THRESHOLDS:
            homolog_ids = {cid for cid, p in best_id.items() if p >= thr}
            n_hom = len(homolog_ids)

            kept = cands[~cands.index.astype(str).isin(homolog_ids)]
            n_kept = len(kept)
            rate_all = cands["has_defense_keyword"].mean()
            rate_kept = kept["has_defense_keyword"].mean()

            rows.append({
                "species": species,
                "id_threshold": thr,
                "n_candidates_moderate": n_cand,
                "n_anchor_homologs": n_hom,
                "pct_candidates_homolog": round(100 * n_hom / n_cand, 3),
                "fold_enrich_all": round(rate_all / bg_rate, 3),
                "fold_enrich_homologs_removed": round(rate_kept / bg_rate, 3),
                "n_candidates_after_removal": n_kept,
            })
            print(f"{species:22s} id>={thr:.2f}  "
                  f"homologs={n_hom:3d}/{n_cand}  "
                  f"fold {rate_all/bg_rate:.2f} -> {rate_kept/bg_rate:.2f}")

    out = pd.DataFrame(rows)
    out.to_csv("r2_anchor_in_candidates.csv", index=False)
    print("\nWrote r2_anchor_in_candidates.csv")


if __name__ == "__main__":
    main()
