#!/usr/bin/env python3
"""
11_anchor_robustness.py
=======================
Leave-one-anchor-out (LOAO) and leave-three-out robustness analysis for the
ESM-2 defense candidate predictions.  (Reviewer 1 Comment 4; Reviewer 3
major comment on anchor-set bias / candidate stability.)

The selection pipeline (step 03) is reproduced exactly from the CACHED
embeddings (no GPU / no re-embedding needed): category centroids are the
L2-normalized mean of the anchor embeddings in each category; proteins are
scored by cosine similarity; per-category Z-scores and percentile ranks are
computed; and the moderate-tier candidate set is the union of the 99.0th
percentile threshold and the top-50-per-category rule.

For each leave-out run we recompute the moderate candidate set and compare
it to the full-anchor (baseline) run:
    - Jaccard overlap of candidate sets
    - fraction of baseline candidates retained
    - number of candidates added / removed
    - Spearman correlation of the proteome-wide max-Z ranking

Run from the project root:
    python 11_anchor_robustness.py

Outputs -> results/11_anchor_robustness/
    loao_per_anchor_<species>.csv   - effect of removing each single anchor
    robustness_summary.csv          - per-species summary (LOAO + leave-3-out)
    summary.txt
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy import stats as sp_stats

SPECIES = {
    "arabidopsis_thaliana": "Arabidopsis thaliana",
    "oryza_sativa":         "Oryza sativa",
    "vitis_vinifera":       "Vitis vinifera",
}

PCTL_MODERATE = 99.0
TOP_N = 50
N_LEAVE3 = 50          # leave-three-out repetitions
SEED = 42

RESULTS = Path("results")
OUTDIR = RESULTS / "11_anchor_robustness"
OUTDIR.mkdir(parents=True, exist_ok=True)


def load_npz(path):
    d = np.load(path, allow_pickle=True)
    return d["embeddings"], list(d["protein_ids"])


def moderate_set(proteome_norm, anchor_emb, anchor_cat, keep_mask):
    """Reproduce step-03 moderate candidate set + max-Z ranking for a subset
    of anchors (keep_mask over anchors). Returns (candidate_id_index_set, z_max array)."""
    cats = sorted(set(anchor_cat))
    cents = []
    for cat in cats:
        idx = [i for i in range(len(anchor_cat))
               if anchor_cat[i] == cat and keep_mask[i]]
        c = anchor_emb[idx].mean(axis=0)
        c = c / np.linalg.norm(c)
        cents.append(c)
    C = np.stack(cents)                      # (n_cat, D)
    sim = proteome_norm @ C.T                # (N, n_cat)

    # per-category z and percentile
    z = sp_stats.zscore(sim, axis=0)
    z_max = z.max(axis=1)
    N = sim.shape[0]
    pctl = np.apply_along_axis(lambda col: sp_stats.rankdata(col) / N * 100, 0, sim)
    pctl_max = pctl.max(axis=1)

    pctl_hit = pctl_max >= PCTL_MODERATE
    top_union = np.zeros(N, dtype=bool)
    for j in range(sim.shape[1]):
        top_idx = np.argpartition(sim[:, j], -TOP_N)[-TOP_N:]
        top_union[top_idx] = True
    moderate = pctl_hit | top_union
    return set(np.where(moderate)[0]), z_max


def compare(base_set, base_rank, new_set, new_rank):
    inter = base_set & new_set
    union = base_set | new_set
    jac = len(inter) / len(union) if union else 1.0
    retained = len(inter) / len(base_set) if base_set else 1.0
    added = len(new_set - base_set)
    removed = len(base_set - new_set)
    rho = sp_stats.spearmanr(base_rank, new_rank).correlation
    return jac, retained, added, removed, rho


def main():
    rng = np.random.default_rng(SEED)
    lines = []
    def log(s=""):
        print(s); lines.append(s)

    summary_rows = []
    log("=" * 70)
    log("Anchor robustness: leave-one-out and leave-three-out")
    log("=" * 70)

    for folder, display in SPECIES.items():
        emb_dir = RESULTS / folder / "02_embed_proteins"
        a_path = emb_dir / "anchor_embeddings.npz"
        p_path = emb_dir / "proteome_embeddings.npz"
        meta_path = RESULTS / folder / "01_fetch_anchors" / "anchors_metadata.csv"
        if not (a_path.exists() and p_path.exists() and meta_path.exists()):
            log(f"\n[SKIP] {display}: missing embeddings or anchor metadata")
            continue

        log(f"\n### {display}")
        anchor_emb, anchor_ids = load_npz(a_path)
        proteome_emb, proteome_ids = load_npz(p_path)
        meta = pd.read_csv(meta_path)
        id2cat = dict(zip(meta["uniprot_id"], meta["category"]))
        anchor_cat = [id2cat.get(str(a).split("|")[0]) for a in anchor_ids]
        if any(c is None for c in anchor_cat):
            missing = [a for a, c in zip(anchor_ids, anchor_cat) if c is None]
            log(f"  [WARN] anchors without category mapping: {missing}")
            anchor_cat = [c if c is not None else "UNKNOWN" for c in anchor_cat]

        pn = proteome_emb / np.where(
            np.linalg.norm(proteome_emb, axis=1, keepdims=True) == 0, 1.0,
            np.linalg.norm(proteome_emb, axis=1, keepdims=True))

        n_anchor = len(anchor_ids)
        keep_all = np.ones(n_anchor, dtype=bool)
        base_set, base_rank = moderate_set(pn, anchor_emb, anchor_cat, keep_all)
        log(f"  anchors: {n_anchor} | baseline moderate candidates: {len(base_set):,}")

        # sanity vs the saved pipeline output
        vr = RESULTS / folder / "04_validate_annotations" / "validated_results.csv"
        if vr.exists():
            df = pd.read_csv(vr)
            if "defense_moderate" in df.columns:
                log(f"  (pipeline moderate count for reference: "
                    f"{int(df['defense_moderate'].sum()):,})")

        # ---- leave-one-anchor-out ----
        rows = []
        for i in range(n_anchor):
            km = keep_all.copy(); km[i] = False
            s, r = moderate_set(pn, anchor_emb, anchor_cat, km)
            jac, ret, add, rem, rho = compare(base_set, base_rank, s, r)
            rows.append({
                "removed_anchor": str(anchor_ids[i]).split("|")[0],
                "category": anchor_cat[i],
                "n_candidates": len(s),
                "jaccard": round(jac, 4),
                "retained_frac": round(ret, 4),
                "added": add, "removed": rem,
                "spearman_rank": round(rho, 4),
            })
        loao = pd.DataFrame(rows).sort_values("jaccard")
        loao.to_csv(OUTDIR / f"loao_per_anchor_{folder}.csv", index=False)

        worst = loao.iloc[0]
        log(f"  LOAO: Jaccard mean {loao['jaccard'].mean():.3f} "
            f"(min {loao['jaccard'].min():.3f}); "
            f"retained mean {loao['retained_frac'].mean():.3f} "
            f"(min {loao['retained_frac'].min():.3f}); "
            f"Spearman mean {loao['spearman_rank'].mean():.4f} "
            f"(min {loao['spearman_rank'].min():.4f})")
        log(f"  LOAO: most influential anchor = {worst['removed_anchor']} "
            f"({worst['category']}), Jaccard {worst['jaccard']:.3f}, "
            f"retained {worst['retained_frac']:.3f}")

        # ---- leave-three-out (random) ----
        l3 = []
        for _ in range(N_LEAVE3):
            km = keep_all.copy()
            drop = rng.choice(n_anchor, size=3, replace=False)
            km[drop] = False
            s, r = moderate_set(pn, anchor_emb, anchor_cat, km)
            jac, ret, add, rem, rho = compare(base_set, base_rank, s, r)
            l3.append((jac, ret, rho))
        l3 = np.array(l3)
        log(f"  Leave-3-out ({N_LEAVE3} runs): Jaccard mean {l3[:,0].mean():.3f} "
            f"(min {l3[:,0].min():.3f}); retained mean {l3[:,1].mean():.3f} "
            f"(min {l3[:,1].min():.3f}); Spearman mean {l3[:,2].mean():.4f}")

        summary_rows.append({
            "species": display,
            "n_anchors": n_anchor,
            "baseline_moderate": len(base_set),
            "loao_jaccard_mean": round(loao["jaccard"].mean(), 4),
            "loao_jaccard_min": round(loao["jaccard"].min(), 4),
            "loao_retained_mean": round(loao["retained_frac"].mean(), 4),
            "loao_retained_min": round(loao["retained_frac"].min(), 4),
            "loao_spearman_mean": round(loao["spearman_rank"].mean(), 4),
            "loao_spearman_min": round(loao["spearman_rank"].min(), 4),
            "most_influential_anchor": worst["removed_anchor"],
            "most_influential_category": worst["category"],
            "leave3_jaccard_mean": round(float(l3[:, 0].mean()), 4),
            "leave3_jaccard_min": round(float(l3[:, 0].min()), 4),
            "leave3_retained_mean": round(float(l3[:, 1].mean()), 4),
        })

    if summary_rows:
        pd.DataFrame(summary_rows).to_csv(OUTDIR / "robustness_summary.csv", index=False)
    (OUTDIR / "summary.txt").write_text("\n".join(lines), encoding="utf-8")
    log(f"\nOutputs -> {OUTDIR}")


if __name__ == "__main__":
    if not RESULTS.exists():
        sys.exit("ERROR: run from the project root (no 'results/' folder here).")
    main()
