#!/usr/bin/env python3
"""
14_truncation_sensitivity.py
============================
Does retaining the C-terminal regions of long proteins change the ESM-2
candidate sets, the proteome-wide ranking, or the recall of long immune
receptors?  (Reviewer 1 Comment 3, Option A.)

Compares three scoring arms on identical selection rules:

    A  truncated       - the published pipeline: every sequence >1,022 aa is
                         cut to its first 1,022 residues (step 02).
    B  windowed-mean   - long sequences represented by the length-weighted
                         mean of their per-window embeddings (step 13).
    C  windowed-max    - long sequences scored by the MAXIMUM cosine
                         similarity across their windows, per category, so a
                         defense signal confined to one region is not diluted
                         by the rest of the protein.

Category centroids in arms B and C are rebuilt from the window-mean
embeddings of the 6 truncated anchors; the 27 short anchors are unchanged.
(Max-over-windows is a scoring rule for candidates, not a way to define a
centroid, so arm C uses the same centroids as arm B.)

No GPU is used: everything runs from the cached step-02 and step-13
embeddings. Curated GO and InterPro family sets are reused from the caches
written by 09_benchmark_curated.py and 10_benchmark_families.py.

Run from the project root, after steps 09, 10 and 13:
    python 14_truncation_sensitivity.py

Outputs -> results/14_truncation_sensitivity/
    centroid_shift.csv              - cosine(baseline centroid, windowed centroid)
    arm_summary.csv                 - candidate counts, Jaccard, Spearman, AUPRC
    tier_stability.csv              - per tier x arm set comparison
    long_protein_changes_<sp>.csv   - per long protein: percentile / tier changes
    family_recall_arms.csv          - family recall per tier x arm
    recovered_false_negatives.csv   - baseline moderate-tier misses recovered
    summary.txt
"""

from pathlib import Path
from urllib.parse import quote
import urllib.request
import sys
import re

import numpy as np
import pandas as pd
from scipy import stats as sp_stats
from sklearn.metrics import average_precision_score, roc_auc_score

SPECIES = {
    "arabidopsis_thaliana": ("Arabidopsis thaliana", 3702),
    "oryza_sativa":         ("Oryza sativa Japonica Group", 39947),
    "vitis_vinifera":       ("Vitis vinifera", 29760),
}

# selection rules - identical to step 03 / config.yaml
TIERS = {"strict": 99.5, "moderate": 99.0, "lenient": 97.0}
TOP_N = 50

GO_TERMS = ["go:0006952", "go:0002376"]
FAMILIES = {
    "NLR (NB-ARC)":              ("or",  ["IPR002182"]),
    "PR proteins":               ("or",  ["IPR000726", "IPR000490", "IPR001938", "IPR014044"]),
    "RLK (LRR receptor kinase)": ("and", ["IPR001611", "IPR000719"]),
}

RESULTS = Path("results")
OUTDIR = RESULTS / "14_truncation_sensitivity"
OUTDIR.mkdir(parents=True, exist_ok=True)
BENCH_CACHE = RESULTS / "09_benchmark_curated" / "_cache"

REFSEQ_RE = re.compile(r"\b([NXY]P_\d+)")
strip_ver = lambda s: re.sub(r"\.\d+$", "", str(s))

ARMS = ["truncated", "windowed_mean", "windowed_max"]


# ── curated set retrieval (reuses the step 09/10 caches) ────────────────────

def _uniprot_tsv(url: str, cache_file: Path, log) -> str:
    if cache_file.exists():
        return cache_file.read_text(encoding="utf-8")
    log(f"    UniProt query (not cached): {cache_file.name}")
    req = urllib.request.Request(url, headers={"User-Agent": "PlantDefenseESM-benchmark"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        text = resp.read().decode("utf-8")
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    cache_file.write_text(text, encoding="utf-8")
    return text


def _refseq_ids(text: str) -> set:
    ids = set()
    for line in text.splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        for m in REFSEQ_RE.findall(parts[1]):
            ids.add(strip_ver(m))
    return ids


def fetch_curated_refseq(taxid, log) -> set:
    cache_file = BENCH_CACHE / f"uniprot_curated_{taxid}.tsv"
    go_q = " OR ".join(GO_TERMS)
    query = f"(organism_id:{taxid}) AND ({go_q})"
    url = ("https://rest.uniprot.org/uniprotkb/stream"
           f"?query={quote(query)}&fields=accession,xref_refseq&format=tsv")
    return _refseq_ids(_uniprot_tsv(url, cache_file, log))


def fetch_family_refseq(taxid, mode, interpros, log) -> set:
    tag = "_".join(interpros)
    cache_file = BENCH_CACHE / f"uniprot_{taxid}_{tag}.tsv"
    joiner = " AND " if mode == "and" else " OR "
    ip_q = joiner.join(f"xref:interpro-{ip}" for ip in interpros)
    query = f"(organism_id:{taxid}) AND ({ip_q})"
    url = ("https://rest.uniprot.org/uniprotkb/stream"
           f"?query={quote(query)}&fields=accession,xref_refseq&format=tsv")
    return _refseq_ids(_uniprot_tsv(url, cache_file, log))


# ── embedding / scoring helpers ─────────────────────────────────────────────

def load_npz(path):
    d = np.load(path, allow_pickle=True)
    return d, [str(x) for x in d["protein_ids"]]


def l2norm(mat):
    n = np.linalg.norm(mat, axis=1, keepdims=True)
    return mat / np.where(n == 0, 1.0, n)


def window_mean(emb, prot_ids, weights):
    """Length-weighted mean of window embeddings per protein.
    Windows are all `max_seq_len` long, so this reduces to the unweighted
    mean; the weighting is kept for generality."""
    df = pd.DataFrame({"pid": prot_ids, "w": weights.astype(np.float64)})
    out = {}
    for pid, idx in df.groupby("pid").groups.items():
        idx = np.asarray(idx)
        w = df.loc[idx, "w"].values
        out[str(pid)] = np.average(emb[idx], axis=0, weights=w)
    return out


def build_centroids(anchor_emb, anchor_ids, id2cat, override=None):
    """L2-normalized mean of the anchor embeddings per category.
    `override` maps anchor_id -> replacement embedding (windowed anchors)."""
    cats = sorted(set(id2cat[a] for a in anchor_ids if a in id2cat))
    cents = {}
    for cat in cats:
        vecs = []
        for i, aid in enumerate(anchor_ids):
            if id2cat.get(aid) != cat:
                continue
            v = override[aid] if (override and aid in override) else anchor_emb[i]
            vecs.append(v)
        c = np.mean(vecs, axis=0)
        cents[cat] = c / np.linalg.norm(c)
    return cats, np.stack([cents[c] for c in cats])


def select(sim, cats):
    """Reproduce step 03: per-category Z and percentile, tier flags."""
    N = sim.shape[0]
    z = sp_stats.zscore(sim, axis=0)
    z_max = z.max(axis=1)
    best = np.array(cats)[z.argmax(axis=1)]
    pctl = np.apply_along_axis(lambda col: sp_stats.rankdata(col) / N * 100, 0, sim)
    pctl_max = pctl.max(axis=1)

    top_union = np.zeros(N, dtype=bool)
    for j in range(sim.shape[1]):
        top_union[np.argpartition(sim[:, j], -TOP_N)[-TOP_N:]] = True

    flags = {t: (pctl_max >= p) | top_union for t, p in TIERS.items()}
    return {"z_max": z_max, "pctl_max": pctl_max, "best": best, "flags": flags}


def set_compare(base_mask, new_mask):
    b, n = set(np.where(base_mask)[0]), set(np.where(new_mask)[0])
    inter, union = b & n, b | n
    return {
        "n_baseline": len(b), "n_arm": len(n),
        "jaccard": round(len(inter) / len(union), 4) if union else 1.0,
        "retained_frac": round(len(inter) / len(b), 4) if b else 1.0,
        "added": len(n - b), "removed": len(b - n),
    }


# ── main ────────────────────────────────────────────────────────────────────

def main():
    lines = []
    def log(s=""):
        print(s); lines.append(s)

    centroid_rows, arm_rows, tier_rows, fam_rows, fn_rows = [], [], [], [], []

    log("=" * 74)
    log("Truncation sensitivity: truncated vs windowed embeddings")
    log("=" * 74)

    for folder, (display, taxid) in SPECIES.items():
        emb_dir = RESULTS / folder / "02_embed_proteins"
        win_dir = RESULTS / folder / "13_windowed_embed"
        meta_p = RESULTS / folder / "01_fetch_anchors" / "anchors_metadata.csv"
        need = [emb_dir / "anchor_embeddings.npz", emb_dir / "proteome_embeddings.npz",
                win_dir / "windowed_proteome_embeddings.npz",
                win_dir / "windowed_anchor_embeddings.npz", meta_p]
        if not all(p.exists() for p in need):
            log(f"\n[SKIP] {display}: run steps 02 and 13 first")
            continue

        log(f"\n{'=' * 74}")
        log(f"### {display}")
        log("=" * 74)

        # ---- embeddings -------------------------------------------------
        a_d, a_ids = load_npz(emb_dir / "anchor_embeddings.npz")
        p_d, p_ids = load_npz(emb_dir / "proteome_embeddings.npz")
        anchor_emb, proteome_emb = a_d["embeddings"], p_d["embeddings"]

        wa_d, wa_prot = load_npz(win_dir / "windowed_anchor_embeddings.npz")
        wp_d, wp_prot = load_npz(win_dir / "windowed_proteome_embeddings.npz")

        meta = pd.read_csv(meta_p)
        id2cat = dict(zip(meta["uniprot_id"].astype(str), meta["category"]))

        a_wmean = window_mean(wa_d["embeddings"], wa_prot,
                              wa_d["win_end"] - wa_d["win_start"])
        p_wmean = window_mean(wp_d["embeddings"], wp_prot,
                              wp_d["win_end"] - wp_d["win_start"])

        row_of = {pid: i for i, pid in enumerate(p_ids)}
        long_ids = [pid for pid in p_wmean if pid in row_of]
        long_idx = np.array([row_of[p] for p in long_ids])
        log(f"  proteome {len(p_ids):,} | long (windowed) {len(long_idx):,} "
            f"({len(long_idx)/len(p_ids):.2%}) | anchors windowed {len(a_wmean)}")
        if len(long_ids) != len(p_wmean):
            log(f"  [WARN] {len(p_wmean) - len(long_ids)} windowed proteins not "
                f"found in the step-02 proteome; check for a stale step-13 run")

        # ---- centroids --------------------------------------------------
        cats, C_base = build_centroids(anchor_emb, a_ids, id2cat)
        _, C_win = build_centroids(anchor_emb, a_ids, id2cat, override=a_wmean)
        log("\n  Centroid shift (cosine, baseline vs windowed anchors):")
        for j, cat in enumerate(cats):
            cs = float(C_base[j] @ C_win[j])
            n_trunc = sum(1 for a in a_wmean if id2cat.get(a) == cat)
            log(f"    {cat:20s} {cs:.5f}   ({n_trunc} of its anchors windowed)")
            centroid_rows.append({"species": display, "category": cat,
                                  "cosine_base_vs_windowed": round(cs, 5),
                                  "n_anchors_windowed": n_trunc})

        # ---- three scoring arms -----------------------------------------
        P_norm = l2norm(proteome_emb)
        sims = {"truncated": P_norm @ C_base.T}

        patched = proteome_emb.copy()
        for pid, i in zip(long_ids, long_idx):
            patched[i] = p_wmean[pid]
        sims["windowed_mean"] = l2norm(patched) @ C_win.T

        sim_max = (P_norm @ C_win.T).copy()
        wsim = l2norm(wp_d["embeddings"]) @ C_win.T
        wdf = pd.DataFrame(wsim, columns=cats)
        wdf["pid"] = wp_prot
        gmax = wdf.groupby("pid")[cats].max()
        keep = [p for p in long_ids if p in gmax.index]
        sim_max[np.array([row_of[p] for p in keep])] = gmax.loc[keep, cats].values
        sims["windowed_max"] = sim_max

        res = {arm: select(sims[arm], cats) for arm in ARMS}
        base = res["truncated"]

        # sanity check against the published pipeline output
        val_p = RESULTS / folder / "04_validate_annotations" / "validated_results.csv"
        val = None
        if val_p.exists():
            val = pd.read_csv(val_p)
            if "protein_id" not in val.columns:
                val = val.rename(columns={val.columns[0]: "protein_id"})
            if "defense_moderate" in val.columns:
                log(f"\n  baseline moderate recomputed: "
                    f"{int(base['flags']['moderate'].sum()):,}  "
                    f"(pipeline: {int(val['defense_moderate'].sum()):,})")

        # ---- stability ---------------------------------------------------
        log("\n  Set stability vs the truncated baseline:")
        for arm in ARMS[1:]:
            rho = sp_stats.spearmanr(base["z_max"], res[arm]["z_max"]).correlation
            rho_long = sp_stats.spearmanr(base["z_max"][long_idx],
                                          res[arm]["z_max"][long_idx]).correlation
            for tier in TIERS:
                cmp = set_compare(base["flags"][tier], res[arm]["flags"][tier])
                tier_rows.append({"species": display, "arm": arm, "tier": tier, **cmp})
                if tier == "moderate":
                    log(f"    {arm:14s} moderate: {cmp['n_arm']:,} candidates, "
                        f"Jaccard {cmp['jaccard']:.3f}, retained {cmp['retained_frac']:.3f}, "
                        f"+{cmp['added']} / -{cmp['removed']}")
            log(f"    {arm:14s} Spearman(z_max) all {rho:.4f} | "
                f"long proteins only {rho_long:.4f}")
            arm_rows.append({"species": display, "arm": arm,
                             "spearman_all": round(float(rho), 4),
                             "spearman_long_only": round(float(rho_long), 4),
                             "n_long": len(long_idx)})

        # ---- curated benchmark (AUPRC / ROC-AUC) -------------------------
        curated = fetch_curated_refseq(taxid, log)
        pid_base = pd.Series(p_ids).map(strip_ver)
        y = pid_base.isin(curated).astype(int).values
        n_pos = int(y.sum())
        log(f"\n  curated GO positives mapped: {n_pos:,} "
            f"({n_pos/len(p_ids):.2%} prevalence)")
        if n_pos >= 20:
            for arm in ARMS:
                ap = average_precision_score(y, res[arm]["z_max"])
                ra = roc_auc_score(y, res[arm]["z_max"])
                m = res[arm]["flags"]["moderate"]
                tp = int((m & (y == 1)).sum())
                prec, rec = tp / m.sum(), tp / n_pos
                log(f"    {arm:14s} AUPRC {ap:.4f}  ROC-AUC {ra:.4f}  "
                    f"moderate precision {prec:.3f} recall {rec:.3f}")
                arm_rows.append({"species": display, "arm": arm,
                                 "AUPRC": round(float(ap), 4),
                                 "ROC_AUC": round(float(ra), 4),
                                 "moderate_precision": round(prec, 4),
                                 "moderate_recall": round(rec, 4),
                                 "n_curated_positives": n_pos})
        else:
            log("    [skip] too few curated positives for benchmark metrics")

        # ---- family recall + recovered false negatives --------------------
        is_long = np.zeros(len(p_ids), dtype=bool)
        is_long[long_idx] = True
        for fam, (mode, ips) in FAMILIES.items():
            curated_fam = fetch_family_refseq(taxid, mode, ips, log)
            in_fam = pid_base.isin(curated_fam).values
            n_fam = int(in_fam.sum())
            if n_fam == 0:
                continue
            log(f"\n  {fam}: {n_fam} in proteome ({int((in_fam & is_long).sum())} >1,022 aa)")
            for tier in TIERS:
                cells = {}
                for arm in ARMS:
                    r = int((in_fam & res[arm]["flags"][tier]).sum())
                    cells[arm] = r / n_fam
                    fam_rows.append({"species": display, "family": fam, "tier": tier,
                                     "arm": arm, "n_curated_in_proteome": n_fam,
                                     "recovered": r, "recall": round(r / n_fam, 4)})
                log(f"    recall@{tier:9s} " + "  ".join(
                    f"{a}={cells[a]:.3f}" for a in ARMS))

            # baseline moderate-tier misses, and their fate under windowing
            miss = in_fam & ~base["flags"]["moderate"]
            n_miss = int(miss.sum())
            row = {"species": display, "family": fam,
                   "n_curated_in_proteome": n_fam,
                   "n_missed_baseline_moderate": n_miss,
                   "n_missed_and_long": int((miss & is_long).sum())}
            for arm in ARMS[1:]:
                rec_mask = miss & res[arm]["flags"]["moderate"]
                row[f"recovered_by_{arm}"] = int(rec_mask.sum())
                row[f"recovered_by_{arm}_long"] = int((rec_mask & is_long).sum())
                row[f"newly_lost_{arm}"] = int(
                    (in_fam & base["flags"]["moderate"] & ~res[arm]["flags"]["moderate"]).sum())
            fn_rows.append(row)
            log(f"    baseline moderate misses: {n_miss} "
                f"({row['n_missed_and_long']} of them >1,022 aa)  -> recovered: "
                + ", ".join(f"{a} {row[f'recovered_by_{a}']}" for a in ARMS[1:]))

        # ---- per-long-protein detail -------------------------------------
        detail = pd.DataFrame({
            "protein_id": [p_ids[i] for i in long_idx],
            "is_curated_GO": y[long_idx].astype(bool),
        })
        for arm in ARMS:
            detail[f"pctl_{arm}"] = res[arm]["pctl_max"][long_idx].round(3)
            detail[f"z_{arm}"] = res[arm]["z_max"][long_idx].round(4)
            detail[f"moderate_{arm}"] = res[arm]["flags"]["moderate"][long_idx]
            detail[f"category_{arm}"] = res[arm]["best"][long_idx]
        if val is not None and "description" in val.columns:
            detail = detail.merge(
                val[["protein_id", "description"]], on="protein_id", how="left")
        detail["pctl_delta_mean"] = (detail["pctl_windowed_mean"]
                                     - detail["pctl_truncated"]).round(3)
        detail["pctl_delta_max"] = (detail["pctl_windowed_max"]
                                    - detail["pctl_truncated"]).round(3)
        detail.sort_values("pctl_delta_max", ascending=False).to_csv(
            OUTDIR / f"long_protein_changes_{folder}.csv", index=False)

        for arm in ARMS[1:]:
            d = detail[f"pctl_delta_{arm.split('_')[1]}"]
            gained = int((~detail["moderate_truncated"] & detail[f"moderate_{arm}"]).sum())
            lost = int((detail["moderate_truncated"] & ~detail[f"moderate_{arm}"]).sum())
            log(f"\n  long proteins under {arm}: median percentile shift "
                f"{d.median():+.3f} (IQR {d.quantile(.25):+.3f} to {d.quantile(.75):+.3f}); "
                f"{gained} crossed into the moderate tier, {lost} dropped out")

    # ---- write tables -----------------------------------------------------
    for rows, name in [(centroid_rows, "centroid_shift.csv"),
                       (arm_rows, "arm_summary.csv"),
                       (tier_rows, "tier_stability.csv"),
                       (fam_rows, "family_recall_arms.csv"),
                       (fn_rows, "recovered_false_negatives.csv")]:
        if rows:
            pd.DataFrame(rows).to_csv(OUTDIR / name, index=False)
    (OUTDIR / "summary.txt").write_text("\n".join(lines), encoding="utf-8")
    log(f"\nOutputs -> {OUTDIR}")
    log("Done OK")


if __name__ == "__main__":
    if not RESULTS.exists():
        sys.exit("ERROR: run from the project root (no 'results/' folder here).")
    main()
