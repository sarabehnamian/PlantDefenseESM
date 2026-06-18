#!/usr/bin/env python3
"""
09_benchmark_curated.py
=======================
Benchmark the ESM-2 defense predictions against an INDEPENDENT, curated
positive set of known defense genes, and report precision / recall / F1,
AUPRC, ROC-AUC, and a precision-recall curve.  (Reviewer 1 Comment 3;
also addresses Reviewer 2 and Reviewer 3 benchmarking requests.)

Curated positive set
--------------------
Proteins annotated under Gene Ontology "defense response" (GO:0006952)
and "immune system process" (GO:0002376) are retrieved from UniProt for
each species and mapped to the proteome's RefSeq accessions.  GO
annotations are curated from literature / experiments and are independent
of (a) the ESM-2 sequence embeddings used for prediction and (b) the
RefSeq description keywords used in the earlier validation step.

NOTE ON PRECISION
-----------------
Curated defense-gene sets are known to be incomplete (this is the premise
of the study).  Precision against such a set is therefore a CONSERVATIVE
LOWER BOUND: many predicted candidates without a GO defense annotation may
still be genuine, currently-unannotated defense genes.  Recall (fraction of
KNOWN defense genes recovered) is the cleaner, directly interpretable metric.

Run from the project root:
    python 09_benchmark_curated.py

Outputs -> results/09_benchmark_curated/
    benchmark_metrics.csv          - precision/recall/F1/fold per species x tier
    benchmark_overall.csv          - AUPRC, ROC-AUC, set sizes per species
    pr_curve_<species>.png         - precision-recall curve with tier points
    summary.txt                    - human-readable summary (also printed)
"""

from pathlib import Path
from urllib.parse import quote
import urllib.request
import sys
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import average_precision_score, roc_auc_score, precision_recall_curve

# species folder -> (display name, NCBI taxonomy id used by UniProt)
SPECIES = {
    "arabidopsis_thaliana": ("Arabidopsis thaliana", 3702),
    "oryza_sativa":         ("Oryza sativa Japonica Group", 39947),
    "vitis_vinifera":       ("Vitis vinifera", 29760),
}

# Curated GO terms defining the positive set (UniProt search is ancestor-aware)
GO_TERMS = ["go:0006952", "go:0002376"]   # defense response; immune system process

RESULTS = Path("results")
OUTDIR = RESULTS / "09_benchmark_curated"
OUTDIR.mkdir(parents=True, exist_ok=True)

REFSEQ_RE = re.compile(r"\b([NXY]P_\d+)")        # NP_ / XP_ / YP_ accessions
strip_ver = lambda s: re.sub(r"\.\d+$", "", str(s))   # NP_199333.1 -> NP_199333


def fetch_curated_refseq(taxid: int, cache: Path, log) -> set:
    """Return set of versionless RefSeq protein accessions in the curated GO set."""
    cache_file = cache / f"uniprot_curated_{taxid}.tsv"
    if cache_file.exists():
        text = cache_file.read_text(encoding="utf-8")
    else:
        go_q = " OR ".join(GO_TERMS)
        query = f"(organism_id:{taxid}) AND ({go_q})"
        url = (
            "https://rest.uniprot.org/uniprotkb/stream"
            f"?query={quote(query)}"
            "&fields=accession,xref_refseq"
            "&format=tsv"
        )
        log(f"  querying UniProt: {query}")
        req = urllib.request.Request(url, headers={"User-Agent": "PlantDefenseESM-benchmark"})
        with urllib.request.urlopen(req, timeout=120) as resp:
            text = resp.read().decode("utf-8")
        cache_file.write_text(text, encoding="utf-8")

    ids = set()
    for line in text.splitlines()[1:]:           # skip header
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        for m in REFSEQ_RE.findall(parts[1]):
            ids.add(strip_ver(m))
    return ids


def metrics_at(mask_called, y_true):
    """precision/recall/F1 for a boolean candidate mask vs y_true (0/1 array)."""
    tp = int(((mask_called) & (y_true == 1)).sum())
    fp = int(((mask_called) & (y_true == 0)).sum())
    fn = int(((~mask_called) & (y_true == 1)).sum())
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    return tp, fp, fn, prec, rec, f1


def main():
    cache = OUTDIR / "_cache"
    cache.mkdir(exist_ok=True)
    lines = []
    def log(s=""):
        print(s); lines.append(s)

    metric_rows = []
    overall_rows = []

    log("=" * 70)
    log("Benchmark vs curated GO defense set  (defense response + immune system)")
    log("=" * 70)

    for folder, (display, taxid) in SPECIES.items():
        val = RESULTS / folder / "04_validate_annotations" / "validated_results.csv"
        if not val.exists():
            log(f"\n[SKIP] {val} not found")
            continue

        log(f"\n### {display}  (taxid {taxid})")
        df = pd.read_csv(val)
        if "protein_id" not in df.columns:
            df = df.rename(columns={df.columns[0]: "protein_id"})

        # curated positive set -> RefSeq ids present in this proteome
        curated = fetch_curated_refseq(taxid, cache, log)
        prot_base = df["protein_id"].map(strip_ver)
        y_true = prot_base.isin(curated).astype(int).values
        n_pos = int(y_true.sum())
        n_total = len(df)
        prevalence = n_pos / n_total

        log(f"  curated GO proteins (RefSeq, unique): {len(curated):,}")
        log(f"  mapped into proteome (positives)    : {n_pos:,} / {n_total:,} "
            f"({prevalence:.2%} prevalence)")

        if n_pos == 0:
            log("  [WARN] no curated positives mapped; skipping metrics")
            continue

        score = df["z_max"].values
        auprc = average_precision_score(y_true, score)
        rocauc = roc_auc_score(y_true, score)
        log(f"  AUPRC = {auprc:.3f}   (random baseline = prevalence = {prevalence:.3f})")
        log(f"  ROC-AUC = {rocauc:.3f}")

        overall_rows.append({
            "species": display, "taxid": taxid,
            "n_proteome": n_total, "n_curated_in_proteome": n_pos,
            "prevalence": round(prevalence, 4),
            "AUPRC": round(auprc, 3), "ROC_AUC": round(rocauc, 3),
        })

        # metrics at each stringency tier
        log(f"  {'tier':9s} {'called':>7s} {'TP':>6s} {'precision':>10s} "
            f"{'recall':>8s} {'F1':>6s} {'fold':>6s}")
        for tier in ("strict", "moderate", "lenient"):
            col = f"defense_{tier}"
            if col not in df.columns:
                continue
            mask = df[col].values.astype(bool)
            tp, fp, fn, prec, rec, f1 = metrics_at(mask, y_true)
            fold = (prec / prevalence) if prevalence else float("nan")
            log(f"  {tier:9s} {int(mask.sum()):>7d} {tp:>6d} {prec:>10.3f} "
                f"{rec:>8.3f} {f1:>6.3f} {fold:>6.2f}")
            metric_rows.append({
                "species": display, "tier": tier,
                "n_called": int(mask.sum()), "n_curated_positives": n_pos,
                "TP": tp, "FP": fp, "FN": fn,
                "precision": round(prec, 4), "recall": round(rec, 4),
                "F1": round(f1, 4), "fold_vs_prevalence": round(fold, 3),
            })

        # precision-recall curve + tier points
        prec_c, rec_c, _ = precision_recall_curve(y_true, score)
        plt.figure(figsize=(5, 4))
        plt.plot(rec_c, prec_c, lw=2, label=f"ESM-2 z_max (AUPRC={auprc:.2f})")
        plt.axhline(prevalence, ls="--", c="grey", lw=1, label=f"random ({prevalence:.2f})")
        for tier, mk in [("strict", "o"), ("moderate", "s"), ("lenient", "^")]:
            col = f"defense_{tier}"
            if col in df.columns:
                m = df[col].values.astype(bool)
                _, _, _, p, r, _ = metrics_at(m, y_true)
                plt.scatter([r], [p], marker=mk, s=60, zorder=5, label=tier)
        plt.xlabel("Recall (known defense genes recovered)")
        plt.ylabel("Precision")
        plt.title(f"{display}: ESM-2 vs curated GO defense set")
        plt.legend(fontsize=7, loc="upper right")
        plt.tight_layout()
        fig_path = OUTDIR / f"pr_curve_{folder}.png"
        plt.savefig(fig_path, dpi=200); plt.close()
        log(f"  PR curve -> {fig_path.name}")

    # write tables
    if metric_rows:
        pd.DataFrame(metric_rows).to_csv(OUTDIR / "benchmark_metrics.csv", index=False)
    if overall_rows:
        pd.DataFrame(overall_rows).to_csv(OUTDIR / "benchmark_overall.csv", index=False)
    (OUTDIR / "summary.txt").write_text("\n".join(lines), encoding="utf-8")

    log(f"\nOutputs -> {OUTDIR}")
    log("Done. ")


if __name__ == "__main__":
    if not RESULTS.exists():
        sys.exit("ERROR: run from the project root (no 'results/' folder here).")
    main()
