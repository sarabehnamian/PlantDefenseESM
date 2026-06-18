#!/usr/bin/env python3
"""
10_benchmark_families.py
========================
Family-specific recall against curated, reviewer-named defense gene sets
(NLR catalogs, PR proteins, RLK immune receptors), complementing the
overall precision/recall benchmark in 09_benchmark_curated.py.
(Reviewer 1 Comment 3; Reviewer 3 major comment on curated positive sets.)

Curated family sets
-------------------
Each family is defined by its diagnostic InterPro domain signature(s),
retrieved from UniProt per species and mapped to the proteome's RefSeq
accessions.  InterPro domain assignments are independent of the ESM-2
sequence embeddings used for prediction and of the RefSeq description
keywords used in the earlier validation step.

    NLR (NB-ARC)                 : IPR002182
    PR proteins                  : IPR000726 (GH19 chitinase, PR-3)
                                   IPR000490 (GH17 glucanase, PR-2)
                                   IPR001938 (thaumatin, PR-5)
                                   IPR014044 (CAP, PR-1)            [OR]
    RLK (LRR receptor kinase)    : IPR001611 (LRR) AND IPR000719 (kinase)

For each family we report RECALL = fraction of that family present in the
proteome that is recovered among the ESM-2 candidates at each tier.

Run from the project root:
    python 10_benchmark_families.py

Outputs -> results/09_benchmark_curated/
    families_recall.csv
    families_summary.txt
"""

from pathlib import Path
from urllib.parse import quote
import urllib.request
import sys
import re

import pandas as pd

SPECIES = {
    "arabidopsis_thaliana": ("Arabidopsis thaliana", 3702),
    "oryza_sativa":         ("Oryza sativa Japonica Group", 39947),
    "vitis_vinifera":       ("Vitis vinifera", 29760),
}

# family -> (combine mode, list of InterPro IDs)
FAMILIES = {
    "NLR (NB-ARC)":              ("or",  ["IPR002182"]),
    "PR proteins":               ("or",  ["IPR000726", "IPR000490", "IPR001938", "IPR014044"]),
    "RLK (LRR receptor kinase)": ("and", ["IPR001611", "IPR000719"]),
}

RESULTS = Path("results")
OUTDIR = RESULTS / "09_benchmark_curated"
OUTDIR.mkdir(parents=True, exist_ok=True)

REFSEQ_RE = re.compile(r"\b([NXY]P_\d+)")
strip_ver = lambda s: re.sub(r"\.\d+$", "", str(s))


def fetch_family_refseq(taxid, mode, interpros, cache, log):
    tag = "_".join(interpros)
    cache_file = cache / f"uniprot_{taxid}_{tag}.tsv"
    if cache_file.exists():
        text = cache_file.read_text(encoding="utf-8")
    else:
        joiner = " AND " if mode == "and" else " OR "
        ip_q = joiner.join(f"xref:interpro-{ip}" for ip in interpros)
        query = f"(organism_id:{taxid}) AND ({ip_q})"
        url = (
            "https://rest.uniprot.org/uniprotkb/stream"
            f"?query={quote(query)}&fields=accession,xref_refseq&format=tsv"
        )
        log(f"    UniProt: {query}")
        req = urllib.request.Request(url, headers={"User-Agent": "PlantDefenseESM-benchmark"})
        with urllib.request.urlopen(req, timeout=120) as resp:
            text = resp.read().decode("utf-8")
        cache_file.write_text(text, encoding="utf-8")

    ids = set()
    for line in text.splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        for m in REFSEQ_RE.findall(parts[1]):
            ids.add(strip_ver(m))
    return ids


def main():
    cache = OUTDIR / "_cache"
    cache.mkdir(exist_ok=True)
    lines = []
    def log(s=""):
        print(s); lines.append(s)

    rows = []
    log("=" * 70)
    log("Family-specific recall vs curated InterPro defense sets")
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
        prot_base = set(df["protein_id"].map(strip_ver))
        masks = {t: df[f"defense_{t}"].values.astype(bool)
                 for t in ("strict", "moderate", "lenient") if f"defense_{t}" in df.columns}
        base_series = df["protein_id"].map(strip_ver)

        for fam, (mode, ips) in FAMILIES.items():
            curated = fetch_family_refseq(taxid, mode, ips, cache, log)
            present = curated & prot_base
            n_present = len(present)
            log(f"  {fam}: curated={len(curated):,}  present_in_proteome={n_present:,}")
            if n_present == 0:
                rows.append({"species": display, "family": fam,
                             "curated_total": len(curated), "present_in_proteome": 0})
                continue
            in_present = base_series.isin(present).values
            rec = {}
            for t, m in masks.items():
                recovered = int((in_present & m).sum())
                rec[t] = recovered / n_present
                log(f"      recall@{t:8s} = {rec[t]:.3f}  ({recovered}/{n_present})")
            rows.append({
                "species": display, "family": fam,
                "curated_total": len(curated), "present_in_proteome": n_present,
                "recall_strict": round(rec.get("strict", float("nan")), 3),
                "recall_moderate": round(rec.get("moderate", float("nan")), 3),
                "recall_lenient": round(rec.get("lenient", float("nan")), 3),
            })

    if rows:
        pd.DataFrame(rows).to_csv(OUTDIR / "families_recall.csv", index=False)
    (OUTDIR / "families_summary.txt").write_text("\n".join(lines), encoding="utf-8")
    log(f"\nOutputs -> {OUTDIR}")


if __name__ == "__main__":
    if not RESULTS.exists():
        sys.exit("ERROR: run from the project root (no 'results/' folder here).")
    main()
