#!/usr/bin/env python3
"""
r3_false_negatives.py
=====================
Concrete false-negative examples for the Arabidopsis benchmark
(Reviewer 3, Major comment 1: "a careful Arabidopsis-only benchmark with
precision, recall, and false-negative examples would make the paper much
stronger").

Runs from cached step-04 outputs + the same curated-set retrieval used in
09_benchmark_curated.py / 10_benchmark_families.py. No GPU, no re-embedding.

A FALSE NEGATIVE here = a curated defense protein that IS present in the
A. thaliana RefSeq proteome but was NOT selected at the moderate tier
(defense_moderate == False). For each missed protein we report its
max Z-score, best category, percentile, and RefSeq description, so the
manuscript can give named examples and explain why they scored low
(e.g., truncation of long LRR receptors, divergent sequences, categories
under-represented in the anchor set).

Two curated sources (both independent of the ESM-2 embeddings and the
RefSeq description keywords):
  - GO defense set: GO:0006952 (defense response) + GO:0002376 (immune
    system process), retrieved from UniProt and mapped to RefSeq.
  - InterPro structural families: NLR (NB-ARC), PR proteins, LRR-RK.

Run from the project root (the folder that contains `results/`):
    python3 r3_false_negatives.py

Outputs (written next to this script):
    r3_false_negatives_GO.csv        - all missed curated GO proteins (ranked)
    r3_false_negatives_families.csv  - missed members of each structural family
    r3_false_negatives_summary.csv   - counts: present / recovered / missed
"""

import sys
import re
import urllib.request
from urllib.parse import quote
from pathlib import Path

import pandas as pd

FOLDER = "arabidopsis_thaliana"
DISPLAY = "Arabidopsis thaliana"
TAXID = 3702

GO_TERMS = ["go:0006952", "go:0002376"]
FAMILIES = {
    "NLR (NB-ARC)":              ("or",  ["IPR002182"]),
    "PR proteins":               ("or",  ["IPR000726", "IPR000490",
                                          "IPR001938", "IPR014044"]),
    "RLK (LRR receptor kinase)": ("and", ["IPR001611", "IPR000719"]),
}

RESULTS = Path("results")
REFSEQ_RE = re.compile(r"\b([NXY]P_\d+)")
strip_ver = lambda s: re.sub(r"\.\d+$", "", str(s))


def uniprot_refseq(query: str, cache_file: Path, log) -> set:
    if cache_file.exists():
        text = cache_file.read_text(encoding="utf-8")
    else:
        url = ("https://rest.uniprot.org/uniprotkb/stream"
               f"?query={quote(query)}&fields=accession,xref_refseq&format=tsv")
        log(f"    UniProt: {query}")
        req = urllib.request.Request(
            url, headers={"User-Agent": "PlantDefenseESM-benchmark"})
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
    cache = Path(__file__).resolve().parent / "_uniprot_fn_cache"
    cache.mkdir(exist_ok=True)
    here = Path(__file__).resolve().parent
    lines = []
    def log(s=""):
        print(s); lines.append(s)

    val = RESULTS / FOLDER / "04_validate_annotations" / "validated_results.csv"
    if not val.exists():
        sys.exit(f"ERROR: {val} not found")

    log("=" * 70)
    log(f"False-negative examples  ({DISPLAY}, taxid {TAXID})")
    log("=" * 70)

    df = pd.read_csv(val)
    if "protein_id" not in df.columns:
        df = df.rename(columns={df.columns[0]: "protein_id"})
    df["base"] = df["protein_id"].map(strip_ver)
    base_to_row = df.set_index("base")
    mod = df.set_index("base")["defense_moderate"].astype(bool)

    cols_out = ["protein_id", "z_max", "best_category", "description"]
    if "pctl_max" in df.columns:
        cols_out.insert(3, "pctl_max")

    summary_rows = []

    # ---- GO defense set ----
    go_q = f"(organism_id:{TAXID}) AND (" + " OR ".join(GO_TERMS) + ")"
    go_ids = uniprot_refseq(go_q, cache / f"go_{TAXID}.tsv", log)
    present = go_ids & set(df["base"])
    missed = [b for b in present if b in mod.index and not mod.loc[b]]
    log(f"\nGO defense set: curated={len(go_ids):,}  present_in_proteome="
        f"{len(present):,}  recovered={len(present)-len(missed):,}  "
        f"missed(false neg)={len(missed):,}")

    go_miss = base_to_row.loc[missed, cols_out].copy()
    go_miss = go_miss.sort_values("z_max", ascending=False)
    go_miss.to_csv(here / "r3_false_negatives_GO.csv", index=False)

    log("\n  Highest-scoring missed GO proteins (just below threshold):")
    for _, r in go_miss.head(8).iterrows():
        d = str(r["description"])[:70]
        log(f"    {r['protein_id']:>16s}  z={r['z_max']:.2f}  "
            f"{r['best_category']:<18s} {d}")

    summary_rows.append({
        "set": "GO defense", "curated_in_proteome": len(present),
        "recovered_moderate": len(present) - len(missed),
        "missed_false_negative": len(missed),
        "recall": round((len(present) - len(missed)) / len(present), 4)
        if present else 0.0,
    })

    # ---- structural families ----
    fam_rows = []
    for fam, (mode, ips) in FAMILIES.items():
        joiner = " AND " if mode == "and" else " OR "
        q = (f"(organism_id:{TAXID}) AND ("
             + joiner.join(f"xref:interpro-{ip}" for ip in ips) + ")")
        ids = uniprot_refseq(q, cache / f"fam_{TAXID}_{'_'.join(ips)}.tsv", log)
        pres = ids & set(df["base"])
        miss = [b for b in pres if b in mod.index and not mod.loc[b]]
        log(f"\n{fam}: present={len(pres):,}  recovered={len(pres)-len(miss):,}  "
            f"missed={len(miss):,}")
        sub = base_to_row.loc[miss, cols_out].copy()
        sub.insert(0, "family", fam)
        sub = sub.sort_values("z_max", ascending=False)
        fam_rows.append(sub)
        for _, r in sub.head(5).iterrows():
            d = str(r["description"])[:65]
            log(f"    {r['protein_id']:>16s}  z={r['z_max']:.2f}  {d}")
        summary_rows.append({
            "set": fam, "curated_in_proteome": len(pres),
            "recovered_moderate": len(pres) - len(miss),
            "missed_false_negative": len(miss),
            "recall": round((len(pres) - len(miss)) / len(pres), 4)
            if pres else 0.0,
        })

    if fam_rows:
        pd.concat(fam_rows, ignore_index=True).to_csv(
            here / "r3_false_negatives_families.csv", index=False)

    pd.DataFrame(summary_rows).to_csv(
        here / "r3_false_negatives_summary.csv", index=False)
    log(f"\nWrote -> r3_false_negatives_GO.csv, "
        f"r3_false_negatives_families.csv, r3_false_negatives_summary.csv")


if __name__ == "__main__":
    if not RESULTS.exists():
        sys.exit("ERROR: run from the project root (no 'results/' folder here).")
    main()
