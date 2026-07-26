#!/usr/bin/env python3
"""
r2_novel_domain_support.py
==========================
Orthogonal (non-keyword) support for the keyword-negative candidates
(Reviewer 2: the ">50% novel" claim rests only on keyword absence;
Reviewer 3: tone down "novel defense genes").

Runs from cached step-04 outputs + UniProt InterPro lookups. No GPU, no
re-embedding. Reuses the same InterPro retrieval approach as
10_benchmark_families.py.

Rationale
---------
A moderate-tier candidate is called keyword-negative when its RefSeq DESCRIPTION
carries none of the defense keywords. That is a statement about the text
annotation, not about the protein. As an independent line of evidence, we
ask how many of these keyword-negative candidates nonetheless carry a
recognized defense-related INTERPRO DOMAIN. Domain assignments come from
InterPro/Pfam signatures and are independent of (a) the ESM-2 embeddings
used for prediction and (b) the RefSeq description keywords used to define
annotation support. A keyword-negative candidate that still carries an NB-ARC, LRR,
kinase, chitinase, etc. domain is supported by orthogonal evidence rather
than being a threshold artefact.

For each species it reports, at the moderate tier:
    - n keyword-negative candidates
    - how many carry >=1 defense-related InterPro domain (orthogonal support)
    - the same for keyword-positive candidates (sanity reference)
    - a per-domain breakdown among the keyword-negative set

Run from the project root (the folder that contains `results/`):
    python3 r2_novel_domain_support.py
    # or custom dirs:
    python3 r2_novel_domain_support.py results/arabidopsis_thaliana

Outputs (written next to this script):
    r2_novel_domain_support_<species>.csv   - per-domain counts (keyword-neg set)
    r2_novel_domain_support_summary.csv     - combined per-species summary
"""

import sys
import re
import urllib.request
from urllib.parse import quote
from pathlib import Path

import pandas as pd

# species folder -> (display name, UniProt taxonomy id)   [same as step 10]
SPECIES = {
    "arabidopsis_thaliana": ("Arabidopsis thaliana", 3702),
    "oryza_sativa":         ("Oryza sativa Japonica Group", 39947),
    "vitis_vinifera":       ("Vitis vinifera", 29760),
}

# Defense-related InterPro domains used as orthogonal evidence.
# (Superset of the family-benchmark signatures in 10_benchmark_families.py.)
DEFENSE_DOMAINS = {
    "IPR002182": "NB-ARC (NLR)",
    "IPR032675": "LRR domain",
    "IPR001611": "LRR (receptor)",
    "IPR000719": "Protein kinase",
    "IPR011009": "Kinase-like",
    "IPR000726": "GH19 chitinase (PR-3)",
    "IPR000490": "GH17 glucanase (PR-2)",
    "IPR001938": "Thaumatin (PR-5)",
    "IPR014044": "CAP / PR-1",
    "IPR008680": "NPR1-like",
    "IPR003657": "WRKY domain",
    "IPR001128": "Cytochrome P450",
    "IPR044851": "PAL / histidase",
}

RESULTS = Path("results")
REFSEQ_RE = re.compile(r"\b([NXY]P_\d+)")
strip_ver = lambda s: re.sub(r"\.\d+$", "", str(s))


def fetch_domain_refseq(taxid, interpro, cache, log) -> set:
    """Versionless RefSeq accessions in this organism carrying one InterPro domain."""
    cache_file = cache / f"uniprot_{taxid}_{interpro}.tsv"
    if cache_file.exists():
        text = cache_file.read_text(encoding="utf-8")
    else:
        query = f"(organism_id:{taxid}) AND (xref:interpro-{interpro})"
        url = (
            "https://rest.uniprot.org/uniprotkb/stream"
            f"?query={quote(query)}&fields=accession,xref_refseq&format=tsv"
        )
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
    cache = Path(__file__).resolve().parent / "_uniprot_domain_cache"
    cache.mkdir(exist_ok=True)
    here = Path(__file__).resolve().parent

    lines = []
    def log(s=""):
        print(s); lines.append(s)

    dirs = sys.argv[1:] if len(sys.argv) > 1 else [
        f"results/{k}" for k in SPECIES]

    summary_rows = []
    log("=" * 70)
    log("Orthogonal domain support for keyword-negative candidates")
    log("=" * 70)

    for d in dirs:
        folder = Path(d).name
        if folder not in SPECIES:
            log(f"\n[skip] {d}: unknown species folder")
            continue
        display, taxid = SPECIES[folder]
        val = RESULTS / folder / "04_validate_annotations" / "validated_results.csv"
        if not val.exists():
            log(f"\n[skip] {val} not found")
            continue

        log(f"\n### {display}  (taxid {taxid})")
        df = pd.read_csv(val)
        if "protein_id" not in df.columns:
            df = df.rename(columns={df.columns[0]: "protein_id"})
        base = df["protein_id"].map(strip_ver)

        mod = df["defense_moderate"].astype(bool)
        kwneg = mod & (~df["has_defense_keyword"].astype(bool))
        kwpos = mod & (df["has_defense_keyword"].astype(bool))
        n_neg = int(kwneg.sum())
        n_pos = int(kwpos.sum())
        log(f"  moderate candidates: keyword-negative={n_neg:,}  "
            f"keyword-positive={n_pos:,}")

        # Build the union of all defense-domain-carrying RefSeq ids
        any_domain_ids = set()
        per_domain_rows = []
        neg_ids = set(base[kwneg])
        for ip, name in DEFENSE_DOMAINS.items():
            dom_ids = fetch_domain_refseq(taxid, ip, cache, log)
            any_domain_ids |= dom_ids
            n_neg_dom = len(neg_ids & dom_ids)
            per_domain_rows.append({
                "species": display, "interpro": ip, "domain": name,
                "kwneg_with_domain": n_neg_dom,
            })

        n_neg_supported = len(neg_ids & any_domain_ids)
        pos_ids = set(base[kwpos])
        n_pos_supported = len(pos_ids & any_domain_ids)

        frac_neg = n_neg_supported / n_neg if n_neg else 0.0
        frac_pos = n_pos_supported / n_pos if n_pos else 0.0

        log(f"  keyword-negative with >=1 defense InterPro domain: "
            f"{n_neg_supported:,} / {n_neg:,} ({frac_neg:.1%})")
        log(f"  keyword-positive with >=1 defense InterPro domain: "
            f"{n_pos_supported:,} / {n_pos:,} ({frac_pos:.1%})  [reference]")

        pd.DataFrame(per_domain_rows).sort_values(
            "kwneg_with_domain", ascending=False).to_csv(
            here / f"r2_novel_domain_support_{folder}.csv", index=False)

        summary_rows.append({
            "species": display,
            "moderate_kwneg": n_neg,
            "kwneg_with_defense_domain": n_neg_supported,
            "kwneg_domain_supported_frac": round(frac_neg, 4),
            "moderate_kwpos": n_pos,
            "kwpos_with_defense_domain": n_pos_supported,
            "kwpos_domain_supported_frac": round(frac_pos, 4),
        })

    if summary_rows:
        out = here / "r2_novel_domain_support_summary.csv"
        pd.DataFrame(summary_rows).to_csv(out, index=False)
        log(f"\nWrote combined summary -> {out}")
        log("Per-species tables -> r2_novel_domain_support_<species>.csv")


if __name__ == "__main__":
    if not RESULTS.exists():
        sys.exit("ERROR: run from the project root (no 'results/' folder here).")
    main()
