#!/usr/bin/env python3
"""
r2_novel_pfam_support.py
========================
Orthogonal (non-keyword, sequence-based) support for the keyword-negative
candidates (Reviewer 2: the ">50% novel" claim rests only on keyword
absence; Reviewer 3: tone down "novel defense genes").

Unlike the earlier UniProt cross-reference attempt (which barely mapped to
the rice/grapevine RefSeq proteomes), this scans the candidate protein
SEQUENCES DIRECTLY against the Pfam-A HMM library with HMMER (via pyhmmer,
pip-installable, no system install, no GPU). Domain assignments are
therefore measured on our own RefSeq proteins and are independent of
(a) the ESM-2 embeddings used for prediction and (b) the RefSeq
description keywords used to define novelty.

For each species, among the MODERATE-tier candidates it reports:
    - keyword-negative ("novel") candidates carrying >=1 defense-related
      Pfam domain (orthogonal support)
    - the same for keyword-positive candidates (reference)
    - a per-domain breakdown for the keyword-negative set

Requirements
------------
    pip install pyhmmer
    Pfam-A.hmm (auto-downloaded & pressed on first run if absent; ~1.5 GB)

Run from the project root (the folder that contains `results/`):
    python3 r2_novel_pfam_support.py
    # or custom dirs / Pfam path:
    python3 r2_novel_pfam_support.py --pfam D:/data/Pfam-A.hmm results/arabidopsis_thaliana

Outputs (written next to this script):
    r2_novel_pfam_support_<species>.csv   - per-domain counts (keyword-neg set)
    r2_novel_pfam_support_summary.csv     - combined per-species summary
"""

import sys
import gzip
import shutil
import argparse
import urllib.request
from pathlib import Path

import pandas as pd

try:
    import pyhmmer
except ImportError:
    sys.exit("ERROR: pyhmmer not installed. Run:  pip install pyhmmer")

SPECIES = {
    "arabidopsis_thaliana": "Arabidopsis thaliana",
    "oryza_sativa":         "Oryza sativa",
    "vitis_vinifera":       "Vitis vinifera",
}

PFAM_URL = ("https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/"
            "Pfam-A.hmm.gz")

# Defense-related Pfam accessions used as orthogonal evidence. Matched by
# the leading Pfam accession (version-stripped) of each profile name/acc.
DEFENSE_PFAM = {
    "PF00931": "NB-ARC (NLR)",
    "PF13306": "LRR (LRR_5)",
    "PF00560": "LRR_1",
    "PF13855": "LRR_8",
    "PF08263": "LRR_N (LRRNT_2)",
    "PF00069": "Protein kinase",
    "PF07714": "Tyr kinase (PK_Tyr_Ser-Thr)",
    "PF00187": "Chitin-binding (PR-4)",
    "PF00182": "GH19 chitinase (PR-3)",
    "PF00332": "GH17 glucanase (PR-2)",
    "PF00314": "Thaumatin (PR-5)",
    "PF00188": "CAP / PR-1 (SCP)",
    "PF00954": "S-locus glycoprotein",
    "PF01453": "B-lectin (D-mannose)",
    "PF00704": "GH18 chitinase",
    "PF03106": "WRKY DNA-binding",
    "PF00067": "Cytochrome P450",
    "PF00221": "PAL / histidase (Lyase_aromatic)",
    "PF01419": "Jacalin lectin",
    "PF00190": "Cupin / germin",
    "PF00141": "Peroxidase",
}


def strip_acc(acc) -> str:
    """PF00069.27 -> PF00069 ; bytes/str/None safe."""
    if acc is None:
        return ""
    if isinstance(acc, bytes):
        acc = acc.decode()
    return acc.split(".")[0] if acc else ""


def ensure_pfam(pfam_arg, log) -> Path:
    """Locate or download+press Pfam-A.hmm. Returns path to the .hmm file."""
    if pfam_arg:
        p = Path(pfam_arg)
        if not p.exists():
            sys.exit(f"ERROR: --pfam path not found: {p}")
        return p

    here = Path(__file__).resolve().parent
    hmm = here / "Pfam-A.hmm"
    if hmm.exists():
        log(f"  Using Pfam-A.hmm -> {hmm}")
        return hmm

    gz = here / "Pfam-A.hmm.gz"
    if not gz.exists():
        log(f"  Downloading Pfam-A.hmm.gz (~1.5 GB) from EBI ...")
        urllib.request.urlretrieve(PFAM_URL, gz)
    log("  Decompressing Pfam-A.hmm.gz ...")
    with gzip.open(gz, "rb") as fi, open(hmm, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    return hmm


def load_candidate_protein_to_domains(fasta_path: Path, pfam_path: Path,
                                       defense_accs: set, log):
    """
    hmmscan candidate sequences vs Pfam-A; return
    {protein_id: set(defense Pfam accessions hit}}.
    Only domain hits meeting the per-profile gathering threshold are kept.
    """
    alphabet = pyhmmer.easel.Alphabet.amino()

    # Read candidate sequences (digital) keyed by protein_id (header up to '|')
    seqs = []
    with pyhmmer.easel.SequenceFile(str(fasta_path), digital=True,
                                    alphabet=alphabet) as sf:
        for s in sf:
            name = s.name.decode() if isinstance(s.name, bytes) else s.name
            pid = name.split("|")[0]
            s.name = pid.encode()
            seqs.append(s)
    log(f"  candidate sequences read: {len(seqs):,}")

    prot2dom = {}
    n_profiles = 0
    with pyhmmer.plan7.HMMFile(str(pfam_path)) as hmm_file:
        # hmmsearch over all profiles against the candidate sequence block,
        # using each profile's gathering cutoff (Pfam standard).
        for hits in pyhmmer.hmmsearch(hmm_file, seqs,
                                      bit_cutoffs="gathering"):
            n_profiles += 1
            q = hits.query
            acc = strip_acc(q.accession or q.name)
            if acc not in defense_accs:
                continue
            for hit in hits:
                if not hit.included:
                    continue
                pid = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
                prot2dom.setdefault(pid, set()).add(acc)
    log(f"  Pfam profiles scanned: {n_profiles:,}")
    return prot2dom


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pfam", default=None,
                    help="Path to Pfam-A.hmm (else auto-download next to script)")
    ap.add_argument("dirs", nargs="*",
                    default=[f"results/{k}" for k in SPECIES])
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    lines = []
    def log(s=""):
        print(s); lines.append(s)

    log("=" * 70)
    log("Orthogonal Pfam-domain support for keyword-negative candidates")
    log("=" * 70)

    defense_accs = set(DEFENSE_PFAM)
    pfam_path = ensure_pfam(args.pfam, log)

    dirs = args.dirs if args.dirs else [f"results/{k}" for k in SPECIES]
    summary_rows = []

    for d in dirs:
        folder = Path(d).name
        if folder not in SPECIES:
            log(f"\n[skip] {d}: unknown species folder")
            continue
        display = SPECIES[folder]
        base = Path(d)
        val = base / "04_validate_annotations" / "validated_results.csv"
        fasta = base / "05_extract_candidates" / "candidate_sequences.fasta"
        if not val.exists() or not fasta.exists():
            log(f"\n[skip] {display}: missing validated_results.csv or "
                f"candidate_sequences.fasta")
            continue

        log(f"\n### {display}")
        df = pd.read_csv(val)
        if "protein_id" not in df.columns:
            df = df.rename(columns={df.columns[0]: "protein_id"})
        mod = df["defense_moderate"].astype(bool)
        kwneg_ids = set(df.loc[mod & ~df["has_defense_keyword"].astype(bool),
                               "protein_id"])
        kwpos_ids = set(df.loc[mod & df["has_defense_keyword"].astype(bool),
                               "protein_id"])
        log(f"  moderate: keyword-negative={len(kwneg_ids):,}  "
            f"keyword-positive={len(kwpos_ids):,}")

        prot2dom = load_candidate_protein_to_domains(
            fasta, pfam_path, defense_accs, log)

        neg_sup = {p for p in kwneg_ids if p in prot2dom}
        pos_sup = {p for p in kwpos_ids if p in prot2dom}
        frac_neg = len(neg_sup) / len(kwneg_ids) if kwneg_ids else 0.0
        frac_pos = len(pos_sup) / len(kwpos_ids) if kwpos_ids else 0.0

        log(f"  keyword-negative with >=1 defense Pfam domain: "
            f"{len(neg_sup):,} / {len(kwneg_ids):,} ({frac_neg:.1%})")
        log(f"  keyword-positive with >=1 defense Pfam domain: "
            f"{len(pos_sup):,} / {len(kwpos_ids):,} ({frac_pos:.1%})  [reference]")

        # per-domain breakdown among keyword-negative
        rows = []
        for acc, name in DEFENSE_PFAM.items():
            n = sum(1 for p in kwneg_ids if acc in prot2dom.get(p, set()))
            rows.append({"species": display, "pfam": acc, "domain": name,
                         "kwneg_with_domain": n})
        pd.DataFrame(rows).sort_values(
            "kwneg_with_domain", ascending=False).to_csv(
            here / f"r2_novel_pfam_support_{folder}.csv", index=False)

        summary_rows.append({
            "species": display,
            "moderate_kwneg": len(kwneg_ids),
            "kwneg_with_defense_domain": len(neg_sup),
            "kwneg_domain_supported_frac": round(frac_neg, 4),
            "moderate_kwpos": len(kwpos_ids),
            "kwpos_with_defense_domain": len(pos_sup),
            "kwpos_domain_supported_frac": round(frac_pos, 4),
        })

    if summary_rows:
        out = here / "r2_novel_pfam_support_summary.csv"
        pd.DataFrame(summary_rows).to_csv(out, index=False)
        log(f"\nWrote combined summary -> {out}")
        log("Per-species tables -> r2_novel_pfam_support_<species>.csv")


if __name__ == "__main__":
    if not Path("results").exists():
        sys.exit("ERROR: run from the project root (no 'results/' folder here).")
    main()
