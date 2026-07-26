#!/usr/bin/env python3
"""
12_benchmark_baselines.py
=========================
Head-to-head comparison of the ESM-2 pipeline against standard, non-embedding
baselines on the SAME independent curated positive set, to show whether the
ESM-2 enrichment reflects a real advantage or just generic embedding behaviour.
(Reviewer 2: lack of comparative benchmarking; Reviewer 3: method needs
baseline comparisons.)

Methods compared (all evaluated against the curated GO defense set used in
09_benchmark_curated.py, so the ground truth is identical for every method):

  1. ESM-2 (moderate tier)   - this pipeline (defense_moderate flag, ranked by z_max)
  2. BLAST-to-anchor         - blastp of the proteome vs the 33 anchor proteins;
                               per protein = best bitscore to any anchor.
                               Candidate set = top-N by bitscore, N matched to
                               the ESM-2 moderate candidate count (fair, equal-N).
  3. InterPro-domain         - proteins carrying any diagnostic defense domain
                               (NB-ARC / PR / LRR+kinase; same domains as Table 4).
  4. RefSeq-keyword          - proteins whose RefSeq description contains a
                               defense keyword (the has_defense_keyword flag).

Also reports a RECOVERED-BY breakdown among the curated positives:
how many curated defense genes are recovered by ESM-2, by BLAST, by domains,
and by which combinations (the "ESM-2 only / standard only / both" picture).

Focus species: Arabidopsis thaliana (the only species with a sufficiently
large curated set; rice and grapevine map <20 curated proteins). The loop
runs any species that has results, but baselines are only meaningful where
enough curated positives exist.

Requirements:
  - results/ produced by steps 00-04 (validated_results.csv, proteome.fasta)
  - results/01_fetch_anchors/anchors.fasta
  - BLAST+ (`makeblastdb`, `blastp`) on PATH for the BLAST baseline.
    If BLAST+ is not found, that one method is skipped and the rest still run.

Run from the project root:
    python 12_benchmark_baselines.py

Outputs -> results/12_benchmark_baselines/
    baseline_metrics.csv      - method x precision/recall/F1/AUPRC/fold
    recovered_breakdown.csv   - Venn-style counts among curated positives
    summary.txt               - human-readable summary (also printed)
"""

from pathlib import Path
from urllib.parse import quote
import urllib.request
import subprocess
import shutil
import tempfile
import sys
import re

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

# ── species -> (display name, UniProt/NCBI taxonomy id) ──────────────────────
SPECIES = {
    "arabidopsis_thaliana": ("Arabidopsis thaliana", 3702),
    "oryza_sativa":         ("Oryza sativa Japonica Group", 39947),
    "vitis_vinifera":       ("Vitis vinifera", 29760),
}

GO_TERMS = ["go:0006952", "go:0002376"]   # defense response; immune system process

# Diagnostic InterPro domains for the "domain-based" baseline (same as Table 4)
DOMAIN_FAMILIES = {
    "NLR (NB-ARC)":              ("or",  ["IPR002182"]),
    "PR proteins":               ("or",  ["IPR000726", "IPR000490",
                                          "IPR001938", "IPR014044"]),
    "RLK (LRR receptor kinase)": ("and", ["IPR001611", "IPR000719"]),
}

MIN_CURATED = 50          # don't compute baselines if fewer positives than this

RESULTS = Path("results")
OUTDIR = RESULTS / "12_benchmark_baselines"
OUTDIR.mkdir(parents=True, exist_ok=True)

REFSEQ_RE = re.compile(r"\b([NXY]P_\d+)")
strip_ver = lambda s: re.sub(r"\.\d+$", "", str(s))


# ── UniProt fetch helpers (mirror 09 / 10) ───────────────────────────────────
def _uniprot_stream(query: str, cache_file: Path, log) -> str:
    if cache_file.exists():
        return cache_file.read_text(encoding="utf-8")
    url = ("https://rest.uniprot.org/uniprotkb/stream"
           f"?query={quote(query)}&fields=accession,xref_refseq&format=tsv")
    log(f"  querying UniProt: {query}")
    req = urllib.request.Request(url, headers={"User-Agent": "PlantDefenseESM-baselines"})
    with urllib.request.urlopen(req, timeout=180) as resp:
        text = resp.read().decode("utf-8")
    cache_file.write_text(text, encoding="utf-8")
    return text


def _refseq_from_tsv(text: str) -> set:
    ids = set()
    for line in text.splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        for m in REFSEQ_RE.findall(parts[1]):
            ids.add(strip_ver(m))
    return ids


def fetch_curated_refseq(taxid: int, cache: Path, log) -> set:
    go_q = " OR ".join(GO_TERMS)
    text = _uniprot_stream(f"(organism_id:{taxid}) AND ({go_q})",
                           cache / f"uniprot_curated_{taxid}.tsv", log)
    return _refseq_from_tsv(text)


def fetch_domain_refseq(taxid: int, cache: Path, log) -> set:
    """Union of the diagnostic-domain families = domain-based predictor."""
    domain_ids = set()
    for fam, (mode, iprs) in DOMAIN_FAMILIES.items():
        if mode == "or":
            q = "(" + " OR ".join(f"xref:interpro-{ip}" for ip in iprs) + ")"
        else:  # and
            q = "(" + " AND ".join(f"xref:interpro-{ip}" for ip in iprs) + ")"
        text = _uniprot_stream(
            f"(organism_id:{taxid}) AND {q}",
            cache / f"uniprot_domain_{taxid}_{fam.split()[0]}.tsv", log)
        domain_ids |= _refseq_from_tsv(text)
    return domain_ids


# ── metric helper ────────────────────────────────────────────────────────────
def metrics(mask_called, y_true, prevalence):
    mask = np.asarray(mask_called, bool)
    tp = int((mask & (y_true == 1)).sum())
    fp = int((mask & (y_true == 0)).sum())
    fn = int((~mask & (y_true == 1)).sum())
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    fold = (prec / prevalence) if prevalence else float("nan")
    return dict(n_called=int(mask.sum()), TP=tp, FP=fp, FN=fn,
                precision=round(prec, 4), recall=round(rec, 4),
                F1=round(f1, 4), fold_vs_background=round(fold, 3))


# ── BLAST-to-anchor baseline ─────────────────────────────────────────────────
def blast_best_bitscore(proteome_fasta: Path, anchor_fasta: Path, log) -> dict:
    """Return {versionless_protein_id: best_bitscore_to_any_anchor}."""
    if not (shutil.which("makeblastdb") and shutil.which("blastp")):
        log("  [BLAST] makeblastdb/blastp not found on PATH - skipping BLAST baseline")
        return None
    if not proteome_fasta.exists() or not anchor_fasta.exists():
        log(f"  [BLAST] missing FASTA ({proteome_fasta} / {anchor_fasta}) - skipping")
        return None

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        db = tmp / "anchors"
        log("  [BLAST] building anchor database ...")
        subprocess.run(["makeblastdb", "-in", str(anchor_fasta),
                        "-dbtype", "prot", "-out", str(db)],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        out_tsv = tmp / "hits.tsv"
        log("  [BLAST] running blastp (proteome vs anchors) ...")
        subprocess.run(
            ["blastp", "-query", str(proteome_fasta), "-db", str(db),
             "-evalue", "10", "-max_target_seqs", "5", "-num_threads", "4",
             "-outfmt", "6 qseqid sseqid bitscore evalue", "-out", str(out_tsv)],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        best = {}
        with open(out_tsv) as fh:
            for line in fh:
                q, s, bits, ev = line.rstrip("\n").split("\t")
                q = strip_ver(q.split()[0])
                b = float(bits)
                if b > best.get(q, -1):
                    best[q] = b
    log(f"  [BLAST] proteins with >=1 anchor hit: {len(best):,}")
    return best


# ── pure-Python homology baseline (no BLAST binary needed) ───────────────────
# Smith-Waterman local alignment (Biopython, BLOSUM62) to the anchor proteins,
# with a k-mer seed prefilter so only proteins that share words with an anchor
# are aligned -- the same seed-then-extend logic BLAST uses. Reported as
# "alignment-to-anchor", an alignment/homology baseline (not NCBI BLAST).
_VALID_AA = set("ACDEFGHIKLMNPQRSTVWYXBZ*")
_KMER = 5            # seed word length
_MIN_SHARED = 4     # min shared k-mers with an anchor to bother aligning


def _clean(seq: str) -> str:
    s = str(seq).upper()
    return "".join(c if c in _VALID_AA else "X" for c in s)


def _kmers(seq: str, k: int = _KMER):
    return {seq[i:i + k] for i in range(len(seq) - k + 1)} if len(seq) >= k else set()


def alignment_best_score(proteome_fasta: Path, anchor_fasta: Path, log) -> dict:
    """Return {versionless_protein_id: best local-alignment score to any anchor}."""
    try:
        from Bio import SeqIO
        from Bio.Align import PairwiseAligner, substitution_matrices
    except Exception as e:
        log(f"  [ALIGN] Biopython not available ({e}) - skipping homology baseline")
        return None
    if not proteome_fasta.exists() or not anchor_fasta.exists():
        log(f"  [ALIGN] missing FASTA - skipping")
        return None

    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.mode = "local"
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1

    # load anchors, build per-anchor k-mer sets
    anchors = []
    for rec in SeqIO.parse(str(anchor_fasta), "fasta"):
        s = _clean(rec.seq)
        anchors.append((s, _kmers(s)))
    log(f"  [ALIGN] {len(anchors)} anchors loaded; aligning seeded proteins "
        f"(BLOSUM62 local, k={_KMER})")

    best = {}
    n_seen = n_aligned = 0
    for rec in SeqIO.parse(str(proteome_fasta), "fasta"):
        n_seen += 1
        pid = strip_ver(rec.id.split()[0])
        s = _clean(rec.seq)
        ks = _kmers(s)
        if not ks:
            best[pid] = 0.0
            continue
        # pick the single anchor sharing the most k-mers (seed prefilter)
        shared, idx = 0, -1
        for i, (_, aks) in enumerate(anchors):
            c = len(ks & aks)
            if c > shared:
                shared, idx = c, i
        if shared < _MIN_SHARED or idx < 0:
            best[pid] = 0.0
            continue
        best[pid] = float(aligner.score(s, anchors[idx][0]))
        n_aligned += 1
        if n_aligned % 2000 == 0:
            log(f"    aligned {n_aligned:,} seeded proteins "
                f"({n_seen:,} scanned) ...")
    log(f"  [ALIGN] scanned {n_seen:,} proteins; aligned {n_aligned:,} "
        f"(rest scored 0 = no seed match)")
    return best


# ── main ─────────────────────────────────────────────────────────────────────
def main():
    cache = OUTDIR / "_cache"
    cache.mkdir(exist_ok=True)
    lines = []
    def log(s=""):
        print(s); lines.append(s)

    metric_rows = []
    recovered_rows = []

    log("=" * 72)
    log("Baseline comparison vs curated GO defense set")
    log("ESM-2 (moderate)  |  BLAST-to-anchor (matched-N)  |  InterPro-domain  |  keyword")
    log("=" * 72)

    for folder, (display, taxid) in SPECIES.items():
        val = RESULTS / folder / "04_validate_annotations" / "validated_results.csv"
        if not val.exists():
            log(f"\n[SKIP] {val} not found")
            continue

        log(f"\n### {display}  (taxid {taxid})")
        df = pd.read_csv(val)
        if "protein_id" not in df.columns:
            df = df.rename(columns={df.columns[0]: "protein_id"})
        df["pid"] = df["protein_id"].map(strip_ver)

        curated = fetch_curated_refseq(taxid, cache, log)
        y_true = df["pid"].isin(curated).astype(int).values
        n_pos = int(y_true.sum())
        n_total = len(df)
        prevalence = n_pos / n_total if n_total else 0.0
        log(f"  curated positives mapped: {n_pos:,} / {n_total:,} ({prevalence:.2%})")

        if n_pos < MIN_CURATED:
            log(f"  [SKIP baselines] fewer than {MIN_CURATED} curated positives "
                f"(annotation gap); baselines not meaningful here.")
            continue

        # --- 1. ESM-2 (moderate tier) ---
        esm_mask = df["defense_moderate"].values.astype(bool)
        n_esm = int(esm_mask.sum())
        m = metrics(esm_mask, y_true, prevalence)
        m.update(species=display, method="ESM-2 (moderate)",
                 AUPRC=round(average_precision_score(y_true, df["z_max"].values), 3),
                 ROC_AUC=round(roc_auc_score(y_true, df["z_max"].values), 3))
        metric_rows.append(m)
        log(f"  ESM-2 (moderate)  : N={m['n_called']:>6d}  P={m['precision']:.3f}  "
            f"R={m['recall']:.3f}  F1={m['F1']:.3f}  fold={m['fold_vs_background']:.2f}  "
            f"AUPRC={m['AUPRC']:.3f}")

        # --- 2. Homology-to-anchor (matched-N) ---
        # Prefer NCBI blastp if available; otherwise pure-Python Smith-Waterman.
        blast_mask = None
        def _find(patterns):
            for pat in patterns:
                hits = sorted(RESULTS.glob(pat))
                if hits:
                    return hits[0]
            return Path("___missing___")
        prot_fa = _find([f"{folder}/00_download_proteome/proteome.fasta",
                         f"{folder}/**/proteome.fasta",
                         f"{folder}/**/*.fasta"])
        anchor_fa = _find([f"{folder}/01_fetch_anchors/anchors.fasta",
                           "01_fetch_anchors/anchors.fasta",
                           "**/anchors.fasta"])
        log(f"  proteome FASTA: {prot_fa}")
        log(f"  anchor   FASTA: {anchor_fa}")
        if shutil.which("blastp") and shutil.which("makeblastdb"):
            best, hlabel = blast_best_bitscore(prot_fa, anchor_fa, log), "BLAST-to-anchor"
        else:
            best, hlabel = alignment_best_score(prot_fa, anchor_fa, log), "alignment-to-anchor (SW)"
        if best is not None:
            bit = df["pid"].map(lambda p: best.get(p, 0.0)).values
            # candidate set = top-N by score, N matched to ESM-2 moderate count
            order = np.argsort(-bit, kind="stable")
            blast_mask = np.zeros(n_total, bool)
            blast_mask[order[:n_esm]] = True
            mb = metrics(blast_mask, y_true, prevalence)
            mb.update(species=display, method=f"{hlabel} (top-{n_esm})",
                      AUPRC=round(average_precision_score(y_true, bit), 3),
                      ROC_AUC=round(roc_auc_score(y_true, bit), 3))
            metric_rows.append(mb)
            log(f"  {hlabel:17s}: N={mb['n_called']:>6d}  P={mb['precision']:.3f}  "
                f"R={mb['recall']:.3f}  F1={mb['F1']:.3f}  fold={mb['fold_vs_background']:.2f}  "
                f"AUPRC={mb['AUPRC']:.3f}")

        # --- 3. InterPro-domain ---
        domain_ids = fetch_domain_refseq(taxid, cache, log)
        domain_mask = df["pid"].isin(domain_ids).values
        md = metrics(domain_mask, y_true, prevalence)
        md.update(species=display, method="InterPro-domain", AUPRC=np.nan, ROC_AUC=np.nan)
        metric_rows.append(md)
        log(f"  InterPro-domain   : N={md['n_called']:>6d}  P={md['precision']:.3f}  "
            f"R={md['recall']:.3f}  F1={md['F1']:.3f}  fold={md['fold_vs_background']:.2f}")

        # --- 4. RefSeq-keyword ---
        if "has_defense_keyword" in df.columns:
            kw_mask = df["has_defense_keyword"].values.astype(bool)
            mk = metrics(kw_mask, y_true, prevalence)
            mk.update(species=display, method="RefSeq-keyword", AUPRC=np.nan, ROC_AUC=np.nan)
            metric_rows.append(mk)
            log(f"  RefSeq-keyword    : N={mk['n_called']:>6d}  P={mk['precision']:.3f}  "
                f"R={mk['recall']:.3f}  F1={mk['F1']:.3f}  fold={mk['fold_vs_background']:.2f}")

        # --- recovered-by breakdown among curated positives ---
        pos = (y_true == 1)
        e = esm_mask & pos
        d = domain_mask & pos
        b = (blast_mask & pos) if blast_mask is not None else np.zeros(n_total, bool)
        std = d | b                      # any "standard" method
        rec = {
            "species": display, "curated_positives": n_pos,
            "recovered_ESM2": int(e.sum()),
            "recovered_BLAST": int(b.sum()) if blast_mask is not None else None,
            "recovered_domain": int(d.sum()),
            "ESM2_only_vs_standard": int((e & ~std).sum()),
            "standard_only_vs_ESM2": int((std & ~e).sum()),
            "both_ESM2_and_standard": int((e & std).sum()),
            "recovered_by_neither": int((pos & ~e & ~std).sum()),
        }
        recovered_rows.append(rec)
        log(f"  recovered curated positives: ESM-2={rec['recovered_ESM2']}, "
            f"standard(domain|BLAST)={int(std[pos].sum())}; "
            f"ESM-2-only={rec['ESM2_only_vs_standard']}, "
            f"standard-only={rec['standard_only_vs_ESM2']}, "
            f"both={rec['both_ESM2_and_standard']}, "
            f"neither={rec['recovered_by_neither']}")

    if metric_rows:
        cols = ["species", "method", "n_called", "TP", "FP", "FN",
                "precision", "recall", "F1", "fold_vs_background", "AUPRC", "ROC_AUC"]
        pd.DataFrame(metric_rows)[cols].to_csv(OUTDIR / "baseline_metrics.csv", index=False)
    if recovered_rows:
        pd.DataFrame(recovered_rows).to_csv(OUTDIR / "recovered_breakdown.csv", index=False)
    (OUTDIR / "summary.txt").write_text("\n".join(lines), encoding="utf-8")

    log(f"\nOutputs -> {OUTDIR}")
    log("Done.")


if __name__ == "__main__":
    if not RESULTS.exists():
        sys.exit("ERROR: run from the project root (no 'results/' folder here).")
    main()
