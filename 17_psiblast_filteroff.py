#!/usr/bin/env python3
"""
13_psiblast_filteroff.py
========================
Closes editor comments E4 (PSI-BLAST) and E5 (low-complexity filtering off)
for PlantDefenseESM.

What it does
------------
Runs four (optionally six) anchor-seeded homology-search arms against the
A. thaliana proteome, scores every proteome protein by its best bitscore to
any of the 33 anchors, and evaluates each arm against the curated GO defense
set at a candidate count matched to the ESM-2 moderate tier (n = 2,807):

  1. blastp,   default settings                     (composition-based stats ON)
  2. blastp,   -seg no -comp_based_stats 0          <- E5, filters off
  3. psiblast, 3 iterations, default                <- E4
  4. psiblast, 3 iterations, -seg no -comp_based_stats 0   <- E4 + E5
  5. phmmer,   default                              (optional, needs pyhmmer)
  6. phmmer,   --max --nonull2                      (optional, E5 for HMMER)

It then recomputes the ESM-2-only residual against the UNION of every
sequence-search arm run. That converts the manuscript's current sentence
("the ESM-2-only count is consequently an upper bound") into a measured
number, which is exactly what the editor asked for.

Finally it re-characterises whatever residual survives (length, truncation
status, low-complexity content, repeat content) using the same definitions
already used in the manuscript, so the E6/E7 paragraph can be updated with
the harder numbers.

Usage
-----
    python 13_psiblast_filteroff.py --species arabidopsis_thaliana --threads 8

Point --curated / --domain-set at whatever files 09_benchmark_curated.py and
12_benchmark_baselines.py wrote; the loaders auto-detect the accession column.

Outputs (results/13_psiblast_filteroff/)
----------------------------------------
    arm_scores.csv            best bitscore per protein per arm
    arm_metrics.csv           precision / recall / F1 / AUPRC per arm  -> Table 8
    recovery_matrix.csv       per curated protein, which methods found it
    residual_esm_only.csv     ESM-2-only proteins vs the full search panel
    residual_features.csv     length / truncation / low-complexity / repeats
    summary.md               numbers pre-formatted for manuscript + letter

Runtime
-------
blastp arms: minutes. psiblast arms: ~30-60 min each on 8 threads (33 queries
x 3 iterations vs 48k sequences; psiblast cannot batch multi-FASTA queries
when -num_iterations > 1, so they are looped). phmmer --max is the slowest.

Requires BLAST+ >= 2.12 on PATH. Install with ONE of:
    conda install -c bioconda -c conda-forge blast
    sudo apt-get install ncbi-blast+
    # or the static tarball from
    # https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
"""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from sklearn.metrics import average_precision_score
except ImportError:
    average_precision_score = None


# --------------------------------------------------------------------------
# small utilities
# --------------------------------------------------------------------------

ACC_COLUMN_HINTS = [
    "protein_id", "protein_accession", "accession", "refseq", "refseq_id",
    "protein", "sseqid", "id", "acc", "target", "gene_id",
]


def die(msg: str) -> None:
    print(f"\nERROR: {msg}\n", file=sys.stderr)
    sys.exit(1)


def read_fasta(path: Path) -> dict[str, str]:
    """Minimal FASTA reader. Key = first whitespace-delimited token of header."""
    seqs: dict[str, str] = {}
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def load_accessions(path: Path) -> set[str]:
    """Load a set of accessions from CSV/TSV/TXT, auto-detecting the column."""
    if not path.exists():
        die(f"file not found: {path}")

    if path.suffix.lower() in {".txt", ".list"}:
        return {ln.strip() for ln in open(path) if ln.strip()}

    sep = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    df = pd.read_csv(path, sep=sep)

    if df.shape[1] == 1:
        return set(df.iloc[:, 0].astype(str).str.strip())

    lower = {c.lower(): c for c in df.columns}
    for hint in ACC_COLUMN_HINTS:
        if hint in lower:
            return set(df[lower[hint]].dropna().astype(str).str.strip())

    # fall back: first column that looks like RefSeq protein accessions
    for c in df.columns:
        vals = df[c].dropna().astype(str).head(50)
        if len(vals) and (vals.str.match(r"^[NXWAY]P_\d+").mean() > 0.8):
            return set(df[c].dropna().astype(str).str.strip())

    die(
        f"could not find an accession column in {path}. "
        f"Columns present: {list(df.columns)}. "
        f"Pass a single-column file, or rename the column to 'protein_id'."
    )


def check_binary(name: str) -> None:
    if shutil.which(name) is None:
        die(
            f"'{name}' not found on PATH.\n"
            f"  conda install -c bioconda -c conda-forge blast\n"
            f"  sudo apt-get install ncbi-blast+\n"
            f"  or download from "
            f"https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
        )


def run(cmd: list[str], log: Path | None = None) -> None:
    print("  $ " + " ".join(str(c) for c in cmd), flush=True)
    with open(log, "a") if log else subprocess.DEVNULL as errfh:
        proc = subprocess.run(cmd, stderr=errfh if log else subprocess.PIPE)
    if proc.returncode != 0:
        tail = ""
        if log and Path(log).exists():
            tail = "".join(open(log).readlines()[-20:])
        elif proc.stderr:
            tail = proc.stderr.decode(errors="replace")[-2000:]
        die(f"command failed (exit {proc.returncode}):\n{tail}")


# --------------------------------------------------------------------------
# sequence feature definitions (identical to those used in the manuscript)
# --------------------------------------------------------------------------

def low_complexity_fraction(seq: str, window: int = 20, bits: float = 3.2) -> float:
    """Fraction of `window`-residue windows with Shannon entropy below `bits`."""
    if len(seq) < window:
        return 0.0
    n_low = 0
    n_win = len(seq) - window + 1
    for i in range(n_win):
        counts = Counter(seq[i:i + window])
        h = -sum((c / window) * math.log2(c / window) for c in counts.values())
        if h < bits:
            n_low += 1
    return n_low / n_win


def repeat_fraction(seq: str, k: int = 6) -> float:
    """Fraction of k-mers that occur more than once."""
    if len(seq) < k:
        return 0.0
    kmers = [seq[i:i + k] for i in range(len(seq) - k + 1)]
    counts = Counter(kmers)
    repeated = sum(c for c in counts.values() if c > 1)
    return repeated / len(kmers)


# --------------------------------------------------------------------------
# search arms
# --------------------------------------------------------------------------

def parse_best_bitscores(tab_path: Path) -> dict[str, float]:
    """Max bitscore per subject across all queries and all PSI-BLAST iterations."""
    best: dict[str, float] = defaultdict(float)
    if not tab_path.exists():
        return best
    with open(tab_path) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row or row[0].startswith("#") or len(row) < 12:
                continue
            sseqid, bitscore = row[1], row[11]
            try:
                bs = float(bitscore)
            except ValueError:
                continue
            if bs > best[sseqid]:
                best[sseqid] = bs
    return best


def arm_blastp(anchors_fa, db, out_tab, threads, filters_off, log):
    cmd = [
        "blastp", "-query", str(anchors_fa), "-db", str(db),
        "-outfmt", "6", "-evalue", "10", "-max_target_seqs", "50000",
        "-num_threads", str(threads), "-out", str(out_tab),
    ]
    if filters_off:
        # NOTE: -seg is already 'no' by default for blastp; the correction that
        # actually suppresses hits is composition-based statistics. We disable
        # both so the arm tests what the editor intends.
        cmd += ["-seg", "no", "-comp_based_stats", "0"]
    run(cmd, log)
    return parse_best_bitscores(out_tab)


def arm_psiblast(anchor_seqs, db, workdir, out_tab, threads, filters_off,
                 iterations, log):
    """psiblast rejects multi-FASTA queries when -num_iterations > 1, so loop."""
    if out_tab.exists():
        out_tab.unlink()
    for i, (aid, seq) in enumerate(sorted(anchor_seqs.items()), 1):
        qfa = workdir / f"{aid}.faa"
        qfa.write_text(f">{aid}\n{seq}\n")
        tmp = workdir / f"{aid}.tsv"
        cmd = [
            "psiblast", "-query", str(qfa), "-db", str(db),
            "-num_iterations", str(iterations),
            "-inclusion_ethresh", "0.001",
            "-evalue", "10", "-outfmt", "6",
            "-max_target_seqs", "50000",
            "-num_threads", str(threads), "-out", str(tmp),
        ]
        if filters_off:
            cmd += ["-seg", "no", "-comp_based_stats", "0"]
        print(f"  [{i}/{len(anchor_seqs)}] {aid}", flush=True)
        run(cmd, log)
        if tmp.exists():
            with open(out_tab, "a") as agg, open(tmp) as src:
                agg.writelines(src.readlines())
            tmp.unlink()
        qfa.unlink()
    return parse_best_bitscores(out_tab)


def arm_phmmer(anchor_seqs, proteome_seqs, max_mode: bool):
    """Optional HMMER arm via pyhmmer. --max == F1=F2=F3=1.0, bias filter off."""
    try:
        import pyhmmer
    except ImportError:
        print("  pyhmmer not installed - skipping phmmer arm", flush=True)
        return None

    alpha = pyhmmer.easel.Alphabet.amino()
    targets = [
        pyhmmer.easel.TextSequence(name=acc.encode(), sequence=seq).digitize(alpha)
        for acc, seq in proteome_seqs.items()
    ]
    queries = [
        pyhmmer.easel.TextSequence(name=aid.encode(), sequence=seq).digitize(alpha)
        for aid, seq in anchor_seqs.items()
    ]
    kwargs = dict(bias_filter=False, null2=False, F1=1.0, F2=1.0, F3=1.0) \
        if max_mode else {}

    best: dict[str, float] = defaultdict(float)
    for hits in pyhmmer.hmmer.phmmer(queries, targets, **kwargs):
        for hit in hits:
            acc = hit.name.decode()
            if hit.score > best[acc]:
                best[acc] = hit.score
    return best


# --------------------------------------------------------------------------
# evaluation
# --------------------------------------------------------------------------

def evaluate(scores: dict[str, float], all_accs: list[str],
             curated: set[str], n_matched: int) -> tuple[dict, set[str]]:
    """Matched-count candidate set + precision/recall/F1/AUPRC."""
    vec = np.array([scores.get(a, 0.0) for a in all_accs])
    order = np.argsort(-vec, kind="stable")
    cand = {all_accs[i] for i in order[:n_matched] if vec[i] > 0}

    tp = len(cand & curated)
    prec = tp / len(cand) if cand else 0.0
    rec = tp / len(curated) if curated else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0

    auprc = float("nan")
    if average_precision_score is not None:
        y = np.array([1 if a in curated else 0 for a in all_accs])
        if y.sum():
            auprc = float(average_precision_score(y, vec))

    return (
        {"candidates": len(cand), "recovered_TP": tp,
         "precision": round(prec, 3), "recall": round(rec, 3),
         "F1": round(f1, 3), "AUPRC": round(auprc, 3)},
        cand,
    )


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--species", default="arabidopsis_thaliana")
    p.add_argument("--results-root", default="results", type=Path)
    p.add_argument("--proteome", type=Path,
                   help="default: results/<species>/00_download_proteome/proteome.fasta")
    p.add_argument("--anchors", type=Path,
                   help="default: results/<species>/01_fetch_anchors/anchors.fasta")
    p.add_argument("--curated", type=Path, required=True,
                   help="curated GO defense accessions (from 09_benchmark_curated.py)")
    p.add_argument("--esm-candidates", type=Path,
                   help="default: results/<species>/05_extract_candidates/candidates_moderate.csv")
    p.add_argument("--domain-set", type=Path,
                   help="InterPro defense-domain accessions (from 12_benchmark_baselines.py)")
    p.add_argument("--outdir", type=Path, default=Path("results/13_psiblast_filteroff"))
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--iterations", type=int, default=3)
    p.add_argument("--n-matched", type=int, default=None,
                   help="candidate count to match; default = size of ESM moderate tier")
    p.add_argument("--max-len", type=int, default=1022,
                   help="ESM-2 input limit, for truncation flags")
    p.add_argument("--skip-psiblast", action="store_true")
    p.add_argument("--with-phmmer", action="store_true")
    p.add_argument("--check", action="store_true", help="verify inputs and exit")
    args = p.parse_args()

    sroot = args.results_root / args.species
    proteome_fa = args.proteome or sroot / "00_download_proteome" / "proteome.fasta"
    anchors_fa = args.anchors or sroot / "01_fetch_anchors" / "anchors.fasta"
    esm_cand_f = args.esm_candidates or sroot / "05_extract_candidates" / "candidates_moderate.csv"

    check_binary("makeblastdb")
    check_binary("blastp")
    if not args.skip_psiblast:
        check_binary("psiblast")
    for f in (proteome_fa, anchors_fa, args.curated, esm_cand_f):
        if not Path(f).exists():
            die(f"missing input: {f}")

    print("Loading inputs ...")
    proteome = read_fasta(proteome_fa)
    anchors_raw = read_fasta(anchors_fa)
    # anchor headers are 'Q39214|NBS_LRR|RPM1'; pipes confuse BLAST ID parsing,
    # so rewrite to safe IDs and keep a map for the record.
    anchors = {f"anchor{i:02d}": s for i, (h, s) in enumerate(sorted(anchors_raw.items()))}
    anchor_map = {f"anchor{i:02d}": h for i, (h, _) in enumerate(sorted(anchors_raw.items()))}

    curated = load_accessions(args.curated) & set(proteome)
    esm_cand = load_accessions(esm_cand_f) & set(proteome)
    domain_set = (load_accessions(args.domain_set) & set(proteome)) if args.domain_set else set()
    n_matched = args.n_matched or len(esm_cand)

    print(f"  proteome        : {len(proteome):,} proteins")
    print(f"  anchors         : {len(anchors)}")
    print(f"  curated GO set  : {len(curated):,}")
    print(f"  ESM moderate    : {len(esm_cand):,}  (matched count n = {n_matched:,})")
    print(f"  domain set      : {len(domain_set):,}" if domain_set else "  domain set      : (not supplied)")
    if args.check:
        print("\n--check passed; inputs look sane.")
        return

    args.outdir.mkdir(parents=True, exist_ok=True)
    log = args.outdir / "search.log"
    work = Path(tempfile.mkdtemp(prefix="pdesm_"))
    safe_anchors = work / "anchors_safe.faa"
    safe_anchors.write_text("".join(f">{k}\n{v}\n" for k, v in anchors.items()))

    print("\nBuilding BLAST database ...")
    db = work / "proteome_db"
    run(["makeblastdb", "-in", str(proteome_fa), "-dbtype", "prot",
         "-out", str(db), "-title", "proteome"], log)

    all_accs = sorted(proteome)
    arms: dict[str, dict[str, float]] = {}

    print("\n[1/4] blastp, default")
    arms["blastp_default"] = arm_blastp(
        safe_anchors, db, args.outdir / "blastp_default.tsv",
        args.threads, False, log)

    print("\n[2/4] blastp, low-complexity + composition filters OFF")
    arms["blastp_filteroff"] = arm_blastp(
        safe_anchors, db, args.outdir / "blastp_filteroff.tsv",
        args.threads, True, log)

    if not args.skip_psiblast:
        print(f"\n[3/4] psiblast, {args.iterations} iterations, default")
        arms["psiblast_default"] = arm_psiblast(
            anchors, db, work, args.outdir / "psiblast_default.tsv",
            args.threads, False, args.iterations, log)

        print(f"\n[4/4] psiblast, {args.iterations} iterations, filters OFF")
        arms["psiblast_filteroff"] = arm_psiblast(
            anchors, db, work, args.outdir / "psiblast_filteroff.tsv",
            args.threads, True, args.iterations, log)

    if args.with_phmmer:
        print("\n[extra] phmmer, default")
        r = arm_phmmer(anchors, proteome, max_mode=False)
        if r is not None:
            arms["phmmer_default"] = r
        print("\n[extra] phmmer, --max --nonull2")
        r = arm_phmmer(anchors, proteome, max_mode=True)
        if r is not None:
            arms["phmmer_max"] = r

    # ---------------- metrics ----------------
    print("\nEvaluating arms against the curated GO defense set ...")
    rows, arm_sets = [], {}
    esm_tp = len(esm_cand & curated)
    rows.append({
        "method": "ESM-2 (moderate tier)", "candidates": len(esm_cand),
        "recovered_TP": esm_tp,
        "precision": round(esm_tp / len(esm_cand), 3) if esm_cand else 0,
        "recall": round(esm_tp / len(curated), 3) if curated else 0,
        "F1": round(2 * esm_tp / (len(esm_cand) + len(curated)), 3),
        "AUPRC": "(from 09_benchmark_curated.py)",
    })
    for name, sc in arms.items():
        m, cand = evaluate(sc, all_accs, curated, n_matched)
        m["method"] = name
        rows.append(m)
        arm_sets[name] = cand
        print(f"  {name:22s} TP={m['recovered_TP']:5d}  "
              f"P={m['precision']:.3f}  R={m['recall']:.3f}  AUPRC={m['AUPRC']}")

    metrics = pd.DataFrame(rows)[
        ["method", "candidates", "recovered_TP", "precision", "recall", "F1", "AUPRC"]]
    metrics.to_csv(args.outdir / "arm_metrics.csv", index=False)

    pd.DataFrame(
        {"protein_id": all_accs,
         **{n: [arms[n].get(a, 0.0) for a in all_accs] for n in arms}}
    ).to_csv(args.outdir / "arm_scores.csv", index=False)

    # ---------------- residual vs the full search panel ----------------
    search_union = set().union(*arm_sets.values()) if arm_sets else set()
    if domain_set:
        search_union |= domain_set
    residual = (esm_cand & curated) - search_union

    print(f"\nESM-2-only residual against the full search panel: {len(residual)}")
    print("  (manuscript currently reports 72 against jackhmmer + domain only)")

    rec_rows = []
    for acc in sorted(curated):
        row = {"protein_id": acc, "ESM2": acc in esm_cand}
        for n, s in arm_sets.items():
            row[n] = acc in s
        if domain_set:
            row["interpro_domain"] = acc in domain_set
        row["esm_only_vs_panel"] = acc in residual
        rec_rows.append(row)
    pd.DataFrame(rec_rows).to_csv(args.outdir / "recovery_matrix.csv", index=False)

    # ---------------- characterise the residual ----------------
    print("Characterising the residual (this is the slow bit) ...")

    def describe(accs: set[str]) -> pd.DataFrame:
        return pd.DataFrame([{
            "protein_id": a,
            "length": len(proteome[a]),
            "truncated_by_esm": len(proteome[a]) > args.max_len,
            "low_complexity_frac": round(low_complexity_fraction(proteome[a]), 4),
            "repeat_frac": round(repeat_fraction(proteome[a]), 4),
        } for a in sorted(accs)])

    res_df = describe(residual)
    ref_df = describe(esm_cand)
    res_df.to_csv(args.outdir / "residual_esm_only.csv", index=False)

    feat = pd.DataFrame([
        {"set": "ESM-2-only residual", "n": len(res_df),
         "median_length": res_df["length"].median() if len(res_df) else np.nan,
         "mean_low_complexity": res_df["low_complexity_frac"].mean() if len(res_df) else np.nan,
         "mean_repeat": res_df["repeat_frac"].mean() if len(res_df) else np.nan,
         "n_truncated": int(res_df["truncated_by_esm"].sum()) if len(res_df) else 0},
        {"set": "all ESM-2 moderate-tier candidates", "n": len(ref_df),
         "median_length": ref_df["length"].median(),
         "mean_low_complexity": ref_df["low_complexity_frac"].mean(),
         "mean_repeat": ref_df["repeat_frac"].mean(),
         "n_truncated": int(ref_df["truncated_by_esm"].sum())},
    ])
    feat.to_csv(args.outdir / "residual_features.csv", index=False)

    # ---------------- summary ----------------
    with open(args.outdir / "summary.md", "w") as fh:
        fh.write("# E4 / E5 search arms — results\n\n")
        fh.write(f"Species: {args.species}\n\n")
        fh.write(f"Curated GO defense set: {len(curated):,} proteins\n\n")
        fh.write(f"Matched candidate count: n = {n_matched:,}\n\n")
        fh.write("## Per-arm performance (extends Table 8)\n\n")
        fh.write(metrics.to_markdown(index=False))
        fh.write("\n\n## ESM-2-only residual\n\n")
        fh.write(f"- vs. jackhmmer + domain (as published): 72\n")
        fh.write(f"- vs. the full panel run here ({', '.join(arm_sets)}"
                 f"{' + InterPro domain' if domain_set else ''}): "
                 f"**{len(residual)}**\n\n")
        fh.write("## Residual sequence features\n\n")
        fh.write(feat.to_markdown(index=False))
        fh.write("\n\n## Anchor ID map\n\n")
        for k, v in anchor_map.items():
            fh.write(f"- {k} = {v}\n")

    shutil.rmtree(work, ignore_errors=True)
    print(f"\nDone. Outputs in {args.outdir}/")
    print("Read summary.md first — it is formatted for pasting into the letter.")


if __name__ == "__main__":
    main()
