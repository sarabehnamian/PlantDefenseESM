#!/usr/bin/env python3
"""
16_profile_baselines.py
=======================
Profile-based and low-complexity-filter-off homology baselines, and the
residual set of proteins that ESM-2 recovers but no sequence search does.
(Editor comments E4, E5 and E6.)

12_benchmark_baselines.py compared ESM-2 against a SINGLE-SEQUENCE alignment
baseline (best local-alignment score to any of the 33 anchors). The Editor
notes that single-sequence search is the weakest form of homology detection
and asks for two stronger comparisons before any claim is made about what the
embeddings add:

  E4  profile-based search, which is more sensitive than single-sequence search
  E5  sequence search with low-complexity filtering switched off

Only once those are in place is the "found by PLM, not by homology" set
meaningful -- and that residual set is what E6 asks us to characterise.

Arms
----
All arms score every proteome protein by its best bitscore to any of the 33
anchors, and all use the SAME candidate rule as 12_benchmark_baselines.py:
the top N by score, with N matched to the ESM-2 moderate-tier count. Ground
truth is the same curated GO defense set (GO:0006952 + GO:0002376).

  blastp_fwd      blastp, anchors vs proteome DB, default settings.
                  Bridge to the step-12 arm, which searched the other
                  direction; lets us see whether direction matters at all.
  blastp_nofilter blastp with -seg no and -comp_based_stats 0.          [E5]
                  Composition-based statistics is the filter that actually
                  bites for blastp (SEG is already off by default), so both
                  are switched off here.
  psiblast        psiblast, 3 iterations, anchors vs proteome DB.       [E4]
  jackhmmer       jackhmmer, 3 iterations, anchors vs proteome.         [E4]
                  Iterative profile HMM search; the most sensitive arm.
  phmmer_nofilter phmmer with --max --nonull2, i.e. all heuristic
                  filters and the biased-composition correction off.    [E5]

Every arm is optional. Missing binaries are reported and skipped rather than
aborting the run, so this is usable with whatever is installed:

  BLAST+  (blastp, psiblast, makeblastdb)  -> Windows binaries available
  HMMER   (jackhmmer, phmmer)              -> Linux/macOS/WSL, or `pip install
                                              pyhmmer` for the same algorithms

Outputs
-------
Two things the manuscript needs:

  1. a metrics table putting ESM-2 alongside every homology arm, which is what
     E1/E2/E3 have to be written from; and
  2. the residual proteins ESM-2 recovers that NO sequence method finds, with
     length, truncation status, compositional bias, repeat content, domain and
     keyword status attached -- the raw material for the E6 characterisation
     and the E7 structural-homology discussion.

Run from the project root:
    python 16_profile_baselines.py
    python 16_profile_baselines.py --species arabidopsis_thaliana --iterations 3
    python 16_profile_baselines.py --arms jackhmmer,psiblast --refresh

Outputs -> results/16_profile_baselines/
    profile_metrics_<species>.csv       - every arm: N, TP, P, R, F1, AUPRC, fold
    profile_overlap_<species>.csv       - ESM-2 vs union-of-homology vs domain
    residual_curated_<species>.csv      - curated positives found ONLY by ESM-2
    residual_candidates_<species>.csv   - all moderate candidates found by no
                                          sequence method (the larger set)
    arm_scores_<species>.csv.gz         - per-protein score for every arm
    summary.txt                         - human-readable log (also printed)
"""

import argparse
import gzip
import importlib.util
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

RESULTS = Path("results")
OUTDIR = RESULTS / "16_profile_baselines"

ESM_TRUNCATION_LIMIT = 1022      # the ESM-2 input limit used in step 02
DEFAULT_ITERATIONS = 3           # psiblast / jackhmmer iterations
DEFAULT_EVALUE = 10.0            # permissive: we rank by score, not by E cutoff
DEFAULT_INCLUSION_E = 0.001      # profile inclusion threshold per iteration
MIN_ARM_HITS = 50                # an arm producing fewer real hits than this is
                                 # treated as FAILED, not as a weak result

ALL_ARMS = ["blastp_fwd", "blastp_nofilter", "psiblast",
            "jackhmmer", "phmmer_nofilter"]


# ── reuse step 12 verbatim ───────────────────────────────────────────────────

def load_step12():
    """Import 12_benchmark_baselines.py (module name starts with a digit)."""
    path = Path(__file__).resolve().parent / "12_benchmark_baselines.py"
    if not path.exists():
        raise FileNotFoundError(
            f"Cannot find {path} (run this from the project root)")
    spec = importlib.util.spec_from_file_location("_step12_baselines", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def find_fasta(folder: str, patterns) -> Path:
    """Same lookup logic as 12_benchmark_baselines.py / 15_baseline_overlap.py."""
    for pat in patterns:
        hits = sorted(RESULTS.glob(pat))
        if hits:
            return hits[0]
    return Path("___missing___")


# ── homology arms ────────────────────────────────────────────────────────────
# Every arm returns {versionless_protein_id: best bitscore to any anchor}.
# Anchors are the QUERY and the proteome is the DATABASE, so that the profile
# arms build their profile around each anchor -- which is the whole point of
# asking for profile search. Step 12 searched the other direction; blastp_fwd
# is included so that difference is measurable rather than assumed.

def _tabular_best(path: Path, strip_ver, qcol: int, scol: int) -> dict:
    """Best score per SUBJECT (proteome protein) from a tabular hit file."""
    best = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) <= max(qcol, scol):
                continue
            sid = strip_ver(parts[qcol].split()[0])
            try:
                score = float(parts[scol])
            except ValueError:
                continue
            if score > best.get(sid, -1.0):
                best[sid] = score
    return best


def _make_blast_db(proteome_fasta: Path, tmp: Path, log) -> Path:
    db = tmp / "proteome"
    log("    building proteome BLAST database ...")
    subprocess.run(["makeblastdb", "-in", str(proteome_fasta),
                    "-dbtype", "prot", "-out", str(db)],
                   check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)
    return db


def run_blastp(proteome_fasta, anchor_fasta, strip_ver, log,
               no_filter=False, threads=4, evalue=DEFAULT_EVALUE):
    """blastp, anchors vs proteome. no_filter switches off SEG and comp-stats."""
    if not (shutil.which("blastp") and shutil.which("makeblastdb")):
        log("    [SKIP] blastp/makeblastdb not on PATH")
        return None
    with tempfile.TemporaryDirectory() as tmpd:
        tmp = Path(tmpd)
        db = _make_blast_db(proteome_fasta, tmp, log)
        out_tsv = tmp / "hits.tsv"
        cmd = ["blastp", "-query", str(anchor_fasta), "-db", str(db),
               "-evalue", str(evalue), "-max_target_seqs", "100000",
               "-num_threads", str(threads),
               "-outfmt", "6 sseqid qseqid bitscore evalue",
               "-out", str(out_tsv)]
        if no_filter:
            # SEG is already off by default for blastp; composition-based
            # statistics is the correction that actually suppresses hits.
            cmd += ["-seg", "no", "-comp_based_stats", "0"]
        log(f"    running blastp{' (filters off)' if no_filter else ''} ...")
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        return _tabular_best(out_tsv, strip_ver, qcol=0, scol=2)


def run_psiblast(proteome_fasta, anchor_fasta, strip_ver, log,
                 iterations=DEFAULT_ITERATIONS, threads=4,
                 evalue=DEFAULT_EVALUE, inclusion=DEFAULT_INCLUSION_E):
    """psiblast, iterated per anchor against the proteome database."""
    if not (shutil.which("psiblast") and shutil.which("makeblastdb")):
        log("    [SKIP] psiblast/makeblastdb not on PATH")
        return None
    from Bio import SeqIO

    with tempfile.TemporaryDirectory() as tmpd:
        tmp = Path(tmpd)
        db = _make_blast_db(proteome_fasta, tmp, log)

        # psiblast takes ONE query at a time when iterating: a multi-query run
        # only iterates on the first. Split the anchors out and loop.
        best = {}
        anchors = list(SeqIO.parse(str(anchor_fasta), "fasta"))
        log(f"    running psiblast, {iterations} iterations, "
            f"{len(anchors)} anchors ...")
        for i, rec in enumerate(anchors, 1):
            q = tmp / "query.fasta"
            SeqIO.write([rec], str(q), "fasta")
            out_tsv = tmp / "hits.tsv"
            subprocess.run(
                ["psiblast", "-query", str(q), "-db", str(db),
                 "-num_iterations", str(iterations),
                 "-evalue", str(evalue),
                 "-inclusion_ethresh", str(inclusion),
                 "-max_target_seqs", "100000",
                 "-num_threads", str(threads),
                 "-outfmt", "6 sseqid qseqid bitscore evalue",
                 "-out", str(out_tsv)],
                check=True, stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL)
            for pid, score in _tabular_best(out_tsv, strip_ver,
                                            qcol=0, scol=2).items():
                if score > best.get(pid, -1.0):
                    best[pid] = score
            if i % 5 == 0:
                log(f"      {i}/{len(anchors)} anchors done "
                    f"({len(best):,} proteins hit so far)")
        return best


def run_hmmer(proteome_fasta, anchor_fasta, strip_ver, log, tool="jackhmmer",
              iterations=DEFAULT_ITERATIONS, threads=4, extra=None):
    """jackhmmer (iterative profile) or phmmer (single-sequence) via binary."""
    if not shutil.which(tool):
        log(f"    [SKIP] {tool} not on PATH")
        return None
    with tempfile.TemporaryDirectory() as tmpd:
        tmp = Path(tmpd)
        tbl = tmp / "hits.tbl"
        cmd = [tool, "--tblout", str(tbl), "-o", str(tmp / "stdout.txt"),
               "--cpu", str(threads), "-E", str(DEFAULT_EVALUE)]
        if tool == "jackhmmer":
            cmd += ["-N", str(iterations)]
        if extra:
            cmd += list(extra)
        cmd += [str(anchor_fasta), str(proteome_fasta)]
        log(f"    running {tool}"
            f"{' ' + ' '.join(extra) if extra else ''} ...")
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        # --tblout columns: target, acc, query, acc, full-E, full-score, ...
        return _tabular_best(tbl, strip_ver, qcol=0, scol=5)


def run_pyhmmer(proteome_fasta, anchor_fasta, strip_ver, log,
                iterations=DEFAULT_ITERATIONS, threads=4):
    """pyhmmer fallback for jackhmmer, for environments without the binary."""
    try:
        import pyhmmer
    except ImportError:
        log("    [SKIP] pyhmmer not installed "
            "(pip install pyhmmer) and jackhmmer not on PATH")
        return None

    # pyhmmer moved jackhmmer between the package root and pyhmmer.hmmer, and
    # switched hit names between bytes and str across versions. Handle both.
    jack = getattr(pyhmmer, "jackhmmer", None)
    if jack is None:
        jack = getattr(getattr(pyhmmer, "hmmer", None), "jackhmmer", None)
    if jack is None:
        log("    [SKIP] this pyhmmer build exposes no jackhmmer entry point")
        return None

    def _name(x):
        return x.decode() if isinstance(x, (bytes, bytearray)) else str(x)

    try:
        alphabet = pyhmmer.easel.Alphabet.amino()
        with pyhmmer.easel.SequenceFile(str(proteome_fasta), digital=True,
                                        alphabet=alphabet) as f:
            targets = f.read_block()
        with pyhmmer.easel.SequenceFile(str(anchor_fasta), digital=True,
                                        alphabet=alphabet) as f:
            queries = list(f)
    except Exception as exc:
        log(f"    [FAIL] pyhmmer could not read the FASTA input: {exc}")
        return None

    log(f"    running pyhmmer jackhmmer, {iterations} iterations, "
        f"{len(queries)} anchors ...")
    best = {}
    first_error = None
    n_failed = 0
    for i, q in enumerate(queries, 1):
        try:
            hits = None
            # jackhmmer wants an ITERABLE of queries, not a bare sequence
            for iteration in jack([q], targets, max_iterations=iterations,
                                  cpus=threads):
                hits = iteration.hits
            if hits is None:
                n_failed += 1
                continue
            for hit in hits:
                pid = strip_ver(_name(hit.name).split()[0])
                score = float(hit.score)
                if score > best.get(pid, -1.0):
                    best[pid] = score
        except Exception as exc:
            n_failed += 1
            if first_error is None:
                first_error = f"{type(exc).__name__}: {exc}"
        if i % 10 == 0:
            log(f"      {i}/{len(queries)} anchors done "
                f"({len(best):,} proteins hit so far)")

    if n_failed:
        log(f"    [WARN] {n_failed}/{len(queries)} anchors failed; "
            f"first error was {first_error}")
    if len(best) < MIN_ARM_HITS:
        log(f"    [FAIL] pyhmmer returned only {len(best)} hits - "
            f"treating this arm as unavailable rather than reporting it")
        return None
    return best


# ── sequence features for the residual characterisation (E6 / E7) ────────────

def _low_complexity_fraction(seq: str, window: int = 20, thresh: float = 3.2):
    """Fraction of the sequence in windows of low Shannon entropy (bits).

    A 20-residue window of typical globular protein sits near 3.8-4.2 bits, so
    the 3.2-bit cut flags genuinely biased windows. SEG's default is stricter
    still; this is a permissive proxy, meant for ranking, not for masking.
    """
    if len(seq) < window:
        return 0.0
    n_low = 0
    n_win = 0
    for i in range(0, len(seq) - window + 1, window // 2):
        w = seq[i:i + window]
        counts = Counter(w)
        ent = -sum((c / window) * np.log2(c / window) for c in counts.values())
        n_win += 1
        if ent < thresh:
            n_low += 1
    return round(n_low / n_win, 4) if n_win else 0.0


def _repeat_fraction(seq: str, k: int = 6):
    """Fraction of k-mers that occur more than once -- a crude repeat proxy."""
    if len(seq) < k * 2:
        return 0.0
    kmers = [seq[i:i + k] for i in range(len(seq) - k + 1)]
    counts = Counter(kmers)
    repeated = sum(c for c in counts.values() if c > 1)
    return round(repeated / len(kmers), 4)


def sequence_features(proteome_fasta: Path, wanted: set, strip_ver) -> dict:
    """{pid: feature dict} for the proteins we need to characterise."""
    from Bio import SeqIO
    feats = {}
    for rec in SeqIO.parse(str(proteome_fasta), "fasta"):
        pid = strip_ver(rec.id.split()[0])
        if pid not in wanted:
            continue
        seq = str(rec.seq).upper()
        feats[pid] = {
            "length": len(seq),
            "truncated_by_esm": len(seq) > ESM_TRUNCATION_LIMIT,
            "low_complexity_frac": _low_complexity_fraction(seq),
            "repeat_frac_k6": _repeat_fraction(seq),
        }
    return feats


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Profile-based and filter-off homology baselines (E4/E5/E6)")
    ap.add_argument("--species", default=None,
                    help="restrict to one results folder, "
                         "e.g. arabidopsis_thaliana")
    ap.add_argument("--arms", default=",".join(ALL_ARMS),
                    help=f"comma-separated subset of: {','.join(ALL_ARMS)}")
    ap.add_argument("--iterations", type=int, default=DEFAULT_ITERATIONS,
                    help="psiblast / jackhmmer iterations (default 3)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--refresh", action="store_true",
                    help="recompute all arms instead of using cached scores")
    args = ap.parse_args()

    arms = [a.strip() for a in args.arms.split(",") if a.strip()]
    unknown = [a for a in arms if a not in ALL_ARMS]
    if unknown:
        sys.exit(f"ERROR: unknown arm(s) {unknown}; choose from {ALL_ARMS}")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    cache = OUTDIR / "_cache"
    cache.mkdir(exist_ok=True)

    lines = []
    def log(s=""):
        print(s)
        lines.append(s)

    step12 = load_step12()

    log("=" * 74)
    log("Profile-based and filter-off homology baselines  (Editor E4, E5, E6)")
    log(f"Arms requested: {', '.join(arms)}")
    log("=" * 74)

    # report what is actually runnable before doing any work
    avail = {t: bool(shutil.which(t)) for t in
             ("makeblastdb", "blastp", "psiblast", "jackhmmer", "phmmer")}
    log("\nBinaries on PATH:")
    for t, ok in avail.items():
        log(f"  {t:12s} {'yes' if ok else 'NO'}")
    if not avail["jackhmmer"]:
        try:
            import pyhmmer  # noqa: F401
            log("  pyhmmer      yes (will substitute for jackhmmer)")
        except ImportError:
            log("  pyhmmer      NO")

    species = step12.SPECIES
    if args.species:
        species = {k: v for k, v in species.items() if k == args.species}
        if not species:
            sys.exit(f"ERROR: unknown species folder '{args.species}'")

    for folder, (display, taxid) in species.items():
        val = RESULTS / folder / "04_validate_annotations" / "validated_results.csv"
        if not val.exists():
            log(f"\n[SKIP] {val} not found")
            continue

        log(f"\n{'=' * 74}")
        log(f"### {display}  (taxid {taxid})")
        log("=" * 74)

        df = pd.read_csv(val)
        if "protein_id" not in df.columns:
            df = df.rename(columns={df.columns[0]: "protein_id"})
        df["pid"] = df["protein_id"].map(step12.strip_ver)

        curated = step12.fetch_curated_refseq(taxid, cache, log)
        y_true = df["pid"].isin(curated).astype(int).values
        n_pos = int(y_true.sum())
        n_total = len(df)
        prevalence = n_pos / n_total if n_total else 0.0
        log(f"  curated positives mapped: {n_pos:,} / {n_total:,} "
            f"({prevalence:.2%})")

        if n_pos < step12.MIN_CURATED:
            log(f"  [SKIP] fewer than {step12.MIN_CURATED} curated positives; "
                f"baselines are not interpretable here.")
            continue

        prot_fa = find_fasta(folder,
                             [f"{folder}/00_download_proteome/proteome.fasta",
                              f"{folder}/**/proteome.fasta"])
        anchor_fa = find_fasta(folder,
                               [f"{folder}/01_fetch_anchors/anchors.fasta",
                                "01_fetch_anchors/anchors.fasta",
                                "**/anchors.fasta"])
        log(f"  proteome FASTA: {prot_fa}")
        log(f"  anchor   FASTA: {anchor_fa}")
        if not prot_fa.exists() or not anchor_fa.exists():
            log("  [SKIP] missing FASTA input")
            continue

        # --- ESM-2 reference arm, unchanged from step 12 ---
        esm_mask = df["defense_moderate"].values.astype(bool)
        n_esm = int(esm_mask.sum())
        metric_rows = []
        m = step12.metrics(esm_mask, y_true, prevalence)
        m.update(species=display, method="ESM-2 (moderate)",
                 AUPRC=round(average_precision_score(y_true,
                                                     df["z_max"].values), 3),
                 ROC_AUC=round(roc_auc_score(y_true, df["z_max"].values), 3))
        metric_rows.append(m)
        log(f"\n  {'ESM-2 (moderate)':22s} N={m['n_called']:>6d}  "
            f"P={m['precision']:.3f}  R={m['recall']:.3f}  "
            f"F1={m['F1']:.3f}  AUPRC={m['AUPRC']:.3f}")

        # --- homology arms ---
        arm_masks = {}
        arm_scores = {}

        for arm in arms:
            cache_file = cache / f"scores_{folder}_{arm}.csv"
            if cache_file.exists() and not args.refresh:
                log(f"\n  [{arm}] using cached scores")
                cached = pd.read_csv(cache_file)
                best = dict(zip(cached["pid"], cached["score"]))
            else:
                log(f"\n  [{arm}]")
                if arm == "blastp_fwd":
                    best = run_blastp(prot_fa, anchor_fa, step12.strip_ver, log,
                                      no_filter=False, threads=args.threads)
                elif arm == "blastp_nofilter":
                    best = run_blastp(prot_fa, anchor_fa, step12.strip_ver, log,
                                      no_filter=True, threads=args.threads)
                elif arm == "psiblast":
                    best = run_psiblast(prot_fa, anchor_fa, step12.strip_ver,
                                        log, iterations=args.iterations,
                                        threads=args.threads)
                elif arm == "jackhmmer":
                    best = run_hmmer(prot_fa, anchor_fa, step12.strip_ver, log,
                                     tool="jackhmmer",
                                     iterations=args.iterations,
                                     threads=args.threads)
                    if best is None:
                        best = run_pyhmmer(prot_fa, anchor_fa,
                                           step12.strip_ver, log,
                                           iterations=args.iterations,
                                           threads=args.threads)
                elif arm == "phmmer_nofilter":
                    best = run_hmmer(prot_fa, anchor_fa, step12.strip_ver, log,
                                     tool="phmmer", threads=args.threads,
                                     extra=["--max", "--nonull2"])
                else:
                    best = None

                if best is None:
                    log(f"    -> arm unavailable, skipped")
                    continue
                pd.DataFrame({"pid": list(best.keys()),
                              "score": list(best.values())}
                             ).to_csv(cache_file, index=False)
                log(f"    cached scores -> {cache_file.name}")

            score = df["pid"].map(lambda p: best.get(p, 0.0)).values
            n_hit = int((score > 0).sum())

            # A search that found (almost) nothing is a FAILED run, not a weak
            # result. Without this guard a top-N cut on an all-zero score
            # vector silently returns the first N proteins in file order and
            # every downstream overlap number becomes meaningless.
            if n_hit < MIN_ARM_HITS:
                log(f"    [FAIL] only {n_hit} proteins scored above zero "
                    f"(minimum {MIN_ARM_HITS}); this arm did not run correctly "
                    f"and is EXCLUDED. Delete its cache file before retrying.")
                continue

            arm_scores[arm] = score

            # identical candidate rule to step 12: top-N, N = ESM-2 moderate,
            # but never call a protein that no search actually hit
            order = np.argsort(-score, kind="stable")
            order = [i for i in order if score[i] > 0][:n_esm]
            mask = np.zeros(n_total, bool)
            mask[order] = True
            arm_masks[arm] = mask
            if len(order) < n_esm:
                log(f"    [NOTE] only {len(order):,} proteins had a nonzero "
                    f"score, fewer than the ESM-2 count of {n_esm:,}; the "
                    f"matched-N cut is limited by the number of real hits")

            ma = step12.metrics(mask, y_true, prevalence)
            auprc = (round(average_precision_score(y_true, score), 3)
                     if score.max() > 0 else float("nan"))
            rocauc = (round(roc_auc_score(y_true, score), 3)
                      if score.max() > 0 else float("nan"))
            ma.update(species=display, method=f"{arm} (top-{n_esm})",
                      AUPRC=auprc, ROC_AUC=rocauc)
            metric_rows.append(ma)
            log(f"    {arm:22s} N={ma['n_called']:>6d}  "
                f"P={ma['precision']:.3f}  R={ma['recall']:.3f}  "
                f"F1={ma['F1']:.3f}  AUPRC={auprc}  "
                f"hits={int((score > 0).sum()):,}")

        if not arm_masks:
            log("\n  [ABORT] no homology arm could be run for this species.")
            continue

        # --- domain retrieval, as in step 12 / step 15 ---
        domain_ids = step12.fetch_domain_refseq(taxid, cache, log)
        domain_mask = df["pid"].isin(domain_ids).values
        md = step12.metrics(domain_mask, y_true, prevalence)
        md.update(species=display, method="InterPro-domain",
                  AUPRC=np.nan, ROC_AUC=np.nan)
        metric_rows.append(md)
        log(f"\n  {'InterPro-domain':22s} N={md['n_called']:>6d}  "
            f"P={md['precision']:.3f}  R={md['recall']:.3f}  F1={md['F1']:.3f}")

        pd.DataFrame(metric_rows).to_csv(
            OUTDIR / f"profile_metrics_{folder}.csv", index=False)

        # --- union of every sequence-search arm ---
        seq_union = np.zeros(n_total, bool)
        for mask in arm_masks.values():
            seq_union |= mask
        log(f"\n  union of {len(arm_masks)} homology arm(s) "
            f"[{', '.join(arm_masks)}]: {int(seq_union.sum()):,} proteins")
        if len(arm_masks) < len([a for a in arms]):
            missing = [a for a in arms if a not in arm_masks]
            log(f"  [CAUTION] these arms did NOT contribute: {', '.join(missing)}. "
                f"The ESM-2-only count below is an UPPER BOUND and must not be "
                f"quoted until every requested arm has run.")

        # --- overlap: ESM-2 vs homology union vs domain ---
        pos = y_true == 1
        e, u, d = esm_mask, seq_union, domain_mask
        cells = [
            ("Recovered by all three", e & u & d),
            ("ESM-2 and homology only", e & u & ~d),
            ("ESM-2 and domain only", e & ~u & d),
            ("Homology and domain only", ~e & u & d),
            ("ESM-2 only", e & ~u & ~d),
            ("Homology only", ~e & u & ~d),
            ("Domain only", ~e & ~u & d),
            ("Missed by all three", ~e & ~u & ~d),
        ]
        overlap = pd.DataFrame([
            {"cell": name, "n": int((mask & pos).sum()),
             "pct_of_curated": round(100 * int((mask & pos).sum()) / n_pos, 1)}
            for name, mask in cells
        ])
        overlap.to_csv(OUTDIR / f"profile_overlap_{folder}.csv", index=False)
        log("")
        log(overlap.to_string(index=False))

        esm_only_n = int((e & ~u & ~d & pos).sum())
        log(f"\n  ESM-2-only among curated positives, against the FULL "
            f"homology panel: {esm_only_n}")
        log(f"  (12_benchmark_baselines.py reported 70 against single-sequence "
            f"alignment + domain alone)")
        if esm_only_n < 70:
            log(f"  -> the stronger searches absorb {70 - esm_only_n} of the 70; "
                f"Table 6 and the Results/Conclusions text must be updated")

        # --- residual sets for E6 / E7 ---
        resid_cur = pos & e & ~u & ~d
        resid_cand = e & ~u & ~d

        wanted = set(df.loc[resid_cand, "pid"])
        log(f"\n  computing sequence features for {len(wanted):,} residual "
            f"candidates ...")
        feats = sequence_features(prot_fa, wanted, step12.strip_ver)

        def build(mask, path):
            sub = df.loc[mask, ["protein_id", "pid"]].copy()
            for col in ("z_max", "best_category", "has_defense_keyword",
                        "description"):
                if col in df.columns:
                    sub[col] = df.loc[mask, col].values
            for arm, score in arm_scores.items():
                sub[f"score_{arm}"] = score[mask]
            for key in ("length", "truncated_by_esm",
                        "low_complexity_frac", "repeat_frac_k6"):
                sub[key] = sub["pid"].map(
                    lambda p, k=key: feats.get(p, {}).get(k))
            sub = sub.sort_values("z_max", ascending=False)
            sub.to_csv(path, index=False)
            return sub

        rc = build(resid_cur, OUTDIR / f"residual_curated_{folder}.csv")
        ra = build(resid_cand, OUTDIR / f"residual_candidates_{folder}.csv")
        log(f"  residual curated positives : {len(rc):,} "
            f"-> residual_curated_{folder}.csv")
        log(f"  residual moderate candidates: {len(ra):,} "
            f"-> residual_candidates_{folder}.csv")

        if len(rc):
            log("\n  Residual curated positives by ESM-2 category:")
            for cat, n in rc["best_category"].value_counts().items():
                log(f"    {cat:22s} {n:>4d}")
            log(f"  median length {rc['length'].median():.0f} aa "
                f"(all candidates: {ra['length'].median():.0f} aa)")
            log(f"  truncated by ESM-2: {int(rc['truncated_by_esm'].sum())} "
                f"/ {len(rc)}")
            log(f"  mean low-complexity fraction {rc['low_complexity_frac'].mean():.3f} "
                f"(all candidates: {ra['low_complexity_frac'].mean():.3f})")
            log(f"  mean repeat fraction {rc['repeat_frac_k6'].mean():.3f} "
                f"(all candidates: {ra['repeat_frac_k6'].mean():.3f})")

        # --- per-protein scores, for anything downstream ---
        scores_out = df[["protein_id", "pid", "z_max"]].copy()
        scores_out["curated_positive"] = pos
        scores_out["esm_moderate"] = esm_mask
        scores_out["domain"] = domain_mask
        for arm, score in arm_scores.items():
            scores_out[f"score_{arm}"] = score
            scores_out[f"called_{arm}"] = arm_masks[arm]
        gz_path = OUTDIR / f"arm_scores_{folder}.csv.gz"
        with gzip.open(gz_path, "wt") as fh:
            scores_out.to_csv(fh, index=False)
        log(f"  per-protein arm scores -> {gz_path.name}")

    (OUTDIR / "summary.txt").write_text("\n".join(lines), encoding="utf-8")
    log(f"\nOutputs -> {OUTDIR}")
    log("Done.")


if __name__ == "__main__":
    if not RESULTS.exists():
        sys.exit("ERROR: run from the project root (no 'results/' folder here).")
    main()
