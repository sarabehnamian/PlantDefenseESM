#!/usr/bin/env python3
"""
00_download_proteome.py
=======================
Download (or locate) the target plant proteome and produce a clean,
filtered FASTA file plus a summary table.

Outputs -> results/00_download_proteome/
    proteome.fasta        - filtered FASTA (>=30 aa, no X-only seqs)
    proteome_stats.csv    - per-protein lengths
    summary.yaml          - counts, length distribution
"""

import gzip
import shutil
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from Bio import SeqIO

from shared import load_config, step_dir, get_logger

# ── NCBI proteome registry ──────────────────────────────────────────────────

SPECIES_DB = {
    "vitis_vinifera": {
        "name": "Vitis vinifera (grapevine)",
        "assembly": "GCF_030704535.1",
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/"
            "030/704/535/GCF_030704535.1_ASM3070453v1/"
            "GCF_030704535.1_ASM3070453v1_protein.faa.gz"
        ),
    },
    "arabidopsis_thaliana": {
        "name": "Arabidopsis thaliana",
        "assembly": "GCF_000001735.4",
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/"
            "000/001/735/GCF_000001735.4_TAIR10.1/"
            "GCF_000001735.4_TAIR10.1_protein.faa.gz"
        ),
    },
    "solanum_lycopersicum": {
        "name": "Solanum lycopersicum (tomato)",
        "assembly": "GCF_000188115.5",
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/"
            "000/188/115/GCF_000188115.5_SL3.1/"
            "GCF_000188115.5_SL3.1_protein.faa.gz"
        ),
    },
    "oryza_sativa": {
        "name": "Oryza sativa (rice)",
        "assembly": "GCF_001433935.1",
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/"
            "001/433/935/GCF_001433935.1_IRGSP-1.0/"
            "GCF_001433935.1_IRGSP-1.0_protein.faa.gz"
        ),
    },
    "triticum_aestivum": {
        "name": "Triticum aestivum (wheat)",
        "assembly": "GCF_018294505.1",
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/"
            "018/294/505/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1/"
            "GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_protein.faa.gz"
        ),
    },
}

MIN_SEQ_LEN = 30  # skip very short fragments


def download_proteome(species: str, cache_dir: Path, logger) -> Path:
    """Return path to uncompressed FASTA, downloading if needed."""
    if species not in SPECIES_DB:
        avail = ", ".join(SPECIES_DB)
        raise ValueError(f"Unknown species '{species}'. Available: {avail}")

    info = SPECIES_DB[species]
    logger.info(f"Species : {info['name']}")
    logger.info(f"Assembly: {info['assembly']}")

    cache_dir.mkdir(parents=True, exist_ok=True)
    cached = cache_dir / f"{species}_protein.faa"
    if cached.exists():
        logger.info(f"Found cached proteome -> {cached}")
        return cached

    gz = cache_dir / f"{species}_protein.faa.gz"
    logger.info(f"Downloading from NCBI ...")
    urllib.request.urlretrieve(info["url"], gz)
    logger.info("Decompressing ...")
    with gzip.open(gz, "rb") as fi, open(cached, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    gz.unlink()
    return cached


def filter_proteome(fasta_in: Path, fasta_out: Path, logger):
    """
    Read raw FASTA -> write filtered FASTA.
    Returns (sequences dict, descriptions dict).
    """
    sequences, descriptions = {}, {}
    skipped = 0

    for rec in SeqIO.parse(str(fasta_in), "fasta"):
        # Remove stop codon (*) and non-standard residues ESM-2 cannot encode
        seq = str(rec.seq)
        for ch in "*XJBZUO":
            seq = seq.replace(ch, "")
        if len(seq) < MIN_SEQ_LEN:
            skipped += 1
            continue
        sequences[rec.id] = seq
        descriptions[rec.id] = rec.description

    # Write filtered FASTA
    with open(fasta_out, "w") as fh:
        for pid, seq in sequences.items():
            # rec.description already starts with rec.id in BioPython
            fh.write(f">{descriptions[pid]}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")

    logger.info(f"Kept {len(sequences):,} proteins (skipped {skipped})")
    return sequences, descriptions


def main():
    cfg = load_config()
    out = step_dir(cfg, "00_download_proteome")
    logger = get_logger("00_download_proteome", out)

    logger.info("=" * 65)
    logger.info("STEP 00  - Download / load proteome")
    logger.info("=" * 65)

    # Locate raw FASTA
    if cfg.get("proteome_path"):
        raw = Path(cfg["proteome_path"])
        if not raw.exists():
            raise FileNotFoundError(raw)
        logger.info(f"Using user-provided FASTA -> {raw}")
    else:
        cache = Path(cfg["base_output_dir"]) / "_cache"
        raw = download_proteome(cfg["species"], cache, logger)

    # Filter
    fasta_out = out / "proteome.fasta"
    seqs, descs = filter_proteome(raw, fasta_out, logger)

    # Stats
    lengths = np.array([len(s) for s in seqs.values()])
    summary = {
        "species": cfg.get("species", "custom"),
        "n_proteins": int(len(seqs)),
        "length_min": int(lengths.min()),
        "length_max": int(lengths.max()),
        "length_median": float(np.median(lengths)),
        "length_mean": float(np.mean(lengths)),
    }
    with open(out / "summary.yaml", "w") as fh:
        yaml.dump(summary, fh, default_flow_style=False)

    # Per-protein table
    df = pd.DataFrame([
        {"protein_id": pid, "length": len(seq), "description": descs[pid]}
        for pid, seq in seqs.items()
    ])
    df.to_csv(out / "proteome_stats.csv", index=False)

    logger.info(f"Output  -> {out}")
    logger.info(f"Proteins: {summary['n_proteins']:,}")
    logger.info(
        f"Length  : median {summary['length_median']:.0f}, "
        f"mean {summary['length_mean']:.0f}, "
        f"range [{summary['length_min']}-{summary['length_max']}]"
    )
    logger.info("Done OK")


if __name__ == "__main__":
    main()
