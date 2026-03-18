#!/usr/bin/env python3
"""
01_fetch_anchors.py
===================
Retrieve curated anchor defense proteins from UniProt and save as a
structured FASTA file.  These anchors define the six defense categories
in ESM-2 embedding space.

Outputs -> results/01_fetch_anchors/
    anchors.fasta         - formatted FASTA with >ID|CATEGORY|name headers
    anchors_metadata.csv  - table of all anchors with category, name, length
    summary.yaml          - counts per category
"""

import time
import urllib.request
from io import StringIO
from pathlib import Path

import pandas as pd
import yaml
from Bio import SeqIO

from shared import (
    CORE_ANCHORS,
    DEFENSE_CATEGORIES,
    load_config,
    step_dir,
    get_logger,
)


def fetch_uniprot_fasta(uniprot_id: str, cache_dir: Path, logger) -> tuple:
    """
    Fetch a single protein FASTA from UniProt.
    Returns (sequence, full_description) or (None, None) on failure.
    """
    cache_file = cache_dir / f"{uniprot_id}.fasta"

    if cache_file.exists():
        rec = next(SeqIO.parse(str(cache_file), "fasta"))
        return str(rec.seq), rec.description

    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        resp = urllib.request.urlopen(url, timeout=30)
        text = resp.read().decode("utf-8")
        with open(cache_file, "w") as fh:
            fh.write(text)
        rec = next(SeqIO.parse(StringIO(text), "fasta"))
        return str(rec.seq), rec.description
    except Exception as exc:
        logger.warning(f"  FAIL {uniprot_id}: {exc}")
        return None, None


def main():
    cfg = load_config()
    out = step_dir(cfg, "01_fetch_anchors")
    cache = out / "_uniprot_cache"
    cache.mkdir(exist_ok=True)
    logger = get_logger("01_fetch_anchors", out)

    logger.info("=" * 65)
    logger.info("STEP 01  - Fetch anchor defense proteins from UniProt")
    logger.info("=" * 65)

    rows = []
    fasta_records = []

    for category, entries in CORE_ANCHORS.items():
        desc = DEFENSE_CATEGORIES[category]["description"]
        logger.info(f"\n{category}  -  {desc}")

        for uniprot_id, name in entries:
            seq, full_desc = fetch_uniprot_fasta(uniprot_id, cache, logger)
            if seq is None:
                continue

            logger.info(f"  OK {uniprot_id}  {name}  ({len(seq)} aa)")
            rows.append({
                "uniprot_id": uniprot_id,
                "category": category,
                "name": name,
                "length": len(seq),
                "uniprot_description": full_desc,
            })
            fasta_records.append((uniprot_id, category, name, seq))
            time.sleep(0.25)  # be polite to UniProt

    # Write FASTA
    fasta_path = out / "anchors.fasta"
    with open(fasta_path, "w") as fh:
        for uid, cat, name, seq in fasta_records:
            fh.write(f">{uid}|{cat}|{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")

    # Metadata table
    df = pd.DataFrame(rows)
    df.to_csv(out / "anchors_metadata.csv", index=False)

    # Summary
    cat_counts = df["category"].value_counts().to_dict()
    summary = {
        "total_anchors": len(df),
        "categories": {k: int(v) for k, v in cat_counts.items()},
    }
    with open(out / "summary.yaml", "w") as fh:
        yaml.dump(summary, fh, default_flow_style=False)

    logger.info(f"\nTotal anchors fetched: {len(df)}")
    for cat, n in cat_counts.items():
        logger.info(f"  {cat}: {n}")
    logger.info(f"Output -> {out}")
    logger.info("Done OK")


if __name__ == "__main__":
    main()
