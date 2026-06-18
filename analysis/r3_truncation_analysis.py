"""
r3_truncation_analysis.py

Reviewer 3, Comment 8 (truncation at 1,022 residues).

Reports, for each species:
  - how many proteome proteins exceed the 1,022-residue ESM-2 input limit
    (i.e. are truncated before embedding), and the length distribution
  - how many of the 33 anchors are truncated (and which ones)
  - whether truncation is enriched in particular defense categories,
    measured among the moderate-tier candidates (truncation rate per
    assigned category vs the overall candidate rate)

No GPU required: this only reads sequence lengths and the candidate table.

Run from the project root (the folder that contains `results/`):
    python r3_truncation_analysis.py
"""

from pathlib import Path
import sys
import pandas as pd

MAX_LEN = 1022
RESULTS = Path("results")

# species key  ->  supplementary sheet name
SPECIES = {
    "arabidopsis_thaliana": "Arabidopsis_thaliana",
    "vitis_vinifera":       "Vitis_vinifera",
    "oryza_sativa":         "Oryza_sativa",
}

# supplementary workbook (protein_id + assigned_category per species)
SUPP = Path("results/08_supplementary_tables/PlantDefenseESM_Supplementary_Candidates.xlsx")
if not SUPP.exists():
    # fall back to current directory
    alt = Path("PlantDefenseESM_Supplementary_Candidates.xlsx")
    if alt.exists():
        SUPP = alt


def read_fasta_lengths(path):
    """Return {accession: sequence_length}. Accession = first token of header."""
    lengths, acc, n = {}, None, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if acc is not None:
                    lengths[acc] = n
                acc = line[1:].split()[0]
                n = 0
            else:
                n += len(line.strip())
        if acc is not None:
            lengths[acc] = n
    return lengths


def find_proteome(species):
    p = RESULTS / species / "00_download_proteome" / "proteome.fasta"
    if p.exists():
        return p
    hits = list(RESULTS.glob(f"*{species}*/**/proteome.fasta")) or \
           list(RESULTS.glob("**/proteome.fasta"))
    return hits[0] if hits else None


def find_anchor_metadata():
    for sp in SPECIES:
        p = RESULTS / sp / "01_fetch_anchors" / "anchors_metadata.csv"
        if p.exists():
            return p
    hits = list(RESULTS.glob("**/anchors_metadata.csv"))
    return hits[0] if hits else None


def main():
    print("=" * 70)
    print(f"TRUNCATION ANALYSIS (ESM-2 input limit = {MAX_LEN} residues)")
    print("=" * 70)

    # ---- 1. Anchors (same 33 across species) -------------------------------
    amd = find_anchor_metadata()
    print("\n[1] ANCHORS")
    if amd is None:
        print("    anchors_metadata.csv not found under results/ -- skipping.")
    else:
        a = pd.read_csv(amd)
        lc = next((c for c in a.columns if c.lower() == "length"), None)
        if lc is None:
            sc = next((c for c in a.columns if c.lower() in ("seq", "sequence")), None)
            a["length"] = a[sc].str.len() if sc else None
            lc = "length"
        catc = next((c for c in a.columns if c.lower() == "category"), None)
        namec = next((c for c in a.columns if c.lower() in ("name", "protein", "id", "uniprot_id")), None)
        a["truncated"] = a[lc] > MAX_LEN
        print(f"    source: {amd}")
        print(f"    total anchors: {len(a)} | truncated (>{MAX_LEN} aa): {int(a['truncated'].sum())}")
        if catc:
            for cat, g in a.groupby(catc):
                t = int(g["truncated"].sum())
                if t:
                    names = ", ".join(f"{r[namec]} ({int(r[lc])} aa)" for _, r in g[g["truncated"]].iterrows())
                    print(f"      {cat}: {t}/{len(g)} truncated  -> {names}")

    # ---- 2 & 3. Per-species proteome + per-category candidate enrichment ----
    supp = None
    if SUPP.exists():
        supp = pd.ExcelFile(SUPP)
    else:
        print(f"\n[!] supplementary workbook not found ({SUPP}); category enrichment will be skipped.")

    for species, sheet in SPECIES.items():
        print("\n" + "-" * 70)
        print(species.upper())
        fa = find_proteome(species)
        if fa is None:
            print("    proteome.fasta not found -- skipping.")
            continue
        L = read_fasta_lengths(fa)
        s = pd.Series(L)
        n_tot = len(s)
        n_tr = int((s > MAX_LEN).sum())
        print(f"    proteome: {fa}")
        print(f"    proteins: {n_tot} | truncated (>{MAX_LEN} aa): {n_tr} "
              f"({100*n_tr/n_tot:.2f}%) | max length: {int(s.max())} aa | median: {int(s.median())} aa")

        if supp is None or sheet not in supp.sheet_names:
            continue
        cand = supp.parse(sheet)
        idc = next((c for c in cand.columns if c.lower() in ("protein_id", "id")), cand.columns[0])
        catc = next((c for c in cand.columns if "categor" in c.lower()), None)
        cand["length"] = cand[idc].map(L)
        miss = cand["length"].isna().sum()
        cand = cand.dropna(subset=["length"])
        cand["truncated"] = cand["length"] > MAX_LEN
        overall = 100 * cand["truncated"].mean()
        print(f"    candidates: {len(cand)} (unmatched ids: {miss}) | "
              f"overall truncation rate: {overall:.1f}%")
        if catc:
            print("    truncation rate by assigned category (rate | enrichment vs overall):")
            rows = []
            for cat, g in cand.groupby(catc):
                rate = 100 * g["truncated"].mean()
                enr = rate / overall if overall else float("nan")
                rows.append((cat, len(g), int(g["truncated"].sum()), rate, enr))
            for cat, n, t, rate, enr in sorted(rows, key=lambda r: -r[3]):
                print(f"      {cat:18} n={n:5} trunc={t:5}  {rate:5.1f}%   {enr:4.2f}x")

    print("\n" + "=" * 70)
    print("done")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("ERROR:", e)
        sys.exit(1)
