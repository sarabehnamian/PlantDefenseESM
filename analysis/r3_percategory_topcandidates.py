#!/usr/bin/env python3
"""
r3_percategory_topcandidates.py  (publication rebuild)

Manuscript Figure 4: the top distinct candidates WITHIN each defense category,
so every category is represented and the weaker ones can be assessed. Figure 3
ranks by maximum Z-score and is therefore dominated by NBS-LRR.

Reviewer 1, comment 7: the previous version of this figure was drawn at
14 x 32 inches (356 x 823 mm) and then scaled down to fit the page, which is
why the labels were unreadable and the figure "vertically compressed". It is
now drawn AT FINAL SIZE to BMC's 170 mm full-page width, with the species
stacked one per panel.

Because 6 candidates x 6 categories x 3 species = 108 rows cannot be shown
legibly inside BMC's 225 mm height limit, the main-text figure shows
N_PER_CAT_MAIN per category and the complete N_PER_CAT_FULL version is written
separately for the supplementary, exactly as was done for Figure 3.

Run from the project root:
    python r3_percategory_topcandidates.py

Outputs -> results/r3_percategory/
    fig_percategory_topcandidates.png          - manuscript Figure 4
    figS_percategory_full_<species>.png        - full version, for supplementary
    r3_percategory_top_candidates.csv          - the backing table
"""

import re
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Configuration ────────────────────────────────────────────────────────────

SPECIES = [
    {"label": "A. thaliana", "run_dir": Path("results/arabidopsis_thaliana"),
     "panel": "A", "short": "arabidopsis"},
    {"label": "V. vinifera", "run_dir": Path("results/vitis_vinifera"),
     "panel": "B", "short": "vitis"},
    {"label": "O. sativa",   "run_dir": Path("results/oryza_sativa"),
     "panel": "C", "short": "rice"},
]

OUT_DIR = Path("results/r3_percategory")

N_PER_CAT_MAIN = 3         # candidates per category in the manuscript figure
N_PER_CAT_FULL = 6         # candidates per category in the supplementary version
MAKE_SUPPLEMENTARY = True  # False = manuscript figure only
SAVE_PDF = False           # True also writes the vector PDF BMC prefers
# Row labels and the heatmap share the 170 mm width, so the two trade off
# directly: at 7 pt a character costs about 1.55 mm, so 60 characters leaves
# roughly 10 mm per heatmap column, 72 characters about 8 mm, and 85 characters
# only 5 mm. Dropping ROW_LABEL_PT to 6 buys roughly 12 more characters.
DESC_CHARS = 90            # row-label truncation
ROW_LABEL_PT = 7           # drop to 6 to fit the longest RefSeq names
ABBREVIATE = True          # apply the standard shortenings below to row labels

CATEGORY_ORDER = [
    "NBS_LRR", "PR_proteins", "RLK_defense",
    "defense_signaling", "antimicrobial", "cell_death_HR",
]

CAT_SHORT = {
    "NBS_LRR": "NBS-LRR",
    "PR_proteins": "PR proteins",
    "RLK_defense": "RLK defense",
    "defense_signaling": "Defense signalling",
    "antimicrobial": "Antimicrobial",
    "cell_death_HR": "Cell death / HR",
}

# ── Publication style — BMC Bioinformatics ───────────────────────────────────
# 170 mm full-page width, max 225 mm high, ~300 dpi at final size, legible at
# the 600 px web width, lines above 0.25 pt, fonts embedded, Arial in the
# graphic, no title inside the image.

MM = 1 / 25.4
FULL_W = 170 * MM
MAX_H = 225 * MM


def _pub_font():
    from matplotlib import font_manager
    have = {f.name for f in font_manager.fontManager.ttflist}
    for name in ("Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"):
        if name in have:
            return name
    return "sans-serif"


PUB_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": [_pub_font()],
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": ROW_LABEL_PT,
    "legend.fontsize": 7,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "figure.dpi": 600,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s | %(levelname)-8s | %(message)s",
                    datefmt="%H:%M:%S")
logger = logging.getLogger("r3")


def _quiet_fonttools():
    for name in ("fontTools", "fontTools.subset", "fontTools.ttLib"):
        logging.getLogger(name).setLevel(logging.WARNING)


def _italic(name):
    return r"$\it{" + name.replace(" ", r"\ ") + "}$"


# ── data ─────────────────────────────────────────────────────────────────────

def load_validated(run_dir: Path) -> pd.DataFrame:
    path = run_dir / "04_validate_annotations" / "validated_results.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing: {path}")
    return pd.read_csv(path, index_col="protein_id")


def _clean_desc(d, pid):
    d = str(d) if pd.notna(d) else ""
    d = d.strip()
    if d.startswith(pid):
        d = d[len(pid):].strip()
    d = re.sub(r"\s*\[[^\]]*\]\s*$", "", d).strip()
    return d if d else "(no description)"



# Standard, meaning-preserving shortenings for RefSeq/TAIR descriptions. Some
# of these names run past 100 characters, which cannot be shown in full at a
# legible size inside a 170 mm figure; these rules recover most of them. The
# untouched descriptions are in r3_percategory_top_candidates.csv.
ABBREVIATIONS = [
    ("LOW QUALITY PROTEIN: ", ""),
    ("serine/threonine-protein kinase", "Ser/Thr kinase"),
    ("Serine/threonine-protein kinase", "Ser/Thr kinase"),
    ("leucine-rich repeat", "LRR"),
    ("Leucine-rich repeat", "LRR"),
    ("domain-containing protein", "domain protein"),
    ("domains-containing protein", "domain protein"),
    ("coenzyme A", "CoA"),
    (" isoform X", " X"),
]


def _abbreviate(desc):
    """Shorten a RefSeq description without changing which protein it names."""
    if not ABBREVIATE:
        return desc
    # a description made of synonyms separated by ' / ' keeps the first one
    if " / " in desc:
        desc = desc.split(" / ")[0].strip()
    for old, new in ABBREVIATIONS:
        desc = desc.replace(old, new)
    return re.sub(r"\s+", " ", desc).strip()


def _stem(desc):
    s = desc.lower()
    s = re.sub(r"low quality protein:\s*", "", s)
    s = re.sub(r"\b(isoform|precursor|partial|putative|probable|class)\b", "", s)
    s = re.sub(r"[-_]?\b\d+\b", "", s)
    s = re.sub(r"[^a-z ]", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s.replace("phenyl alanine", "phenylalanine")


def per_category_rows(df, n_per_cat):
    """Top-N distinct candidates within each category (isoforms collapsed).
    Returns a frame indexed by row label, with only the z_<cat> columns."""
    pool = df[df["defense_moderate"]].copy() if "defense_moderate" in df.columns \
        else df.copy()
    z_cols = [f"z_{c}" for c in CATEGORY_ORDER]
    picked = []
    for cat in CATEGORY_ORDER:
        zc = f"z_{cat}"
        ordered = pool.sort_values(zc, ascending=False)
        seen, kept = set(), 0
        for pid, r in ordered.iterrows():
            desc = _clean_desc(r.get("description", ""), str(pid))
            st = _stem(desc)
            if st in seen:
                continue
            seen.add(st)
            star = "" if r.get("has_defense_keyword", False) else " *"
            label = _abbreviate(desc)
            if len(label) > DESC_CHARS:
                label = label[:DESC_CHARS - 1].rstrip() + "\u2026"
            row = r[z_cols].copy()
            row.name = label + star
            picked.append(row)
            kept += 1
            if kept >= n_per_cat:
                break
        if kept < n_per_cat:
            logger.warning("only %d distinct in %s", kept, cat)
    return pd.DataFrame(picked)


# ── drawing ──────────────────────────────────────────────────────────────────

def _heatmap_panel(ax, frame, vmin, vmax, n_per_cat, show_xticks=True):
    im = ax.imshow(frame.values, cmap="RdYlBu_r", aspect="auto",
                   vmin=vmin, vmax=vmax, interpolation="nearest")
    ax.set_xticks(range(frame.shape[1]))
    if show_xticks:
        ax.set_xticklabels(
            [CAT_SHORT.get(c.replace("z_", ""), c) for c in frame.columns],
            rotation=35, ha="right", rotation_mode="anchor")
    else:
        ax.set_xticklabels([])
    ax.set_yticks(range(frame.shape[0]))
    ax.set_yticklabels(list(frame.index))
    ax.tick_params(axis="y", length=0, pad=1.5)
    ax.tick_params(axis="x", pad=1.5)
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.set_xticks(np.arange(-0.5, frame.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, frame.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.4)
    ax.tick_params(which="minor", length=0)

    # separator between category blocks, so the grouping is readable
    for k in range(1, len(CATEGORY_ORDER)):
        ax.axhline(k * n_per_cat - 0.5, color="#444444", linewidth=0.7)
    return im


def _add_colorbar(fig, im, ax_or_axes, shrink):
    cbar = fig.colorbar(im, ax=ax_or_axes, location="right", shrink=shrink,
                        pad=0.015, aspect=18)
    cbar.set_label("Z-score", size=8, labelpad=3)
    cbar.ax.tick_params(labelsize=7.5, width=0.6, length=2.5)
    cbar.outline.set_linewidth(0.6)


def make_main_figure(dfs, species, n_per_cat=N_PER_CAT_MAIN):
    logger.info("Building manuscript Figure 4 (per-category candidates) ...")
    _quiet_fonttools()

    frames = [per_category_rows(df, n_per_cat) for df in dfs]
    vmax = max(float(np.abs(f.values).max()) for f in frames)
    vmin = -vmax

    row_h = 0.118
    rows = len(CATEGORY_ORDER) * n_per_cat
    fig_h = min(MAX_H * 0.88, 3 * rows * row_h + 1.05)

    with plt.rc_context({**PUB_RC, "savefig.bbox": None}):
        fig, axes = plt.subplots(3, 1, figsize=(FULL_W, fig_h),
                                 constrained_layout=True)
        for i, (ax, sp, frame) in enumerate(zip(axes, species, frames)):
            im = _heatmap_panel(ax, frame, vmin, vmax, n_per_cat,
                                show_xticks=(i == len(axes) - 1))
            ax.set_title(r"$\bf{" + sp["panel"] + "}$   " + _italic(sp["label"]),
                         loc="left", pad=3)
        _add_colorbar(fig, im, axes, 0.35)

        png = OUT_DIR / "fig_percategory_topcandidates.png"
        fig.savefig(png)
        written = [png.name]
        if SAVE_PDF:
            pdf = OUT_DIR / "fig_percategory_topcandidates.pdf"
            fig.savefig(pdf)
            written.append(pdf.name)
        plt.close(fig)

    logger.info(f"  Saved {' + '.join(written)}  "
                f"({FULL_W * 25.4:.0f} mm wide, {n_per_cat} per category)")

    if not MAKE_SUPPLEMENTARY:
        return

    # ── full version, one PNG per species, for the supplementary ─────────
    rows_full = len(CATEGORY_ORDER) * N_PER_CAT_FULL
    with plt.rc_context({**PUB_RC, "savefig.bbox": None}):
        for sp, df in zip(species, dfs):
            frame = per_category_rows(df, N_PER_CAT_FULL)
            fig, ax = plt.subplots(
                figsize=(FULL_W, min(MAX_H, rows_full * row_h + 1.4)),
                constrained_layout=True)
            im = _heatmap_panel(ax, frame, vmin, vmax, N_PER_CAT_FULL)
            ax.set_title(r"$\bf{" + sp["panel"] + "}$   " + _italic(sp["label"]),
                         loc="left", pad=3)
            _add_colorbar(fig, im, ax, 0.4)
            fig.savefig(OUT_DIR / f"figS_percategory_full_{sp['short']}.png")
            plt.close(fig)
    logger.info(f"  Saved figS_percategory_full_*.png  "
                f"({N_PER_CAT_FULL} per category — for supplementary)")


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    dfs = [load_validated(sp["run_dir"]) for sp in SPECIES]

    # backing table (always the full N, whatever the figure shows)
    rows = []
    for sp, df in zip(SPECIES, dfs):
        pool = df[df["defense_moderate"]].copy() \
            if "defense_moderate" in df.columns else df.copy()
        for cat in CATEGORY_ORDER:
            zc = f"z_{cat}"
            ordered = pool.sort_values(zc, ascending=False)
            seen, kept = set(), 0
            for pid, r in ordered.iterrows():
                desc = _clean_desc(r.get("description", ""), str(pid))
                st = _stem(desc)
                if st in seen:
                    continue
                seen.add(st)
                rec = {"species": sp["label"], "category": cat,
                       "protein_id": pid, "clean_description": desc,
                       "ranking_z": float(r[zc]),
                       "has_defense_keyword":
                           bool(r.get("has_defense_keyword", False))}
                for c in CATEGORY_ORDER:
                    rec[f"z_{c}"] = float(r[f"z_{c}"])
                rows.append(rec)
                kept += 1
                if kept >= N_PER_CAT_FULL:
                    break
    pd.DataFrame(rows).to_csv(
        OUT_DIR / "r3_percategory_top_candidates.csv", index=False)
    logger.info("Wrote %s (%d rows)",
                OUT_DIR / "r3_percategory_top_candidates.csv", len(rows))

    make_main_figure(dfs, SPECIES)


if __name__ == "__main__":
    main()
