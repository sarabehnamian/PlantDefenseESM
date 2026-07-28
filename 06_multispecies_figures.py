#!/usr/bin/env python3
"""
06b_multispecies_figures.py
============================
Generate multi-panel publication figures comparing defense gene
classification results across all three species (A. thaliana,
V. vinifera, O. sativa).

Each figure has three side-by-side panels (A / B / C), one per species.

Outputs -> results/06b_multispecies_figures/
    fig1_zscore_distributions.png    - A/B/C: Z-score histograms
    fig2_category_heatmaps.png       - A/B/C: top-50 heatmaps
    fig3_annotation_support_and_breakdown.png   - A/B/C top row: annotation-support bars
                                       D/E/F bottom row: category breakdown
    fig4_tsne_embedding_space.png    - A/B/C: t-SNE projections

Usage (run from the project root):
    python 06b_multispecies_figures.py
    
    # To skip slow t-SNE:
    SKIP_TSNE=1 python 06b_multispecies_figures.py

Requirements:
    numpy, pandas, matplotlib, seaborn, scikit-learn
"""

import os
import sys
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

# ── Configuration ─────────────────────────────────────────────────────────────

SPECIES = [
    {
        "label": "A. thaliana",
        "short": "arabidopsis",
        "run_dir": Path("results/arabidopsis_thaliana"),
        "panel": "A",
    },
    {
        "label": "V. vinifera",
        "short": "vitis",
        "run_dir": Path("results/vitis_vinifera"),
        "panel": "B",
    },
    {
        "label": "O. sativa",
        "short": "rice",
        "run_dir": Path("results/oryza_sativa"),
        "panel": "C",
    },
]

OUT_DIR = Path("results/06b_multispecies_figures")
RANDOM_SEED = 42
TSNE_N_SAMPLE = 5000
TSNE_PERPLEXITY = 30

# ── Global plot style (consistent with 07_cross_species_compare) ──────────────

plt.rcParams.update({
    "font.family": "Times New Roman",
    "font.size": 14,
    "axes.labelsize": 16,
    "axes.titlesize": 18,
    "axes.titleweight": "normal",
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "figure.dpi": 350,
    "savefig.dpi": 350,
    "savefig.bbox": "tight",
    "axes.spines.top": False,
    "axes.spines.right": False,
})

PALETTE = {
    "annotation_supported_candidate":  "#2196F3",
    "candidate_without_defense_keyword": "#FF5722",
    "keyword_only":   "#FFC107",
}

ANCHOR_STYLES = {
    "NBS_LRR":           {"marker": "*", "color": "#E53935", "size": 180},
    "PR_proteins":       {"marker": "^", "color": "#1E88E5", "size": 140},
    "RLK_defense":       {"marker": "s", "color": "#43A047", "size": 120},
    "defense_signaling": {"marker": "D", "color": "#FB8C00", "size": 120},
    "antimicrobial":     {"marker": "p", "color": "#8E24AA", "size": 150},
    "cell_death_HR":     {"marker": "H", "color": "#00ACC1", "size": 150},
}

CAT_LABELS = {
    "NBS_LRR":           "NBS-LRR (R proteins)",
    "PR_proteins":       "PR proteins",
    "RLK_defense":       "RLK defense (PTI)",
    "defense_signaling": "Defense signaling",
    "antimicrobial":     "Antimicrobial",
    "cell_death_HR":     "Cell death / HR",
}

THRESHOLD_STYLES = [
    ("lenient",  "#FFC107", ":"),
    ("moderate", "#FF9800", "--"),
    ("strict",   "#F44336", "-"),
]

# ── Logging ────────────────────────────────────────────────────────────────────

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-8s | %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("06b")

# ── Helpers ────────────────────────────────────────────────────────────────────

def load_validated(run_dir: Path) -> pd.DataFrame:
    path = run_dir / "04_validate_annotations" / "validated_results.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing: {path}")
    return pd.read_csv(path, index_col="protein_id")


def panel_label(ax, letter, fontsize=15):
    """Add bold A / B / C panel label in top-left corner."""
    ax.text(-0.12, 1.05, letter, transform=ax.transAxes,
            fontsize=fontsize, fontweight="bold", va="top", ha="left")



# ══════════════════════════════════════════════════════════════════════════════
# Publication style — BMC Bioinformatics
# ══════════════════════════════════════════════════════════════════════════════
# BMC final PDF: 85 mm half-page width, 170 mm full-page width, maximum 225 mm
# high for figure + legend, ~300 dpi at final size, legible when the web version
# is rendered 600 px wide, every line wider than 0.25 pt, all fonts embedded,
# Arial or Helvetica inside the graphic, and figure titles/legends in the
# manuscript rather than in the image file.
#
# Every figure below is therefore drawn AT FINAL SIZE: the point sizes set here
# are the point sizes that appear in print. Nothing is scaled down on placement,
# which is what made the previous versions unreadable.

MM = 1 / 25.4
FULL_W = 170 * MM          # 6.69 in — full page width
HALF_W = 85 * MM           # 3.35 in — half page width
MAX_H = 225 * MM           # 8.86 in — maximum height, figure + legend


def _pub_font():
    """First available sans-serif that BMC accepts."""
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
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "lines.linewidth": 0.9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 600,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42,        # TrueType, so fonts are embedded and editable
    "ps.fonttype": 42,
}

# Canonical category order, used for every figure so the six classes always
# appear in the same sequence.
CAT_ORDER = ["NBS_LRR", "PR_proteins", "RLK_defense",
             "defense_signaling", "antimicrobial", "cell_death_HR"]

CAT_SHORT = {
    "NBS_LRR": "NBS-LRR",
    "PR_proteins": "PR proteins",
    "RLK_defense": "RLK defense",
    "defense_signaling": "Defense signalling",
    "antimicrobial": "Antimicrobial",
    "cell_death_HR": "Cell death / HR",
}

SUPPORT_LABELS = {
    "annotation_supported_candidate": "Annotation-\nsupported",
    "candidate_without_defense_keyword": "No defense\nkeyword",
    "keyword_only": "Keyword only",
}

# Lower bound of the shared Z-axis in Figure 2. A small number of V. vinifera
# proteins sit far below this; they are excluded from the drawn range and
# counted in-panel rather than compressing the informative part of the axis.
Z_AXIS_FLOOR = -5.0

# Rows shown per species in the main-text heatmap. The fully labelled 50-row
# version goes to supplementary, which is one of the options the Reviewer
# offered for Figures 3 and 4.
SAVE_PDF = False           # True also writes the vector PDF BMC prefers
MAKE_SUPPLEMENTARY_HEATMAPS = True    # full 50-row version, one PNG per species
SHOW_ROW_DESCRIPTIONS = True   # False = accession only on the heatmap rows

HEATMAP_TOP_N_MAIN = 20
HEATMAP_TOP_N_FULL = 50
DESC_CHARS = 60


def _quiet_fonttools():
    """matplotlib's PDF font subsetting logs at INFO; this script logs at INFO."""
    for name in ("fontTools", "fontTools.subset", "fontTools.ttLib"):
        logging.getLogger(name).setLevel(logging.WARNING)


def _italic(name):
    return r"$\it{" + name.replace(" ", r"\ ") + "}$"


def _panel_title(ax, letter, label=None, **kw):
    """Bold panel letter, optionally followed by an italic species name."""
    txt = r"$\bf{" + letter + "}$"
    if label:
        txt += "   " + _italic(label)
    ax.set_title(txt, loc="left", pad=3, **kw)


def _save(fig, stem):
    """Write a 600 dpi PNG, and the vector PDF too if SAVE_PDF is on."""
    written = []
    png_path = OUT_DIR / f"{stem}.png"
    fig.savefig(png_path)
    written.append(png_path.name)
    if SAVE_PDF:
        pdf_path = OUT_DIR / f"{stem}.pdf"
        fig.savefig(pdf_path)
        written.append(pdf_path.name)
    plt.close(fig)
    w_mm = fig.get_size_inches()[0] * 25.4
    logger.info(f"  Saved {' + '.join(written)}  ({w_mm:.0f} mm wide)")


def _support_series(df):
    """Annotation-support class for every protein.

    Results produced before the R1 terminology sweep carry the old column
    ('novelty', with known_defense / novel_candidate); results produced after it
    carry 'annotation_support'. Rather than depend on either, the class is
    derived from the two flags it is defined by, so the figure is correct
    whichever version of 04_validate_annotations.py wrote the table.
    """
    if "annotation_support" in df.columns:
        return df["annotation_support"]

    legacy = {"known_defense": "annotation_supported_candidate",
              "novel_candidate": "candidate_without_defense_keyword",
              "keyword_only": "keyword_only",
              "non_candidate": "non_candidate"}
    if "novelty" in df.columns:
        return df["novelty"].map(lambda v: legacy.get(v, v))

    cand = df["defense_moderate"].astype(bool)
    kw = df["has_defense_keyword"].astype(bool)
    out = pd.Series("non_candidate", index=df.index)
    out[cand & kw] = "annotation_supported_candidate"
    out[cand & ~kw] = "candidate_without_defense_keyword"
    out[~cand & kw] = "keyword_only"
    return out


def _zcols(df):
    """Per-category Z-score columns, in canonical order."""
    return [f"z_{c}" for c in CAT_ORDER if f"z_{c}" in df.columns]


def _row_labels(frame, descs=None):
    """Accession, optionally followed by a truncated RefSeq description."""
    if descs is None:
        return list(frame.index)
    out = []
    for pid in frame.index:
        d = str(descs.get(pid, "") or "")
        d = d.split(" [")[0]                      # drop the trailing organism
        if d.startswith(str(pid)):                # ... and the repeated accession
            d = d[len(str(pid)):].lstrip()
        if len(d) > DESC_CHARS:
            d = d[:DESC_CHARS - 1].rstrip() + "\u2026"
        out.append(f"{pid}  {d}" if d else str(pid))
    return out


# ── Figure 1 (manuscript Figure 2) — Z-score distributions ───────────────────

def make_fig1(dfs, species):
    logger.info("Building fig1 (manuscript Figure 2): Z-score distributions ...")
    _quiet_fonttools()

    hi = max(float(df["z_max"].max()) for df in dfs)
    lo = max(Z_AXIS_FLOOR, min(float(df["z_max"].min()) for df in dfs))
    xlim = (lo - 0.15, hi + 0.15)

    pct_map = {"strict": 99.5, "moderate": 99.0, "lenient": 97.0}

    with plt.rc_context(PUB_RC):
        fig, axes = plt.subplots(1, 3, figsize=(FULL_W, 2.3), sharex=True)

        for ax, sp, df in zip(axes, species, dfs):
            z = df["z_max"].values
            shown = z[z >= xlim[0]]
            n_below = int(z.size - shown.size)

            counts, _, _ = ax.hist(shown, bins=90, range=xlim,
                                   color="#55636B", linewidth=0, zorder=1)
            ax.set_ylim(0, counts.max() * 1.55)   # headroom for the legend

            for thr_name, color, ls in THRESHOLD_STYLES:
                col = f"defense_{thr_name}"
                if col in df.columns and df[col].sum() > 0:
                    n = int(df[col].sum())
                    thr = df.loc[df[col], "z_max"].min()
                    ax.axvline(thr, color=color, linestyle=ls, linewidth=0.9,
                               zorder=3,
                               label=f"{thr_name} (P\u2265{pct_map[thr_name]}) "
                                     f"n={n:,}")

            ax.set_xlim(xlim)
            _panel_title(ax, sp["panel"], sp["label"])
            ax.yaxis.set_major_locator(plt.MaxNLocator(4))
            ax.yaxis.set_major_formatter(
                plt.FuncFormatter(lambda v, _: f"{v:,.0f}"))
            ax.tick_params(pad=1.5)

            if sp["panel"] == "A":
                ax.set_ylabel("Number of proteins", labelpad=2)

            leg = ax.legend(loc="upper left", frameon=False, handlelength=1.4,
                            handletextpad=0.5, borderpad=0.1, labelspacing=0.35)
            for line in leg.get_lines():
                line.set_linewidth(0.9)

            if n_below:
                ax.text(0.02, 0.04, f"{n_below:,} proteins below axis",
                        transform=ax.transAxes, ha="left", va="bottom",
                        fontsize=6.5, color="#666666")

        try:
            fig.supxlabel("Maximum Z-score across defense categories",
                          fontsize=8, y=0.02)
        except AttributeError:                       # matplotlib < 3.4
            for ax in axes:
                ax.set_xlabel("Maximum Z-score across defense categories")

        fig.tight_layout(rect=(0, 0.04, 1, 1), w_pad=1.2)
        _save(fig, "fig1_zscore_distributions")


# ── Figure 2 (manuscript Figure 3) — Category heatmaps ───────────────────────

def _heatmap_panel(ax, frame, labels, vmin, vmax, show_xticks=True):
    im = ax.imshow(frame.values, cmap="RdYlBu_r", aspect="auto",
                   vmin=vmin, vmax=vmax, interpolation="nearest")
    ax.set_xticks(range(frame.shape[1]))
    if show_xticks:
        ax.set_xticklabels([CAT_SHORT.get(c.replace("z_", ""), c)
                            for c in frame.columns],
                           rotation=35, ha="right", rotation_mode="anchor")
    else:
        ax.set_xticklabels([])
    ax.set_yticks(range(frame.shape[0]))
    ax.set_yticklabels(labels)
    ax.tick_params(axis="y", length=0, pad=1.5)
    ax.tick_params(axis="x", pad=1.5)
    for spine in ax.spines.values():
        spine.set_visible(False)
    # thin white grid, drawn as minor ticks so cells stay crisp
    ax.set_xticks(np.arange(-0.5, frame.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, frame.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.4)
    ax.tick_params(which="minor", length=0)
    return im


def make_fig2(dfs, species, top_n=HEATMAP_TOP_N_MAIN):
    """Main-text heatmap: species stacked vertically so the protein
    identifiers and RefSeq descriptions stay legible at 170 mm."""
    logger.info("Building fig2 (manuscript Figure 3): category heatmaps ...")
    _quiet_fonttools()

    frames, labels = [], []
    for df in dfs:
        cols = _zcols(df)
        top = df.nlargest(top_n, "z_max")
        frames.append(top[cols])
        descs = (top["description"].to_dict()
                 if SHOW_ROW_DESCRIPTIONS and "description" in top.columns
                 else None)
        labels.append(_row_labels(top[cols], descs))

    # Symmetric about zero and shared by all three panels: keeps the warm
    # palette of the original figure while making the species comparable.
    vmax = max(float(np.abs(f.values).max()) for f in frames)
    vmin = -vmax

    row_h = 0.118                      # inches per protein row
    panel_h = top_n * row_h
    # only the bottom panel carries the column labels, so the stack stays well
    # inside the 225 mm limit and leaves room for the caption
    fig_h = min(MAX_H * 0.88, 3 * panel_h + 1.05)

    # Exact page size rather than a tight crop: the colour bar and the rotated
    # column labels can otherwise push the cropped width past BMC's 170 mm.
    with plt.rc_context({**PUB_RC, "ytick.labelsize": 7,
                         "xtick.labelsize": 7.5, "savefig.bbox": None}):
        fig, axes = plt.subplots(3, 1, figsize=(FULL_W, fig_h),
                                 constrained_layout=True)

        for i, (ax, sp, frame, lab) in enumerate(
                zip(axes, species, frames, labels)):
            im = _heatmap_panel(ax, frame, lab, vmin, vmax,
                                show_xticks=(i == len(axes) - 1))
            _panel_title(ax, sp["panel"], sp["label"])

        cbar = fig.colorbar(im, ax=axes, location="right", shrink=0.35,
                            pad=0.015, aspect=18)
        cbar.set_label("Z-score", size=8, labelpad=3)
        cbar.ax.tick_params(labelsize=7.5, width=0.6, length=2.5)
        cbar.outline.set_linewidth(0.6)

        _save(fig, "fig2_category_heatmaps")

    # ── optional: fully labelled 50-row version for supplementary ────────
    if not MAKE_SUPPLEMENTARY_HEATMAPS:
        return

    with plt.rc_context({**PUB_RC, "ytick.labelsize": 7}):
        if True:
            for sp, df in zip(species, dfs):
                cols = _zcols(df)
                top = df.nlargest(HEATMAP_TOP_N_FULL, "z_max")
                frame = top[cols]
                descs = (top["description"].to_dict()
                         if SHOW_ROW_DESCRIPTIONS and "description" in top.columns
                         else None)
                fig, ax = plt.subplots(
                    figsize=(FULL_W,
                             min(MAX_H, HEATMAP_TOP_N_FULL * row_h + 1.4)),
                    constrained_layout=True)
                im = _heatmap_panel(ax, frame, _row_labels(frame, descs),
                                    vmin, vmax)
                _panel_title(ax, sp["panel"], sp["label"])
                cbar = fig.colorbar(im, ax=ax, shrink=0.4, pad=0.015, aspect=18)
                cbar.set_label("Z-score", size=8, labelpad=3)
                cbar.ax.tick_params(labelsize=7.5, width=0.6, length=2.5)
                cbar.outline.set_linewidth(0.6)
                out = OUT_DIR / f"figS_category_heatmaps_full_{sp['short']}.png"
                fig.savefig(out)
                plt.close(fig)
    logger.info(f"  Saved figS_category_heatmaps_full_*.png  (all "
                f"{HEATMAP_TOP_N_FULL} rows per species — for supplementary)")


# ── Figure 3 (manuscript Figure 5) — support classes + category breakdown ────

def make_fig3(dfs, species):
    logger.info("Building fig3 (manuscript Figure 5): annotation support "
                "and category breakdown ...")
    _quiet_fonttools()

    letters = [["A", "B", "C"], ["D", "E", "F"]]

    support_cats = [c for c in ("annotation_supported_candidate",
                                "candidate_without_defense_keyword",
                                "keyword_only")]
    supports = [_support_series(df) for df in dfs]
    top_max = max(max(s.value_counts().get(c, 0) for c in support_cats)
                  for s in supports)
    bot_max = max(df[df["defense_moderate"]]["best_category"]
                  .value_counts().max() for df in dfs)

    with plt.rc_context(PUB_RC):
        fig, axes = plt.subplots(2, 3, figsize=(FULL_W, 4.1))

        for col_i, (sp, df, support) in enumerate(
                zip(species, dfs, supports)):

            # top row — annotation-support classes. Horizontal bars, matching
            # the bottom row: the class names are far too long to sit side by
            # side under a 55 mm panel without colliding.
            ax = axes[0, col_i]
            counts = support.value_counts()
            cats = [c for c in support_cats if c in counts.index][::-1]
            vals = [int(counts.get(c, 0)) for c in cats]
            ax.barh(range(len(cats)), vals,
                    color=[PALETTE[c] for c in cats],
                    edgecolor="white", linewidth=0.4, height=0.62)
            for i, v in enumerate(vals):
                ax.text(v + top_max * 0.02, i, f"{v:,}", va="center",
                        ha="left", fontsize=6.8)
            ax.set_yticks(range(len(cats)))
            ax.set_yticklabels([SUPPORT_LABELS[c] for c in cats]
                               if col_i == 0 else [], linespacing=1.15)
            ax.set_xlim(0, top_max * 1.18)
            ax.xaxis.set_major_locator(plt.MaxNLocator(4))
            ax.xaxis.set_major_formatter(
                plt.FuncFormatter(lambda v, _: f"{v:,.0f}"))
            ax.tick_params(axis="y", length=0, pad=1.5)
            ax.tick_params(axis="x", pad=1.5)
            ax.set_xlabel("Number of proteins", labelpad=2)
            _panel_title(ax, letters[0][col_i], sp["label"])

            # bottom row — candidates per defense category
            ax = axes[1, col_i]
            cand = df[df["defense_moderate"]]
            cc = cand["best_category"].value_counts()
            order = [c for c in CAT_ORDER if c in cc.index][::-1]
            vals = [int(cc[c]) for c in order]
            ax.barh(range(len(order)), vals, color="#4C6EB1",
                    edgecolor="white", linewidth=0.4, height=0.68)
            for i, v in enumerate(vals):
                ax.text(v + bot_max * 0.02, i, f"{v:,}", va="center",
                        ha="left", fontsize=6.8)
            ax.set_yticks(range(len(order)))
            ax.set_yticklabels([CAT_SHORT[c] for c in order]
                               if col_i == 0 else [])
            ax.set_xlim(0, bot_max * 1.16)
            ax.xaxis.set_major_locator(plt.MaxNLocator(4))
            ax.xaxis.set_major_formatter(
                plt.FuncFormatter(lambda v, _: f"{v:,.0f}"))
            ax.tick_params(axis="y", length=0, pad=1.5)
            ax.tick_params(axis="x", pad=1.5)
            ax.set_xlabel("Number of candidates", labelpad=2)
            _panel_title(ax, letters[1][col_i], sp["label"])

        fig.tight_layout(w_pad=1.4, h_pad=1.8)
        _save(fig, "fig3_annotation_support_and_breakdown")


# ── Figure 4 (manuscript Figure 6) — t-SNE ───────────────────────────────────

def make_fig4(dfs, species):
    logger.info("Building fig4 (manuscript Figure 6): t-SNE projections ...")
    _quiet_fonttools()

    from matplotlib.lines import Line2D

    with plt.rc_context(PUB_RC):
        fig, axes = plt.subplots(1, 3, figsize=(FULL_W, 2.5))
        legend_handles = None

        for ax, sp, df in zip(axes, species, dfs):
            run_dir = sp["run_dir"]
            emb_dir = run_dir / "02_embed_proteins"
            meta_path = run_dir / "01_fetch_anchors" / "anchors_metadata.csv"

            prot_data = np.load(emb_dir / "proteome_embeddings.npz",
                                allow_pickle=True)
            anch_data = np.load(emb_dir / "anchor_embeddings.npz",
                                allow_pickle=True)
            prot_emb = prot_data["embeddings"]
            prot_ids = list(prot_data["protein_ids"])
            anch_emb = anch_data["embeddings"]
            anch_ids = list(anch_data["protein_ids"])

            cand_set = set(df[df["defense_moderate"]].index)
            is_cand = np.array([pid in cand_set for pid in prot_ids])

            rng = np.random.RandomState(RANDOM_SEED)
            n_samp = min(TSNE_N_SAMPLE, len(prot_ids))
            bg_idx = rng.choice(len(prot_ids), n_samp, replace=False)
            cand_idx = np.where(is_cand)[0]
            use_idx = np.unique(np.concatenate([bg_idx, cand_idx]))

            combined = np.vstack([prot_emb[use_idx], anch_emb])
            n_prot = len(use_idx)

            logger.info(f"  {sp['label']}: t-SNE on {combined.shape[0]} points "
                        f"({n_prot} proteome + {len(anch_ids)} anchors)")

            pca = PCA(n_components=min(50, combined.shape[0] - 1),
                      random_state=RANDOM_SEED)
            coords = TSNE(n_components=2, perplexity=TSNE_PERPLEXITY,
                          random_state=RANDOM_SEED).fit_transform(
                              pca.fit_transform(combined))

            is_cand_sub = is_cand[use_idx]
            bg_mask = ~is_cand_sub

            ax.scatter(coords[:n_prot][bg_mask, 0], coords[:n_prot][bg_mask, 1],
                       c="#D9D9D9", s=1.2, alpha=0.55, linewidths=0,
                       rasterized=True, zorder=1)
            ax.scatter(coords[:n_prot][is_cand_sub, 0],
                       coords[:n_prot][is_cand_sub, 1],
                       c="#FF7043", s=2.4, alpha=0.75, linewidths=0,
                       rasterized=True, zorder=2)

            meta = pd.read_csv(meta_path)
            id_to_cat = dict(zip(meta["uniprot_id"], meta["category"]))
            cats_plotted = set()
            for i, aid in enumerate(anch_ids):
                cat = id_to_cat.get(str(aid).split("|")[0], "unknown")
                style = ANCHOR_STYLES.get(
                    cat, {"marker": "o", "color": "black", "size": 100})
                ax.scatter(coords[n_prot + i, 0], coords[n_prot + i, 1],
                           c=style["color"], s=style["size"] / 9.0,
                           marker=style["marker"], edgecolors="black",
                           linewidths=0.35, zorder=10)
                cats_plotted.add(cat)

            # t-SNE coordinates are arbitrary: keep the frame, drop the numbers
            ax.set_xticks([])
            ax.set_yticks([])
            for side in ("top", "right"):
                ax.spines[side].set_visible(True)
            for spine in ax.spines.values():
                spine.set_linewidth(0.6)
                spine.set_color("#BBBBBB")
            ax.set_xlabel("t-SNE 1", labelpad=2)
            if sp["panel"] == "A":
                ax.set_ylabel("t-SNE 2", labelpad=2)
            _panel_title(ax, sp["panel"], sp["label"])

            if legend_handles is None:
                legend_handles = [
                    Line2D([], [], marker="o", color="w",
                           markerfacecolor="#FF7043", markersize=3.2,
                           label="Defense candidate", linestyle="None"),
                    Line2D([], [], marker="o", color="w",
                           markerfacecolor="#D9D9D9", markersize=3.2,
                           label="Other protein", linestyle="None"),
                ]
                for cat in [c for c in CAT_ORDER if c in cats_plotted]:
                    style = ANCHOR_STYLES.get(cat, {"marker": "o",
                                                    "color": "black"})
                    legend_handles.append(
                        Line2D([], [], marker=style["marker"], color="w",
                               markerfacecolor=style["color"],
                               markeredgecolor="black", markeredgewidth=0.35,
                               markersize=4.0,
                               label=f"Anchor: {CAT_SHORT.get(cat, cat)}",
                               linestyle="None"))

        fig.legend(handles=legend_handles, loc="lower center", ncol=4,
                   frameon=False, fontsize=6, handletextpad=0.4,
                   columnspacing=1.2, borderaxespad=0.0,
                   bbox_to_anchor=(0.5, -0.16))

        fig.tight_layout(w_pad=1.2)
        _save(fig, "fig4_tsne_embedding_space")



# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Verify all run directories exist
    for sp in SPECIES:
        if not sp["run_dir"].exists():
            logger.error(f"Run directory not found: {sp['run_dir']}")
            sys.exit(1)

    logger.info("Loading validated results for all three species ...")
    dfs = [load_validated(sp["run_dir"]) for sp in SPECIES]
    for sp, df in zip(SPECIES, dfs):
        logger.info(f"  {sp['label']}: {len(df):,} proteins loaded")

    make_fig1(dfs, SPECIES)
    make_fig2(dfs, SPECIES)
    make_fig3(dfs, SPECIES)

    if os.environ.get("SKIP_TSNE", "0") != "1":
        make_fig4(dfs, SPECIES)
    else:
        logger.info("Skipping t-SNE (SKIP_TSNE=1)")

    logger.info(f"\nAll figures saved to: {OUT_DIR}")
    logger.info("Done OK")


if __name__ == "__main__":
    main()
