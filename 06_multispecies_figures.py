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
    fig3_novelty_and_breakdown.png   - A/B/C top row: novelty bars
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
    "known_defense":  "#2196F3",
    "novel_candidate": "#FF5722",
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


# ── Figure 1 — Z-score distributions ──────────────────────────────────────────

def make_fig1(dfs, species):
    logger.info("Building Fig 1: Z-score distributions ...")
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=False)

    for ax, sp, df in zip(axes, species, dfs):
        z = df["z_max"].values
        ax.hist(z, bins=120, color="#455A64", alpha=0.82,
                edgecolor="white", linewidth=0.3)

        for thr_name, color, ls in THRESHOLD_STYLES:
            col = f"defense_{thr_name}"
            if col in df.columns and df[col].sum() > 0:
                n   = int(df[col].sum())
                thr = df.loc[df[col], "z_max"].min()
                pct_map = {"strict": 99.5, "moderate": 99.0, "lenient": 97.0}
                ax.axvline(thr, color=color, linestyle=ls, linewidth=2,
                           label=f"{thr_name} (P≥{pct_map[thr_name]}): {n:,}")

        ax.set_xlabel("Maximum Z-score across defense categories")
        # Only show y-label on first panel (leftmost)
        if sp["panel"] == "A":
            ax.set_ylabel("Number of proteins")
        else:
            ax.set_ylabel("")
        italic_title = "$\\it{" + sp["label"].replace(" ", "\\ ") + "}$"
        ax.set_title(italic_title)
        ax.legend(fontsize=12, loc="upper left")
        panel_label(ax, sp["panel"])

    fig.suptitle("Distribution of defense similarity scores",
                 fontsize=18, y=1.02)
    fig.tight_layout()
    path = OUT_DIR / "fig1_zscore_distributions.png"
    fig.savefig(path)
    plt.close(fig)
    logger.info(f"  Saved {path.name}")


# ── Figure 2 — Category heatmaps ──────────────────────────────────────────────

def make_fig2(dfs, species, top_n=50):
    logger.info("Building Fig 2: Category heatmaps ...")
    fig, axes = plt.subplots(1, 3, figsize=(24, max(10, top_n * 0.30)))
    
    # Nice column labels for heatmap x-axis
    col_labels = {
        "NBS_LRR": "NBS LRR",
        "PR_proteins": "PR proteins",
        "RLK_defense": "RLK defense",
        "defense_signaling": "defense signaling",
        "antimicrobial": "antimicrobial",
        "cell_death_HR": "cell death/HR",
    }

    for ax, sp, df in zip(axes, species, dfs):
        z_cols = [c for c in df.columns if c.startswith("z_") and c != "z_max"]
        top    = df.nlargest(top_n, "z_max")[z_cols].copy()
        top.columns = [col_labels.get(c.replace("z_", ""), c.replace("z_", "").replace("_", " ")) 
                       for c in top.columns]

        sns.heatmap(top, cmap="RdYlBu_r", center=0, ax=ax,
                    xticklabels=True, yticklabels=True,
                    linewidths=0.4, linecolor="white",
                    cbar_kws={"label": "Z-score", "shrink": 0.5})
        italic_title = "$\\it{" + sp["label"].replace(" ", "\\ ") + "}$"
        ax.set_title(italic_title, fontsize=24)
        ax.set_ylabel("")
        ax.tick_params(axis="y", labelsize=14)
        ax.tick_params(axis="x", labelsize=22, rotation=45)
        panel_label(ax, sp["panel"])

    fig.suptitle(f"Top {top_n} defense candidates — Z-scores by category",
                 fontsize=26, y=1.01)
    fig.tight_layout()
    path = OUT_DIR / "fig2_category_heatmaps.png"
    fig.savefig(path)
    plt.close(fig)
    logger.info(f"  Saved {path.name}")


# ── Figure 3 — Novelty bars (top) + Category breakdown (bottom) ───────────────

def make_fig3(dfs, species):
    logger.info("Building Fig 3: Novelty and category breakdown ...")
    fig = plt.figure(figsize=(24, 12))
    gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.55, wspace=0.30,
                            top=0.88, bottom=0.08)

    panel_letters = [["A", "B", "C"], ["D", "E", "F"]]
    
    # Row titles at the top of each row
    row_titles = [
        "Candidate classification (moderate threshold)",
        "Defense candidates by category (moderate)"
    ]

    for col_i, (sp, df) in enumerate(zip(species, dfs)):

        # ── top row: novelty bar chart ──
        ax_top = fig.add_subplot(gs[0, col_i])
        counts = df["novelty"].value_counts()
        plot_cats = [c for c in ["known_defense", "novel_candidate", "keyword_only"]
                     if c in counts.index]
        vals   = [counts.get(c, 0) for c in plot_cats]
        colors = [PALETTE[c] for c in plot_cats]
        labels = [c.replace("_", " ").title() for c in plot_cats]

        bars = ax_top.bar(labels, vals, color=colors,
                          edgecolor="white", linewidth=0.8)
        for bar, v in zip(bars, vals):
            ax_top.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + max(vals) * 0.01,
                        f"{v:,}", ha="center", va="bottom",
                        fontsize=14)
        # Only show y-label on leftmost panel
        if col_i == 0:
            ax_top.set_ylabel("Number of proteins", fontsize=14)
        else:
            ax_top.set_ylabel("")
        italic_title = "$\\it{" + sp["label"].replace(" ", "\\ ") + "}$"
        ax_top.set_title(italic_title, fontsize=20)
        ax_top.tick_params(axis="x", labelsize=14, rotation=25)
        ax_top.tick_params(axis="y", labelsize=14)
        plt.setp(ax_top.xaxis.get_majorticklabels(), ha="right")
        panel_label(ax_top, panel_letters[0][col_i])

        # ── bottom row: category breakdown ──
        ax_bot = fig.add_subplot(gs[1, col_i])
        cands     = df[df["defense_moderate"]]
        cat_counts = cands["best_category"].value_counts().sort_values()
        
        # Format category names (remove underscores)
        cat_labels_nice = {
            "NBS_LRR": "NBS LRR",
            "PR_proteins": "PR proteins",
            "RLK_defense": "RLK defense",
            "defense_signaling": "defense signaling",
            "antimicrobial": "antimicrobial",
            "cell_death_HR": "cell death/HR",
        }
        cat_counts.index = [cat_labels_nice.get(c, c.replace("_", " ")) for c in cat_counts.index]
        
        cat_counts.plot(kind="barh", ax=ax_bot,
                        color="#5C6BC0", edgecolor="white")
        # All panels get x-label
        ax_bot.set_xlabel("Number of candidates", fontsize=14)
        # Remove "best_category" y-axis label from all panels
        ax_bot.set_ylabel("")
        # Only show y-tick labels (category names) on leftmost panel
        if col_i != 0:
            ax_bot.set_yticklabels([])
        italic_title = "$\\it{" + sp["label"].replace(" ", "\\ ") + "}$"
        ax_bot.set_title(italic_title, fontsize=20)
        ax_bot.tick_params(axis="both", labelsize=14)
        panel_label(ax_bot, panel_letters[1][col_i])

    # Row titles at the top of each row (centered)
    fig.text(0.5, 0.95, row_titles[0], ha="center", va="bottom", fontsize=20)
    fig.text(0.5, 0.47, row_titles[1], ha="center", va="bottom", fontsize=20)

    fig.suptitle("Defense candidate classification across species",
                 fontsize=22, y=1.02)

    path = OUT_DIR / "fig3_novelty_and_breakdown.png"
    fig.savefig(path)
    plt.close(fig)
    logger.info(f"  Saved {path.name}")


# ── Figure 4 — t-SNE ──────────────────────────────────────────────────────────

def make_fig4(dfs, species):
    logger.info("Building Fig 4: t-SNE projections ...")
    fig, axes = plt.subplots(1, 3, figsize=(24, 8))

    legend_handles = None   # build once, reuse

    for ax, sp, df in zip(axes, species, dfs):
        run_dir = sp["run_dir"]
        emb_dir = run_dir / "02_embed_proteins"
        meta_path = run_dir / "01_fetch_anchors" / "anchors_metadata.csv"

        # Load embeddings
        prot_data = np.load(emb_dir / "proteome_embeddings.npz",
                            allow_pickle=True)
        anch_data = np.load(emb_dir / "anchor_embeddings.npz",
                            allow_pickle=True)
        prot_emb  = prot_data["embeddings"]
        prot_ids  = list(prot_data["protein_ids"])
        anch_emb  = anch_data["embeddings"]
        anch_ids  = list(anch_data["protein_ids"])

        # Candidate mask
        cand_set = set(df[df["defense_moderate"]].index)
        is_cand  = np.array([pid in cand_set for pid in prot_ids])

        # Subsample
        rng    = np.random.RandomState(RANDOM_SEED)
        n_samp = min(TSNE_N_SAMPLE, len(prot_ids))
        bg_idx   = rng.choice(len(prot_ids), n_samp, replace=False)
        cand_idx = np.where(is_cand)[0]
        use_idx  = np.unique(np.concatenate([bg_idx, cand_idx]))

        combined = np.vstack([prot_emb[use_idx], anch_emb])
        n_prot   = len(use_idx)

        logger.info(f"  {sp['label']}: t-SNE on {combined.shape[0]} points "
                    f"({n_prot} proteome + {len(anch_ids)} anchors)")

        # PCA → t-SNE
        pca    = PCA(n_components=min(50, combined.shape[0] - 1),
                     random_state=RANDOM_SEED)
        coords = TSNE(n_components=2, perplexity=TSNE_PERPLEXITY,
                      random_state=RANDOM_SEED).fit_transform(
                          pca.fit_transform(combined))

        is_cand_sub = is_cand[use_idx]
        bg_mask     = ~is_cand_sub

        # Plot background
        ax.scatter(coords[:n_prot][bg_mask, 0],
                   coords[:n_prot][bg_mask, 1],
                   c="#E0E0E0", s=4, alpha=0.30,
                   rasterized=True, zorder=1)
        # Plot candidates
        ax.scatter(coords[:n_prot][is_cand_sub, 0],
                   coords[:n_prot][is_cand_sub, 1],
                   c="#FF7043", s=16, alpha=0.55,
                   edgecolors="#BF360C", linewidths=0.4,
                   rasterized=True, zorder=2)

        # Plot anchors
        meta     = pd.read_csv(meta_path)
        id_to_cat = dict(zip(meta["uniprot_id"], meta["category"]))
        cats_plotted = set()
        for i, aid in enumerate(anch_ids):
            cat   = id_to_cat.get(aid, "unknown")
            style = ANCHOR_STYLES.get(cat,
                        {"marker": "o", "color": "black", "size": 100})
            ax.scatter(coords[n_prot + i, 0], coords[n_prot + i, 1],
                       c=style["color"], s=style["size"],
                       marker=style["marker"],
                       edgecolors="black", linewidths=0.8, zorder=10)
            cats_plotted.add(cat)

        ax.set_xlabel("t-SNE 1")
        # Only show y-label on first panel (leftmost)
        if sp["panel"] == "A":
            ax.set_ylabel("t-SNE 2")
        else:
            ax.set_ylabel("")
        italic_title = "$\\it{" + sp["label"].replace(" ", "\\ ") + "}$"
        ax.set_title(italic_title)
        panel_label(ax, sp["panel"])

        # Build legend handles once (same anchors across all species)
        if legend_handles is None:
            from matplotlib.lines import Line2D
            legend_handles = [
                Line2D([], [], marker="o", color="w",
                       markerfacecolor="#FF7043",
                       markeredgecolor="#BF360C", markeredgewidth=0.5,
                       markersize=8, label="Defense candidates",
                       linestyle="None"),
                Line2D([], [], marker="o", color="w",
                       markerfacecolor="#E0E0E0", markersize=8,
                       label="Other proteins", linestyle="None"),
            ]
            for cat in sorted(cats_plotted):
                style = ANCHOR_STYLES.get(cat,
                            {"marker": "o", "color": "black"})
                legend_handles.append(
                    Line2D([], [], marker=style["marker"], color="w",
                           markerfacecolor=style["color"],
                           markeredgecolor="black", markeredgewidth=0.6,
                           markersize=11,
                           label=f"Anchor: {CAT_LABELS.get(cat, cat)}",
                           linestyle="None"))

    # Shared legend on the right of the last panel
    axes[-1].legend(handles=legend_handles, fontsize=16,
                    framealpha=0.95, loc="center left",
                    bbox_to_anchor=(1.02, 0.5),
                    title="Protein class", title_fontsize=18,
                    handletextpad=0.8, labelspacing=0.9,
                    frameon=True, edgecolor="#BDBDBD")

    fig.suptitle("ESM-2 embedding space — defense gene candidates",
                 fontsize=18, y=1.01)
    fig.tight_layout()

    path = OUT_DIR / "fig4_tsne_embedding_space.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"  Saved {path.name}")


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
