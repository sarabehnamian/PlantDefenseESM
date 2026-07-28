#!/usr/bin/env python3
"""
07_cross_species_compare.py
===========================
Compare defense gene repertoires across species by loading pre-computed
results from multiple pipeline runs.

Expects that you have run steps 00-05 once per species, each with a
different base_output_dir (e.g., results_arabidopsis, results_vitis, ...).

Usage:
    python 07_cross_species_compare.py \
        --runs results_arabidopsis results_vitis results_tomato \
        --labels arabidopsis vitis tomato

Outputs -> results/07_cross_species_compare/
    cross_species_summary.csv   - side-by-side candidate counts
    category_comparison.csv     - per-category counts across species
    fig_species_comparison.png  - grouped bar chart
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yaml
try:
    from yaml import FullLoader as YAMLLoader
except ImportError:
    from yaml import SafeLoader as YAMLLoader

import logging

from shared import get_logger


# ══════════════════════════════════════════════════════════════════════════════
# Publication style — BMC Bioinformatics
# ══════════════════════════════════════════════════════════════════════════════
# 85 mm half-page / 170 mm full-page width, max 225 mm high, ~300 dpi at final
# size, legible at the 600 px web width, lines above 0.25 pt, fonts embedded,
# Arial or Helvetica in the graphic, and no title inside the image (BMC wants
# figure titles and legends in the manuscript). The figure is drawn AT FINAL
# SIZE, so the point sizes below are the ones that appear in print.

MM = 1 / 25.4
FULL_W = 170 * MM          # 6.69 in
HALF_W = 85 * MM           # 3.35 in

SAVE_PDF = False           # True also writes the vector PDF BMC prefers


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
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}

# Same order as the other figures, so the categories line up across the paper.
CAT_ORDER = ["NBS_LRR", "PR_proteins", "RLK_defense",
             "defense_signaling", "antimicrobial", "cell_death_HR"]

CAT_LABELS_2LINE = {
    "NBS_LRR": "NBS-LRR",
    "PR_proteins": "PR\nproteins",
    "RLK_defense": "RLK\ndefense",
    "defense_signaling": "Defense\nsignalling",
    "antimicrobial": "Anti-\nmicrobial",
    "cell_death_HR": "Cell death\n/ HR",
}


def _quiet_fonttools():
    for name in ("fontTools", "fontTools.subset", "fontTools.ttLib"):
        logging.getLogger(name).setLevel(logging.WARNING)



def load_run_summary(run_dir: Path) -> dict:
    """Load summary.yaml files from a completed pipeline run."""
    info = {}

    # Classification summary
    cls_yaml = run_dir / "03_classify_defense" / "summary.yaml"
    if cls_yaml.exists():
        with open(cls_yaml) as fh:
            info["classification"] = yaml.load(fh, Loader=yaml.UnsafeLoader)

    # Validation summary
    val_yaml = run_dir / "04_validate_annotations" / "summary.yaml"
    if val_yaml.exists():
        with open(val_yaml) as fh:
            info["validation"] = yaml.load(fh, Loader=yaml.UnsafeLoader)

    # Candidate counts
    cand_yaml = run_dir / "05_extract_candidates" / "summary.yaml"
    if cand_yaml.exists():
        with open(cand_yaml) as fh:
            info["candidates"] = yaml.load(fh, Loader=yaml.UnsafeLoader)

    return info


def main():
    parser = argparse.ArgumentParser(
        description="Compare defense gene repertoires across species"
    )
    parser.add_argument("--runs", nargs="+", required=True,
                        help="Paths to completed pipeline result directories")
    parser.add_argument("--labels", nargs="+", required=True,
                        help="Species labels (same order as --runs)")
    parser.add_argument("--output-dir", default="results/07_cross_species_compare")
    args = parser.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    logger = get_logger("07_cross_species_compare", out)

    logger.info("=" * 65)
    logger.info("STEP 07  - Cross-species comparison")
    logger.info("=" * 65)

    if len(args.runs) != len(args.labels):
        raise ValueError("--runs and --labels must have the same length")

    rows = []
    cat_rows = []

    for run_path, label in zip(args.runs, args.labels):
        run_dir = Path(run_path)
        logger.info(f"\n{label}: {run_dir}")

        info = load_run_summary(run_dir)

        # Side-by-side summary
        cls = info.get("classification", {})
        val = info.get("validation", {})
        cand = info.get("candidates", {})

        row = {"species": label, "n_proteins": cls.get("n_proteins", "?")}
        for thr in ("strict", "moderate", "lenient"):
            t_info = cls.get("combined", {}).get(thr, {})
            row[f"{thr}_candidates"] = t_info.get("n_candidates", "?")
            row[f"{thr}_pct"] = t_info.get("pct", "?")

        annotation_support = val.get("annotation_support_moderate", {})
        row["annotation_supported_candidate"] = annotation_support.get("annotation_supported_candidate", "?")
        row["keyword_negative_candidates"] = annotation_support.get("candidate_without_defense_keyword", "?")
        rows.append(row)

        # Category breakdown
        cat_break = cls.get("category_breakdown_moderate", {})
        for cat, n in cat_break.items():
            cat_rows.append({
                "species": label,
                "category": cat,
                "n_candidates": n,
            })

    # Save tables
    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(out / "cross_species_summary.csv", index=False)
    logger.info(f"\n{summary_df.to_string(index=False)}")

    cat_df = pd.DataFrame(cat_rows)
    if not cat_df.empty:
        cat_df.to_csv(out / "category_comparison.csv", index=False)

        # ── Figure (manuscript Figure 7) ─────────────────────────────
        _quiet_fonttools()

        scientific_names = {
            "Arabidopsis": "Arabidopsis thaliana",
            "Rice": "Oryza sativa",
            "Grapevine": "Vitis vinifera",
            "A. thaliana": "Arabidopsis thaliana",
            "V. vinifera": "Vitis vinifera",
            "O. sativa": "Oryza sativa",
            "Arabidopsis thaliana": "Arabidopsis thaliana",
            "Oryza sativa": "Oryza sativa",
            "Vitis vinifera": "Vitis vinifera",
        }
        species_colors = {
            "Arabidopsis thaliana": "#1B9E77",
            "Oryza sativa": "#D95F02",
            "Vitis vinifera": "#7570B3",
        }

        cat_df_plot = cat_df.copy()
        cat_df_plot["species"] = cat_df_plot["species"].map(
            lambda x: scientific_names.get(x, x))

        pivot = cat_df_plot.pivot_table(index="category", columns="species",
                                        values="n_candidates", fill_value=0)

        col_order = [scientific_names.get(lbl, lbl) for lbl in args.labels
                     if scientific_names.get(lbl, lbl) in pivot.columns]
        pivot = pivot[col_order]
        colors = [species_colors.get(sp, "#666666") for sp in col_order]

        row_order = [c for c in CAT_ORDER if c in pivot.index]
        row_order += [c for c in pivot.index if c not in row_order]
        pivot = pivot.loc[row_order]
        labels = [CAT_LABELS_2LINE.get(c, c) for c in pivot.index]

        with plt.rc_context(PUB_RC):
            fig, ax = plt.subplots(figsize=(FULL_W, 2.7))

            n_cat, n_sp = pivot.shape
            x = np.arange(n_cat)
            width = 0.8 / n_sp
            vmax = float(pivot.values.max())

            for j, sp in enumerate(col_order):
                offs = (j - (n_sp - 1) / 2) * width
                vals = pivot[sp].values
                ax.bar(x + offs, vals, width=width * 0.92, color=colors[j],
                       edgecolor="white", linewidth=0.4,
                       label="$\\it{" + sp.replace(" ", "\\ ") + "}$")
                for xi, v in zip(x + offs, vals):
                    ax.text(xi, v + vmax * 0.02, f"{int(v):,}", ha="center",
                            va="bottom", fontsize=5.8, rotation=90)

            ax.set_xticks(x)
            ax.set_xticklabels(labels, linespacing=1.15)
            ax.set_ylim(0, vmax * 1.30)
            ax.set_ylabel("Moderate-tier candidates", labelpad=2)
            ax.yaxis.set_major_locator(plt.MaxNLocator(5))
            ax.yaxis.set_major_formatter(
                plt.FuncFormatter(lambda v, _: f"{v:,.0f}"))
            ax.tick_params(pad=1.5)
            # No title inside the graphic: BMC wants it in the manuscript.
            ax.legend(loc="upper right", frameon=False, handlelength=1.1,
                      handletextpad=0.5, borderpad=0.1, labelspacing=0.3,
                      ncol=3, columnspacing=1.0)

            fig.tight_layout()
            png = out / "fig_species_comparison.png"
            fig.savefig(png)
            written = [png.name]
            if SAVE_PDF:
                pdf = out / "fig_species_comparison.pdf"
                fig.savefig(pdf)
                written.append(pdf.name)
            plt.close(fig)

        logger.info(f"  Saved {' + '.join(written)}  "
                    f"({FULL_W * 25.4:.0f} mm wide)")

    logger.info(f"\nOutput -> {out}")
    logger.info("Done OK")


if __name__ == "__main__":
    main()
