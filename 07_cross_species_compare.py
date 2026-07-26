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

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yaml
try:
    from yaml import FullLoader as YAMLLoader
except ImportError:
    from yaml import SafeLoader as YAMLLoader

from shared import get_logger


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

        # ── Figure ───────────────────────────────────────────────────────
        # Set Times New Roman font
        plt.rcParams["font.family"] = "Times New Roman"
        plt.rcParams["font.size"] = 11
        
        # Map labels to scientific names (handles both short and full labels)
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
        # Species-specific colors (vibrant publication palette)
        species_colors = {
            "Arabidopsis thaliana": "#1B9E77",    # vivid teal-green
            "Oryza sativa": "#D95F02",            # bright orange
            "Vitis vinifera": "#7570B3",          # rich purple
        }
        
        # Rename species to scientific names in dataframe
        cat_df_plot = cat_df.copy()
        cat_df_plot["species"] = cat_df_plot["species"].map(
            lambda x: scientific_names.get(x, x)
        )
        
        pivot = cat_df_plot.pivot_table(
            index="category", columns="species",
            values="n_candidates", fill_value=0,
        )
        
        # Reorder columns and get colors in matching order
        col_order = [scientific_names.get(lbl, lbl) for lbl in args.labels 
                     if scientific_names.get(lbl, lbl) in pivot.columns]
        pivot = pivot[col_order]
        colors = [species_colors.get(sp, "#666666") for sp in col_order]
        
        # Format category names for display
        category_labels = {
            "NBS_LRR": "NBS-LRR",
            "PR_proteins": "PR proteins",
            "RLK_defense": "RLK defense",
            "antimicrobial": "antimicrobial",
            "cell_death_HR": "cell death/HR",
            "defense_signaling": "defense signaling",
        }
        pivot.index = [category_labels.get(cat, cat) for cat in pivot.index]
        
        fig, ax = plt.subplots(figsize=(10, 6))
        pivot.plot(kind="bar", ax=ax, color=colors, edgecolor="white", linewidth=0.5)
        
        # Publication-quality font sizes (consistent for x and y labels)
        label_fontsize = 12
        title_fontsize = 14
        tick_fontsize = 10
        
        ax.set_xlabel("Category", fontsize=label_fontsize)
        ax.set_ylabel("Number of candidates (moderate)", fontsize=label_fontsize)
        ax.set_title("Defense candidates by category across species", fontsize=title_fontsize)
        ax.tick_params(axis='both', labelsize=tick_fontsize)
        
        italic_labels = ["$\\it{" + sp.replace(" ", "\\ ") + "}$" for sp in col_order]
        ax.legend(title="Species", labels=italic_labels, 
                  loc="upper left", bbox_to_anchor=(1.02, 1), frameon=True,
                  fontsize=tick_fontsize, title_fontsize=label_fontsize)
        plt.xticks(rotation=30, ha="right")
        fig.savefig(out / "fig_species_comparison.png",
                    dpi=350, bbox_inches="tight")
        plt.close(fig)
        logger.info(f"  OK fig_species_comparison.png")

    logger.info(f"\nOutput -> {out}")
    logger.info("Done OK")


if __name__ == "__main__":
    main()
