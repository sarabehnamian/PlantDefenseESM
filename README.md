PlantDefenseESM
===============

Computational pipeline for discovering and characterizing plant defense genes using large-scale protein language models (ESM-2) and curated defense-related anchors. The project takes a complete plant proteome, embeds all proteins into ESM-2 representation space, scores their similarity to known defense proteins, and returns ranked candidate defense genes together with rich summary statistics and publication-quality figures.

Project overview
----------------

The repository is organized as a **stepwise, reproducible pipeline**. Each numbered script expects to be run from the project root and writes its outputs into a corresponding subdirectory under `results/`:

- **`00_download_proteome.py`**: Download or load the target plant proteome from NCBI (or from a user-specified FASTA), filter low-quality sequences, and produce:
  - `results/<species>/00_download_proteome/proteome.fasta`
  - `results/<species>/00_download_proteome/proteome_stats.csv`
  - `results/<species>/00_download_proteome/summary.yaml`
- **`01_fetch_anchors.py`**: Fetch experimentally validated defense "anchor" proteins from UniProt, using curated categories such as NBS-LRR, PR proteins, RLK defense, signaling, antimicrobial, and HR/cell death. Outputs:
  - `results/<species>/01_fetch_anchors/anchors.fasta`
  - `results/<species>/01_fetch_anchors/anchors_metadata.csv`
  - `results/<species>/01_fetch_anchors/summary.yaml`
- **`02_embed_proteins.py`**: Use an ESM-2 model (default `esm2_t33_650M_UR50D`) to embed both the full proteome and the anchor set:
  - `results/<species>/02_embed_proteins/proteome_embeddings.npz`
  - `results/<species>/02_embed_proteins/anchor_embeddings.npz`
  - `results/<species>/02_embed_proteins/embedding_stats.yaml`
- **`03_classify_defense.py`**: Score each protein against defense categories using embedding similarity, build per-category Z-scores, and classify defense candidates at multiple thresholds (strict / moderate / lenient). Outputs:
  - `results/<species>/03_classify_defense/similarity_matrix.csv`
  - `results/<species>/03_classify_defense/zscore_matrix.csv`
  - `results/<species>/03_classify_defense/top_per_category.csv`
  - `results/<species>/03_classify_defense/summary.yaml`
- **`04_validate_annotations.py`**: Integrate keyword-based annotation signals from protein descriptions and other metadata to validate candidates and assign novelty labels (known defense / novel candidate / keyword-only). Outputs:
  - `results/<species>/04_validate_annotations/validated_results.csv`
  - `results/<species>/04_validate_annotations/enrichment_tests.csv`
  - `results/<species>/04_validate_annotations/novelty_breakdown.csv`
  - `results/<species>/04_validate_annotations/summary.yaml`
- **`05_extract_candidates.py`**: Extract final candidate tables at different stringency levels and generate FASTA files with the candidate sequences:
  - `results/<species>/05_extract_candidates/candidates_strict.csv`
  - `results/<species>/05_extract_candidates/candidates_moderate.csv`
  - `results/<species>/05_extract_candidates/candidates_lenient.csv`
  - `results/<species>/05_extract_candidates/candidate_sequences.fasta`
  - `results/<species>/05_extract_candidates/novel_candidates.csv`
  - `results/<species>/05_extract_candidates/summary.yaml`
- **`06_multispecies_figures.py`**: Build multi-panel, publication-style figures comparing defense candidate landscapes across multiple species (e.g. *Arabidopsis thaliana*, *Vitis vinifera*, *Oryza sativa*). Outputs in:
  - `results/06b_multispecies_figures/fig*_*.png`
- **`07_cross_species_compare.py`**: Perform cross-species comparison of defense categories and novelty patterns, and produce summary CSVs and logs in:
  - `results/07_cross_species_compare/`

The shared utilities used throughout the pipeline live in:

- **`shared.py`**:
  - `load_config()` – read `config.yaml` and provide defaults (species, model, thresholds, device, etc.).
  - `step_dir()` – create and return per-step result directories under `results/`.
  - `get_logger()` – per-step console + file logging.
  - `DEFENSE_CATEGORIES`, `CORE_ANCHORS`, `VALIDATION_KEYWORDS` – curated biological constants for plant immunity.

Installation
------------

This project is written in **Python 3.9+** and depends on standard scientific and machine learning libraries plus **ESM-2** from Meta AI.

### 1. Create and activate an environment

Using `conda` (recommended):

```bash
conda create -n plant-defense-esm python=3.10 -y
conda activate plant-defense-esm
```

Or with `venv`:

```bash
python -m venv .venv
.\.venv\Scripts\activate  # on Windows
source .venv/bin/activate  # on Linux/macOS
```

### 2. Install Python dependencies

Install core scientific stack, plotting libraries, and BioPython:

```bash
pip install numpy pandas matplotlib seaborn scikit-learn biopython pyyaml
```

Then install ESM-2 and its dependencies (PyTorch, fair-esm or equivalent). For example:

```bash
pip install torch --index-url https://download.pytorch.org/whl/cu121  # pick the right CUDA/CPU build
pip install fair-esm
```

If you are running entirely on CPU, use the CPU-only PyTorch wheels instead. GPU is strongly recommended for the embedding step (`02_embed_proteins.py`).

Configuration
-------------

Pipeline behaviour is controlled via a YAML configuration file at the project root:

- **`config.yaml`** (expected by `shared.load_config()`), with keys such as:
  - `species`: target species identifier (e.g. `vitis_vinifera`, `arabidopsis_thaliana`, `oryza_sativa`).
  - `proteome_path`: optional path to an existing proteome FASTA; if unset, the proteome is downloaded from NCBI.
  - `esm_model`: ESM-2 model name (default `esm2_t33_650M_UR50D`).
  - `batch_size`: batch size for embedding.
  - `max_seq_len`: maximum sequence length (ESM-2 limit, default 1022).
  - `device`: `"auto"`, `"cpu"`, or a specific CUDA device string.
  - `z_threshold_strict`, `z_threshold_moderate`, `z_threshold_lenient`: Z-score cutoffs for candidate calling.
  - `tsne_perplexity`, `tsne_n_sample`: parameters for t-SNE visualisation.
  - `base_output_dir`: base results directory (default `"results"`).

If `config.yaml` is missing, the pipeline falls back to sensible defaults defined in `shared.py`.

Running the pipeline
--------------------

From the project root, after configuring your environment and `config.yaml`, run the steps in order:

```bash
python 00_download_proteome.py
python 01_fetch_anchors.py
python 02_embed_proteins.py
python 03_classify_defense.py
python 04_validate_annotations.py
python 05_extract_candidates.py
```

Each step logs to both the console and a per-step log file under `results/<species>/<step_name>/`.

To generate multi-species figures (assuming you have already run the single-species pipeline for each of the configured species):

```bash
python 06_multispecies_figures.py
```

To perform the final cross-species comparison and summaries:

```bash
python 07_cross_species_compare.py
```

Key outputs
-----------

- **Candidate tables**: CSV files summarising defense candidates at strict, moderate, and lenient thresholds, with per-category Z-scores and novelty labels.
- **FASTA sequences**: filtered proteomes and candidate-only FASTA files suitable for downstream structural or evolutionary analyses.
- **Summary YAML**: per-step YAML files tracking counts, distributions, and parameter settings for provenance.
- **Figures**:
  - Z-score histograms with threshold lines.
  - Category-specific heatmaps of top candidates.
  - Novelty and category breakdown barplots.
  - t-SNE visualisations of ESM-2 embedding space with anchors and candidates highlighted.
  - Multi-species comparison figures for inclusion in publications.

Reproducibility and provenance
------------------------------

- All scripts are **pure Python** and rely only on `config.yaml` and the raw proteome input, making runs easy to reproduce.
- Every step writes an explicit `summary.yaml` and a log file capturing the key parameters, random seeds, and high-level statistics.
- Anchor sets, defense categories, and validation keywords are **curated and versioned** in `shared.py`, so later changes can be tracked via git history.

Citation
--------

If you use this pipeline or derived figures in a publication, please cite this GitHub repository:

> PlantDefenseESM – ESM-2–based pipeline for discovery and comparison of plant defense genes across species.  
> GitHub: https://github.com/sarabehnamian/PlantDefenseESM

You may also wish to cite:

- Meta AI's ESM-2 protein language model.
- Appropriate plant defense biology references describing the anchor genes and categories used here.

License
-------

Specify your preferred license here (e.g. MIT, Apache-2.0, or a more restrictive academic license) and add a corresponding `LICENSE` file in the repository root.