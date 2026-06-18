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

- **`08_make_supplementary_tables.py`**: Assemble the per-protein moderate-tier candidate table for all three species into the supplementary workbook (Additional file 1) — one sheet per species plus a data dictionary and summary:
  - `results/08_supplementary_tables/PlantDefenseESM_Supplementary_Candidates.xlsx`
- **`09_benchmark_curated.py`**: Benchmark the ESM-2 ranking against curated GO defense sets (AUPRC, ROC-AUC, fold-enrichment, precision–recall curves) — backs **Table 3**:
  - `results/09_benchmark_curated/benchmark_metrics.csv`, `benchmark_overall.csv`, `pr_curve_<species>.png`
- **`10_benchmark_families.py`**: Family-level recall against curated NLR / PR / LRR-RK sets — backs **Table 4**:
  - `results/09_benchmark_curated/families_recall.csv`
- **`11_anchor_robustness.py`**: Leave-one-anchor-out and leave-three-out robustness of the moderate-tier predictions — backs **Table 5**:
  - `results/11_anchor_robustness/robustness_summary.csv`, `loao_per_anchor_<species>.csv`
- **`12_benchmark_baselines.py`**: Head-to-head comparison of ESM-2 against alignment (Smith–Waterman), InterPro domain retrieval, and RefSeq keyword retrieval, with the recovered-by-method breakdown — backs the **baseline comparison**:
  - `results/12_benchmark_baselines/baseline_metrics.csv`, `recovered_breakdown.csv`

Peer-review revision analyses live in **`analysis/`** and read the step outputs above (run from the project root):

- **`analysis/r2_param_sensitivity.py`** → `results/r2_param_sensitivity_*.csv` (percentile / top-N threshold sensitivity)
- **`analysis/r2_keyword_specificity.py`** → `results/r2_keyword_enrichment_*.csv`, `r2_keyword_specificity_summary.csv`
- **`analysis/r2_novel_domain_support.py`**, **`analysis/r2_novel_pfam_support.py`** → `results/r2_novel_*_support_*.csv` (orthogonal InterPro/Pfam support for keyword-negative candidates)
- **`analysis/r3_false_negatives.py`** → `results/r3_false_negatives_*.csv`
- **`analysis/r3_truncation_analysis.py`** → per-species 1,022-residue truncation counts and per-category enrichment (printed to console)

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

The exact package versions used to produce the published results are pinned in **`requirements.txt`** at the repository root. To reproduce that environment precisely:

```bash
pip install -r requirements.txt
```

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

Reproducibility
---------------

All stochastic steps use a fixed seed (`random_seed: 42`), so a clean run reproduces the published tables and figures exactly. The pieces needed to reproduce are:

**Software.** Python 3.9+; exact package versions pinned in `requirements.txt` (`pip install -r requirements.txt`). Embedding (`02_embed_proteins.py`) uses GPU; all revision analyses in `analysis/` run on CPU from cached step outputs.

**Input proteomes** (downloaded from NCBI RefSeq by `00_download_proteome.py`; full URLs in that script):

| Species | RefSeq assembly | Proteins (after the ≥30-aa filter) |
|---|---|---|
| *Arabidopsis thaliana* | `GCF_000001735.4` (TAIR10.1) | 48,207 |
| *Oryza sativa* | `GCF_001433935.1` (IRGSP-1.0) | 42,575 |
| *Vitis vinifera* | `GCF_030704535.1` (ASM3070453v1) | 40,632 |

**Fixed parameters** (in `config.yaml` / `shared.py`):

- ESM-2 checkpoint `esm2_t33_650M_UR50D` (33 layers, 1,280-d); per-protein embedding = mean-pool over residues, excluding the BOS token.
- Scoring: cosine similarity to per-category centroids → per-category Z-scores → percentile ranks.
- Stringency tiers: strict ≥ 99.5th percentile, moderate ≥ 99.0th, lenient ≥ 97.0th; `top_n_per_category = 50`.
- `max_seq_len = 1022` (longer proteins truncated); minimum length filter 30 aa; `batch_size = 4`.
- t-SNE: `perplexity = 30`, `n_sample = 5000`.

**Which script produces each table / figure:**

| Manuscript item | Script | Output file(s) |
|---|---|---|
| Table 3 (benchmark vs curated GO) | `09_benchmark_curated.py` | `results/09_benchmark_curated/benchmark_metrics.csv`, `benchmark_overall.csv` |
| Table 4 (family recall) | `10_benchmark_families.py` | `results/09_benchmark_curated/families_recall.csv` |
| Table 5 (anchor robustness) | `11_anchor_robustness.py` | `results/11_anchor_robustness/robustness_summary.csv` |
| Baseline comparison | `12_benchmark_baselines.py` | `results/12_benchmark_baselines/baseline_metrics.csv`, `recovered_breakdown.csv` |
| Additional file 1 (candidate table) | `08_make_supplementary_tables.py` | `PlantDefenseESM_Supplementary_Candidates.xlsx` |
| Figures 1–4 | `06_multispecies_figures.py` | `results/06b_multispecies_figures/fig1–4_*.png` |
| Cross-species comparison | `07_cross_species_compare.py` | `results/07_cross_species_compare/category_comparison.csv`, `fig_species_comparison.png` |
| Threshold sensitivity | `analysis/r2_param_sensitivity.py` | `results/r2_param_sensitivity_*.csv` |
| Benchmark false negatives | `analysis/r3_false_negatives.py` | `results/r3_false_negatives_*.csv` |
| Truncation analysis | `analysis/r3_truncation_analysis.py` | (printed) |

Small output tables and figures backing the manuscript are committed under `results/`; the large per-species ESM-2 embeddings (`results/<species>/02_embed_proteins/proteome_embeddings.npz`, 165–209 MB each) are regenerable and are not tracked.

**Minimal end-to-end example** (one species; repeat for each, then run the cross-species/benchmark steps):

```bash
# 1) Select a species by copying one of the provided presets to config.yaml
copy config_arabidopsis.yaml config.yaml      # Windows  (cp ... on Linux/macOS)

# 2) Per-species pipeline
python 00_download_proteome.py
python 01_fetch_anchors.py
python 02_embed_proteins.py          # GPU recommended
python 03_classify_defense.py
python 04_validate_annotations.py
python 05_extract_candidates.py

# 3) After all three species are processed (these iterate over species internally)
python 06_multispecies_figures.py
python 07_cross_species_compare.py
python 08_make_supplementary_tables.py
python 09_benchmark_curated.py
python 10_benchmark_families.py
python 11_anchor_robustness.py
python 12_benchmark_baselines.py

# 4) Peer-review revision analyses
python analysis/r2_param_sensitivity.py results/arabidopsis_thaliana results/oryza_sativa results/vitis_vinifera
python analysis/r3_truncation_analysis.py
```

Every step also writes a `summary.yaml` and a log file capturing the parameters, random seed, and high-level statistics; anchor sets, defense categories, and validation keywords are curated and versioned in `shared.py` and tracked via git history.

Citation
--------

If you use this pipeline or derived figures in a publication, please cite this GitHub repository:

> PlantDefenseESM – ESM-2–based pipeline for discovery and comparison of plant defense genes across species.  
> GitHub: https://github.com/sarabehnamian/PlantDefenseESM

License
-------

This project is released under the **MIT License**.  
See the `LICENSE` file in the repository root for the full license text.

Associated publication
----------------------

If you use this pipeline in published work, please cite:

> **Behnamian, S., & Boyouk, N.**  
> *Protein language model embeddings enable proteome-wide discovery of plant defense gene networks across species.*  
> Globe Institute, University of Copenhagen, Øster Voldgade 5–7, 1350 Copenhagen K, Denmark (`sara.behnamian@sund.ku.dk`);  
> Kempten University of Applied Sciences, Bahnhofstraße 61, 87435 Kempten, Germany (`naghmeh.boyouk@stud.hs-kempten.de`).
>
> **Abstract**  
> Identifying the full complement of defense genes across plant proteomes remains challenging, particularly for species with incomplete functional annotations. Here we present PlantDefenseESM, a computational pipeline that leverages protein language model embeddings to discover defense gene networks at the proteome scale without requiring species-specific training or curated gene ontology databases. We generated 1,280-dimensional embeddings for all proteins in the proteomes of *Arabidopsis thaliana* (48,207 proteins), *Oryza sativa* (42,575), and *Vitis vinifera* (40,632) using ESM-2, a transformer-based model pre-trained on 250 million protein sequences. Defense candidates were identified by cosine similarity to category centroids defined by 33 experimentally validated anchor proteins spanning six functional classes: NBS-LRR resistance proteins, pathogenesis-related proteins, receptor-like kinases, defense signaling components, antimicrobial enzymes, and hypersensitive response regulators. A multi-tier selection strategy combining percentile-based and rank-based approaches identified 2,807, 2,442, and 2,354 moderate-tier candidates in *A. thaliana*, *O. sativa*, and *V. vinifera*, respectively. Independent validation against RefSeq functional annotations confirmed 3.35–4.22-fold enrichment of defense-annotated proteins among candidates (Fisher's exact test, p < 10⁻¹⁹⁹ in all species). Notably, 55–59% of candidates across all three species lacked any existing defense annotation, representing candidates lacking defense-associated RefSeq keywords. Cross-species comparison revealed a conserved category hierarchy with lineage-specific expansions consistent with known biology, including expanded cell death machinery in grapevine and receptor-like kinase families in rice. The pipeline is species-agnostic, requires only a reference proteome as input, and provides a scalable framework for defense gene discovery in any sequenced plant genome.
>
> **Keywords**  
> protein language model; ESM-2; plant innate immunity; defense gene discovery; proteome-wide classification; NBS-LRR; novel gene prediction; *Arabidopsis thaliana*; *Oryza sativa*; *Vitis vinifera*