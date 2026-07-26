# PlantDefenseESM

Computational pipeline for **prioritizing candidate plant defense proteins** using protein language model embeddings (ESM-2) and a curated panel of defense anchor proteins. The pipeline takes a complete plant proteome, embeds all proteins into ESM-2 representation space, scores their cosine similarity to per-category anchor centroids, and returns ranked candidate lists together with summary statistics and figures.

## Scope and limitations

**PlantDefenseESM prioritizes candidates. It does not confirm defense function.**

- Outputs are **ranked hypotheses for experimental follow-up**, not validated defense genes. No candidate produced by this pipeline has been experimentally tested here.
- A RefSeq description containing a defense-associated keyword provides **annotation-based support**, not independent biological confirmation.
- In benchmarking against a curated *Arabidopsis* defense set, ESM-2 embedding similarity **did not outperform** sequence- or domain-based methods: anchor alignment achieved a higher AUPRC (0.180 vs 0.150) and InterPro domain retrieval was substantially more precise (0.56 vs 0.17). The embeddings recovered a **complementary** subset of curated defense proteins that those baselines missed; that complementarity, not superior discrimination, is the point of the method.
- Curated benchmarking was feasible only in *A. thaliana*. Fewer than 20 curated defense proteins could be mapped to the rice and grapevine RefSeq proteomes, so **equivalent biological performance across species has not been established**. Application to three proteomes demonstrates computational portability, not cross-species accuracy.
- Performance may vary with taxonomic distance from the Arabidopsis-heavy anchor panel, proteome annotation quality, and protein-family composition.
- Proteins longer than 1,022 residues are truncated by ESM-2. This disproportionately affects NBS-LRR proteins and long LRR receptor kinases; sensitivity for these classes is uncertain.

### Terminology used in this repository and in the manuscript

| Term | Meaning |
| --- | --- |
| Curated benchmark protein | Included in a curated GO or structural-family reference set |
| Annotation-supported candidate | Candidate whose RefSeq description contains a defense-associated keyword |
| Domain-supported keyword-negative candidate | Lacks a selected keyword but contains a designated defense-associated domain |
| Candidate lacking defense-associated annotation | Lacks both selected keywords and designated domain evidence |
| Experimentally validated defense protein | Supported by direct experimental evidence in the literature (applies only to the 33 anchor proteins) |

The words *discovery*, *novel*, and *confirmed* are not used to describe pipeline outputs.

## Manuscript version

| | |
| --- | --- |
| Release corresponding to the submitted revision | `vX.Y.Z` |
| Archived DOI | `10.5281/zenodo.XXXXXXX` |
| Manuscript | Behnamian & Boyouk, submitted to *BMC Bioinformatics* (submission ID 642bb21d-50eb-413e-8f0d-48481aa415d4) |

Code, parameters, and committed result tables in that release are the ones underlying the submitted manuscript. Later commits on `main` may diverge.

## Project overview

The repository is organized as a **stepwise, reproducible pipeline**. Each numbered script expects to be run from the project root and writes its outputs into a corresponding subdirectory under `results/`:

- **`00_download_proteome.py`**: Download or load the target plant proteome from NCBI (or from a user-specified FASTA), filter low-quality sequences, and produce:

  * `results/<species>/00_download_proteome/proteome.fasta`
  * `results/<species>/00_download_proteome/proteome_stats.csv`
  * `results/<species>/00_download_proteome/summary.yaml`

- **`01_fetch_anchors.py`**: Fetch experimentally validated defense "anchor" proteins from UniProt, using curated categories such as NBS-LRR, PR proteins, RLK defense, signaling, antimicrobial, and HR/cell death. Outputs:

  * `results/<species>/01_fetch_anchors/anchors.fasta`
  * `results/<species>/01_fetch_anchors/anchors_metadata.csv`
  * `results/<species>/01_fetch_anchors/summary.yaml`

- **`02_embed_proteins.py`**: Use an ESM-2 model (default `esm2_t33_650M_UR50D`) to embed both the full proteome and the anchor set:

  * `results/<species>/02_embed_proteins/proteome_embeddings.npz`
  * `results/<species>/02_embed_proteins/anchor_embeddings.npz`
  * `results/<species>/02_embed_proteins/embedding_stats.yaml`

- **`03_classify_defense.py`**: Score each protein against defense categories using embedding similarity, build per-category Z-scores, and select candidates at multiple stringency tiers (strict / moderate / lenient). Outputs:

  * `results/<species>/03_classify_defense/similarity_matrix.csv`
  * `results/<species>/03_classify_defense/zscore_matrix.csv`
  * `results/<species>/03_classify_defense/top_per_category.csv`
  * `results/<species>/03_classify_defense/summary.yaml`

- **`04_validate_annotations.py`**: Match defense-associated keywords against RefSeq protein descriptions to measure **annotation-based support** for the candidates, and assign each protein an annotation-support class (`annotation_supported_candidate` / `candidate_without_defense_keyword` / `keyword_only` / `non_candidate`). This is an enrichment check against an orthogonal annotation source; it is not biological validation. Outputs:

  * `results/<species>/04_validate_annotations/validated_results.csv`
  * `results/<species>/04_validate_annotations/enrichment_tests.csv`
  * `results/<species>/04_validate_annotations/annotation_support_breakdown.csv`
  * `results/<species>/04_validate_annotations/summary.yaml`

- **`05_extract_candidates.py`**: Extract final candidate tables at each stringency tier and generate FASTA files with the candidate sequences:

  * `results/<species>/05_extract_candidates/candidates_strict.csv`
  * `results/<species>/05_extract_candidates/candidates_moderate.csv`
  * `results/<species>/05_extract_candidates/candidates_lenient.csv`
  * `results/<species>/05_extract_candidates/candidate_sequences.fasta`
  * `results/<species>/05_extract_candidates/keyword_negative_candidates.csv`
  * `results/<species>/05_extract_candidates/summary.yaml`

- **`06_multispecies_figures.py`**: Build multi-panel figures comparing candidate landscapes across species (*Arabidopsis thaliana*, *Vitis vinifera*, *Oryza sativa*). Outputs in:

  * `results/06b_multispecies_figures/fig*_*.png`

- **`07_cross_species_compare.py`**: Cross-species comparison of defense categories and annotation-support patterns, producing summary CSVs and logs in:

  * `results/07_cross_species_compare/`

  Note: because candidate selection is percentile-based, candidate counts are a fixed fraction of each proteome by design. Comparable proportions across species are a property of the selection rule, not evidence of similarly sized defense repertoires.

- **`08_make_supplementary_tables.py`**: Assemble the per-protein moderate-tier candidate table for all three species into the supplementary workbook (Additional file 1) — one sheet per species plus a data dictionary and summary:

  * `results/08_supplementary_tables/PlantDefenseESM_Supplementary_Candidates.xlsx`

- **`09_benchmark_curated.py`**: Benchmark the ESM-2 ranking against curated GO defense sets (AUPRC, ROC-AUC, fold enrichment, precision–recall curves):

  * `results/09_benchmark_curated/benchmark_metrics.csv`, `benchmark_overall.csv`, `pr_curve_<species>.png`

- **`10_benchmark_families.py`**: Family-level recall against curated NLR / PR / LRR-RK sets:

  * `results/09_benchmark_curated/families_recall.csv`

- **`11_anchor_robustness.py`**: Leave-one-anchor-out and leave-three-out robustness of the moderate-tier candidate sets:

  * `results/11_anchor_robustness/robustness_summary.csv`, `loao_per_anchor_<species>.csv`

- **`12_benchmark_baselines.py`**: Head-to-head comparison of ESM-2 against alignment (Smith–Waterman), InterPro domain retrieval, and RefSeq keyword retrieval, with the recovered-by-method breakdown:

  * `results/12_benchmark_baselines/baseline_metrics.csv`, `recovered_breakdown.csv`

Peer-review revision analyses live in **`analysis/`** and read the step outputs above (run from the project root):

- **`analysis/r2_param_sensitivity.py`** → `results/r2_param_sensitivity_*.csv` (percentile / top-N threshold sensitivity)
- **`analysis/r2_keyword_specificity.py`** → `results/r2_keyword_enrichment_*.csv`, `r2_keyword_specificity_summary.csv`
- **`analysis/r2_novel_domain_support.py`**, **`analysis/r2_novel_pfam_support.py`** → `results/r2_novel_*_support_*.csv` (orthogonal InterPro/Pfam support for keyword-negative candidates)
- **`analysis/r3_false_negatives.py`** → `results/r3_false_negatives_*.csv`
- **`analysis/r3_truncation_analysis.py`** → per-species 1,022-residue truncation counts and per-category enrichment (printed to console)

The shared utilities used throughout the pipeline live in:

- **`shared.py`**:
  * `load_config()` – read `config.yaml` and provide defaults (species, model, thresholds, device, etc.).
  * `step_dir()` – create and return per-step result directories under `results/`.
  * `get_logger()` – per-step console + file logging.
  * `DEFENSE_CATEGORIES`, `CORE_ANCHORS`, `VALIDATION_KEYWORDS` – curated biological constants for plant immunity.

## Installation

This project is written in **Python 3.9+** and depends on standard scientific and machine learning libraries plus **ESM-2** from Meta AI.

### 1. Create and activate an environment

Using `conda` (recommended):

```
conda create -n plant-defense-esm python=3.10 -y
conda activate plant-defense-esm
```

Or with `venv`:

```
python -m venv .venv
.\.venv\Scripts\activate  # on Windows
source .venv/bin/activate  # on Linux/macOS
```

### 2. Install Python dependencies

```
pip install numpy pandas matplotlib seaborn scikit-learn biopython pyyaml
```

Then install ESM-2 and its dependencies (PyTorch, fair-esm or equivalent):

```
pip install torch --index-url https://download.pytorch.org/whl/cu121  # pick the right CUDA/CPU build
pip install fair-esm
```

If you are running entirely on CPU, use the CPU-only PyTorch wheels instead. GPU is strongly recommended for the embedding step (`02_embed_proteins.py`).

The exact package versions used to produce the reported results are pinned in **`requirements.txt`**:

```
pip install -r requirements.txt
```

## Configuration

Pipeline behaviour is controlled via a YAML configuration file at the project root:

- **`config.yaml`** (expected by `shared.load_config()`), with keys such as:
  * `species`: target species identifier (e.g. `vitis_vinifera`, `arabidopsis_thaliana`, `oryza_sativa`).
  * `proteome_path`: optional path to an existing proteome FASTA; if unset, the proteome is downloaded from NCBI.
  * `esm_model`: ESM-2 model name (default `esm2_t33_650M_UR50D`).
  * `batch_size`: batch size for embedding.
  * `max_seq_len`: maximum sequence length (ESM-2 limit, default 1022).
  * `device`: `"auto"`, `"cpu"`, or a specific CUDA device string.
  * `z_threshold_strict`, `z_threshold_moderate`, `z_threshold_lenient`: cutoffs for candidate selection.
  * `tsne_perplexity`, `tsne_n_sample`: parameters for t-SNE visualisation.
  * `base_output_dir`: base results directory (default `"results"`).

The stringency tiers are pragmatic prioritization cutoffs, not trained classification boundaries. Because they are percentile-based, the number of candidates is a fixed fraction of each proteome by design.

If `config.yaml` is missing, the pipeline falls back to defaults defined in `shared.py`.

## Running the pipeline

From the project root, after configuring your environment and `config.yaml`, run the steps in order:

```
python 00_download_proteome.py
python 01_fetch_anchors.py
python 02_embed_proteins.py
python 03_classify_defense.py
python 04_validate_annotations.py
python 05_extract_candidates.py
```

Each step logs to both the console and a per-step log file under `results/<species>/<step_name>/`.

To generate multi-species figures (after running the single-species pipeline for each species):

```
python 06_multispecies_figures.py
```

To perform the cross-species comparison and summaries:

```
python 07_cross_species_compare.py
```

## Key outputs

- **Candidate tables**: CSV files listing candidates at strict, moderate, and lenient tiers, with per-category Z-scores and annotation-support class.
- **FASTA sequences**: filtered proteomes and candidate-only FASTA files suitable for downstream structural or evolutionary analyses.
- **Summary YAML**: per-step YAML files tracking counts, distributions, and parameter settings for provenance.
- **Figures**:
  * Z-score histograms with tier thresholds.
  * Category-specific heatmaps of top candidates.
  * Annotation-support and category breakdown barplots.
  * t-SNE visualisations of ESM-2 embedding space with anchors and candidates highlighted. Candidates are selected as high-scoring outliers, so their separation from the background sample is expected and is not independent validation.
  * Multi-species comparison figures.

## Reproducibility

All stochastic steps use a fixed seed (`random_seed: 42`), so a clean run reproduces the reported tables and figures exactly. The pieces needed to reproduce are:

**Software.** Python 3.9+; exact package versions pinned in `requirements.txt` (`pip install -r requirements.txt`). Embedding (`02_embed_proteins.py`) uses GPU; all revision analyses in `analysis/` run on CPU from cached step outputs.

**Input proteomes** (downloaded from NCBI RefSeq by `00_download_proteome.py`; full URLs in that script):

| Species                | RefSeq assembly                  | Proteins (after the ≥30-aa filter) |
| ---------------------- | -------------------------------- | ---------------------------------- |
| *Arabidopsis thaliana* | `GCF_000001735.4` (TAIR10.1)     | 48,207                             |
| *Oryza sativa*         | `GCF_001433935.1` (IRGSP-1.0)    | 42,575                             |
| *Vitis vinifera*       | `GCF_030704535.1` (ASM3070453v1) | 40,632                             |

**Fixed parameters** (in `config.yaml` / `shared.py`):

- ESM-2 checkpoint `esm2_t33_650M_UR50D` (33 layers, 1,280-d); per-protein embedding = mean-pool over residues, excluding the BOS token.
- Scoring: cosine similarity to per-category centroids → per-category Z-scores → percentile ranks.
- Stringency tiers: strict ≥ 99.5th percentile, moderate ≥ 99.0th, lenient ≥ 97.0th; `top_n_per_category = 50`.
- `max_seq_len = 1022` (longer proteins truncated); minimum length filter 30 aa; `batch_size = 4`.
- t-SNE: `perplexity = 30`, `n_sample = 5000`.

**Which script produces each table / figure** (table numbers refer to the release named under *Manuscript version* above):

| Manuscript item                     | Script                               | Output file(s)                                                                           |
| ----------------------------------- | ------------------------------------ | ---------------------------------------------------------------------------------------- |
| Table 3 (benchmark vs curated GO)   | `09_benchmark_curated.py`            | `results/09_benchmark_curated/benchmark_metrics.csv`, `benchmark_overall.csv`            |
| Table 4 (family recall)             | `10_benchmark_families.py`           | `results/09_benchmark_curated/families_recall.csv`                                       |
| Table 5 (anchor robustness)         | `11_anchor_robustness.py`            | `results/11_anchor_robustness/robustness_summary.csv`                                    |
| Table 6 (cross-species summary)     | `07_cross_species_compare.py`        | `results/07_cross_species_compare/category_comparison.csv`                                |
| Baseline comparison and overlap     | `12_benchmark_baselines.py`          | `results/12_benchmark_baselines/baseline_metrics.csv`, `recovered_breakdown.csv`         |
| Additional file 1 (candidate table) | `08_make_supplementary_tables.py`    | `PlantDefenseESM_Supplementary_Candidates.xlsx`                                          |
| Figures 1–4                         | `06_multispecies_figures.py`         | `results/06b_multispecies_figures/fig1–4_*.png`                                          |
| Figure 7 (cross-species bar chart)  | `07_cross_species_compare.py`        | `results/07_cross_species_compare/fig_species_comparison.png`                             |
| Threshold sensitivity               | `analysis/r2_param_sensitivity.py`   | `results/r2_param_sensitivity_*.csv`                                                     |
| Benchmark false negatives           | `analysis/r3_false_negatives.py`     | `results/r3_false_negatives_*.csv`                                                       |
| Truncation analysis                 | `analysis/r3_truncation_analysis.py` | (printed)                                                                                |

Small output tables and figures backing the manuscript are committed under `results/`; the large per-species ESM-2 embeddings (`results/<species>/02_embed_proteins/proteome_embeddings.npz`, 165–209 MB each) are regenerable and are not tracked.

**Minimal end-to-end example** (one species; repeat for each, then run the cross-species/benchmark steps):

```
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

Every step also writes a `summary.yaml` and a log file capturing the parameters, random seed, and high-level statistics; anchor sets, defense categories, and keyword lists are curated and versioned in `shared.py` and tracked via git history.

## Citation

If you use this pipeline or derived figures in a publication, please cite the archived release:

> Behnamian, S., & Boyouk, N. PlantDefenseESM: an ESM-2 embedding pipeline for prioritizing candidate plant defense proteins. Version `vX.Y.Z`. Zenodo. `10.5281/zenodo.XXXXXXX`
> GitHub: <https://github.com/sarabehnamian/PlantDefenseESM>

## License

This project is released under the **MIT License**.
See the `LICENSE` file in the repository root for the full license text.

## Associated publication

> **Behnamian, S., & Boyouk, N.**
> *Protein language model embeddings are less powerful than sequence homology searches for proteome-wide identification of plant defense gene networks across species.*
> Submitted to *BMC Bioinformatics*. DOI to be added on publication.
>
> Lund University, Lund, Sweden; Globe Institute, University of Copenhagen, Copenhagen, Denmark (`sara.behnamian@sund.ku.dk`); Kempten University of Applied Sciences, Kempten, Germany (`naghmeh.boyouk@stud.hs-kempten.de`).

The abstract is deliberately not reproduced here, so that this README cannot drift out of step with the published text. Please read the article for the current abstract, results, and limitations.
