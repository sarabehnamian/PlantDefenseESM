# PlantDefenseESM

Computational pipeline for **prioritizing candidate plant defense proteins** using
protein language model embeddings (ESM-2) and a curated set of defense anchor
proteins. The pipeline takes a reference plant proteome, embeds every protein into
ESM-2 representation space, scores each protein against six defense category
centroids, and returns ranked candidate lists together with summary statistics,
benchmarks and publication-size figures.

## Scope and limitations

Please read this section before using the outputs.

- **PlantDefenseESM prioritizes candidates. It does not confirm defense function.**
  Every output is a ranked hypothesis for experimental follow-up, not a defense
  gene.
- **RefSeq keyword matching provides annotation-based support and enrichment, not
  biological confirmation.** A protein whose description contains a
  defense-associated keyword is an *annotation-supported candidate*; that is a
  statement about the text annotation, not about the protein.
- **ESM-2 did not outperform the sequence-based baselines.** On the curated
  *Arabidopsis* GO defense set, anchor alignment recovered more curated proteins
  than the embeddings at matched candidate counts (674 versus 482), and InterPro
  domain retrieval was more precise though narrower. The embeddings are
  complementary rather than superior: 70 curated proteins were recovered by ESM-2
  alone, raising recovery by the two sequence-based methods combined from 703 to
  773.
- **Curated benchmarking was feasible only in *Arabidopsis*.** Fewer than 20
  curated defense proteins could be mapped to the rice and grapevine RefSeq
  proteomes, so application to three proteomes demonstrates computational
  portability but does not establish equivalent biological performance across
  species.
- **The anchor set is Arabidopsis-heavy** (33 proteins, mostly *A. thaliana*), and
  ESM-2 truncates sequences longer than 1,022 residues, which affects the longest
  immune receptors.

### Levels of evidence

These five terms are used consistently throughout the code, the outputs and the
manuscript, and are not interchangeable:

| Term | Meaning |
| --- | --- |
| Curated benchmark protein | Included in a curated Gene Ontology or structural-family reference set |
| Annotation-supported candidate | Carries a defense-associated RefSeq keyword |
| Domain-supported keyword-negative candidate | No such keyword, but carries a designated defense-associated domain |
| Candidate lacking defense-associated annotation | Neither |
| Experimentally validated defense protein | Supported by direct experimental evidence in the literature |

Only the anchor set and the curated reference sets rest on curated or experimental
evidence.

## Project overview

The repository is a **stepwise, reproducible pipeline**. Each numbered script is
run from the project root and writes into a subdirectory of `results/`.

- **`00_download_proteome.py`** — download or load the proteome from NCBI, filter
  short and non-standard sequences.
  `results/<species>/00_download_proteome/{proteome.fasta, proteome_stats.csv, summary.yaml}`

- **`01_fetch_anchors.py`** — fetch the 33 experimentally validated anchor proteins
  from UniProt across six defense categories (NBS-LRR, PR proteins, RLK defense,
  defense signalling, antimicrobial, HR / cell death).
  `results/<species>/01_fetch_anchors/{anchors.fasta, anchors_metadata.csv, summary.yaml}`

- **`02_embed_proteins.py`** — embed the proteome and the anchors with ESM-2
  (default `esm2_t33_650M_UR50D`), mean-pooled over residues.
  `results/<species>/02_embed_proteins/{proteome_embeddings.npz, anchor_embeddings.npz, embedding_stats.yaml}`

- **`03_classify_defense.py`** — score every protein by cosine similarity to the six
  category centroids, convert to per-category Z-scores and percentile ranks, and
  select candidates at three stringency tiers.
  `results/<species>/03_classify_defense/{similarity_matrix.csv, zscore_matrix.csv, top_per_category.csv, summary.yaml}`

- **`04_validate_annotations.py`** — annotation-based support: scan RefSeq
  descriptions for the 47 defense-associated keywords, test enrichment among
  candidates (Fisher's exact test), and assign each protein an annotation-support
  class (annotation-supported candidate / candidate without defense keyword /
  keyword only / non-candidate).
  `results/<species>/04_validate_annotations/{validated_results.csv, enrichment_tests.csv, annotation_support_breakdown.csv, summary.yaml}`

- **`05_extract_candidates.py`** — candidate tables per tier plus the
  keyword-negative subset and candidate sequences.
  `results/<species>/05_extract_candidates/{candidates_strict.csv, candidates_moderate.csv, candidates_lenient.csv, keyword_negative_candidates.csv, candidate_sequences.fasta, summary.yaml}`

- **`06_multispecies_figures.py`** — manuscript Figures 2, 3, 5 and 6, drawn at
  BMC's 170 mm final size, plus the fully labelled supplementary heatmaps.
  `results/06b_multispecies_figures/`

- **`07_cross_species_compare.py`** — cross-species comparison tables and
  manuscript Figure 7.
  `results/07_cross_species_compare/`

- **`08_make_supplementary_tables.py`** — the per-protein moderate-tier candidate
  workbook (Additional file 1).
  `results/08_supplementary_tables/PlantDefenseESM_Supplementary_Candidates.xlsx`

- **`09_benchmark_curated.py`** — benchmark against the curated GO defense set
  (AUPRC, ROC-AUC, precision/recall, fold enrichment).
  `results/09_benchmark_curated/{benchmark_metrics.csv, benchmark_overall.csv, pr_curve_<species>.png}`

- **`10_benchmark_families.py`** — recall of the curated NLR / PR / LRR-RK
  structural families.
  `results/09_benchmark_curated/families_recall.csv`

- **`11_anchor_robustness.py`** — leave-one-anchor-out and leave-three-out
  robustness of the moderate-tier candidate sets.
  `results/11_anchor_robustness/{robustness_summary.csv, loao_per_anchor_<species>.csv}`

- **`12_benchmark_baselines.py`** — ESM-2 against anchor alignment
  (Smith–Waterman or BLAST), InterPro domain retrieval and RefSeq keyword
  retrieval on the curated set.
  `results/12_benchmark_baselines/{baseline_metrics.csv, recovered_breakdown.csv}`

- **`13_windowed_embed.py`** — windowed ESM-2 re-embedding of the sequences that
  exceed the 1,022-residue input limit, tiled into overlapping windows with the
  final window flush to the C-terminus.
  `results/<species>/13_windowed_embed/{windowed_proteome_embeddings.npz, windowed_anchor_embeddings.npz, windows_manifest.csv, summary.yaml}`

- **`14_truncation_sensitivity.py`** — whether retaining the C-terminal regions of
  long proteins changes the ranking and recall of long immune receptors.
  `results/14_truncation_sensitivity/`

- **`15_baseline_overlap.py`** — the complete method-overlap breakdown among the
  curated defense proteins: the eight mutually exclusive classes formed by ESM-2,
  alignment and domain retrieval, plus per-protein membership.
  `results/15_baseline_overlap/{overlap_cells_<species>.csv, overlap_pairwise_<species>.csv, overlap_table_<species>.csv, positive_membership_<species>.csv}`

Peer-review revision analyses live in **`analysis/`** and read the step outputs
above (run from the project root):

| Script | Purpose |
| --- | --- |
| `analysis/r2_param_sensitivity.py` | Percentile and top-N threshold sensitivity (Additional file 3) |
| `analysis/r2_keyword_specificity.py` | Partition of the 47 keywords into 39 defense-specific and 8 broad terms, and enrichment under each subset (Additional file 2) |
| `analysis/r2_novel_domain_support.py` | InterPro support for keyword-negative candidates |
| `analysis/r2_novel_pfam_support.py` | Pfam (pyhmmer) support for keyword-negative candidates, scanned on our own sequences |
| `analysis/r2_anchor_in_candidates.py` | Whether near-identical anchor homologs enter the candidate sets |
| `analysis/r2_stress_vs_defense.py` | Overlap of candidates with abiotic-stress, developmental and metabolic vocabulary |
| `analysis/r3_false_negatives.py` | Curated defense proteins missed at the moderate tier, with named examples |
| `analysis/r3_truncation_analysis.py` | Truncation counts per proteome, per anchor and per category |
| `analysis/r3_percategory_topcandidates.py` | Manuscript Figure 4 and its supplementary version |
| `analysis/r3_word_boundary_check.py` | Substring versus whole-word keyword matching |
| `analysis/r3_solve_keyword_split.py` | Reconstructs the broad/specific keyword split from the published fold enrichments |

Shared utilities live in **`shared.py`**: `load_config()`, `step_dir()`,
`get_logger()`, and the curated constants `DEFENSE_CATEGORIES`, `CORE_ANCHORS`
and `VALIDATION_KEYWORDS`.

## Installation

Python 3.9+, plus ESM-2 from Meta AI.

```
conda create -n plant-defense-esm python=3.10 -y
conda activate plant-defense-esm
```

or

```
python -m venv .venv
.\.venv\Scripts\activate     # Windows
source .venv/bin/activate    # Linux / macOS
```

```
pip install numpy pandas matplotlib seaborn scikit-learn biopython pyyaml openpyxl
pip install torch --index-url https://download.pytorch.org/whl/cu121   # pick your CUDA/CPU build
pip install fair-esm
pip install pyhmmer                                                    # only for r2_novel_pfam_support.py
```

The exact versions used for the published results are pinned in
`requirements.txt`:

```
pip install -r requirements.txt
```

GPU is strongly recommended for `02_embed_proteins.py` and `13_windowed_embed.py`;
everything else runs on CPU from cached outputs.

## Configuration

Behaviour is controlled by `config.yaml` at the project root, with per-species
presets provided (`config_arabidopsis.yaml`, `config_rice.yaml`,
`config_grapevine.yaml`). Keys include `species`, `proteome_path`, `esm_model`,
`batch_size`, `max_seq_len`, `device`, the percentile cutoffs
(`percentile_strict` 99.5, `percentile_moderate` 99.0, `percentile_lenient` 97.0),
`top_n_per_category` (50), `tsne_perplexity`, `tsne_n_sample`, `random_seed` (42)
and `base_output_dir`. If `config.yaml` is missing, defaults in `shared.py` apply.

## Running the pipeline

```
# 1) select a species
copy config_arabidopsis.yaml config.yaml      # Windows  (cp on Linux/macOS)

# 2) per-species pipeline
python 00_download_proteome.py
python 01_fetch_anchors.py
python 02_embed_proteins.py          # GPU recommended
python 03_classify_defense.py
python 04_validate_annotations.py
python 05_extract_candidates.py

# 3) once all three species are processed (these iterate over species internally)
python 06_multispecies_figures.py
python 07_cross_species_compare.py --runs results\arabidopsis_thaliana results\vitis_vinifera results\oryza_sativa --labels "A. thaliana" "V. vinifera" "O. sativa"
python 08_make_supplementary_tables.py
python 09_benchmark_curated.py
python 10_benchmark_families.py
python 11_anchor_robustness.py
python 12_benchmark_baselines.py
python 15_baseline_overlap.py

# 4) long-sequence re-embedding (per species, GPU)
python 13_windowed_embed.py --config config_arabidopsis.yaml
python 14_truncation_sensitivity.py

# 5) revision analyses
python analysis/r2_param_sensitivity.py
python analysis/r2_keyword_specificity.py
python analysis/r3_percategory_topcandidates.py
python analysis/r3_false_negatives.py
python analysis/r3_truncation_analysis.py
```

Every step writes a `summary.yaml` and a log file recording parameters, the random
seed and high-level statistics.

## Which script produces each manuscript item

| Manuscript item | Script | Output |
| --- | --- | --- |
| Table 2 (RefSeq keyword enrichment by tier) | `04_validate_annotations.py` | `results/<species>/04_validate_annotations/enrichment_tests.csv` |
| Table 3 (anchor robustness) | `11_anchor_robustness.py` | `results/11_anchor_robustness/robustness_summary.csv` |
| Table 4 (benchmark vs curated GO set) | `09_benchmark_curated.py` | `results/09_benchmark_curated/benchmark_metrics.csv`, `benchmark_overall.csv` |
| Table 5 (family recall) | `10_benchmark_families.py` | `results/09_benchmark_curated/families_recall.csv` |
| Table 6 (method-overlap breakdown) | `15_baseline_overlap.py` | `results/15_baseline_overlap/overlap_table_arabidopsis_thaliana.csv` |
| Table 7 (windowed embedding) | `13_windowed_embed.py` → `14_truncation_sensitivity.py` | `results/14_truncation_sensitivity/` |
| Table 8 (cross-species summary) | `07_cross_species_compare.py` | `results/07_cross_species_compare/cross_species_summary.csv` |
| Figure 2 (Z-score distributions) | `06_multispecies_figures.py` | `results/06b_multispecies_figures/fig1_zscore_distributions.png` |
| Figure 3 (category heatmaps) | `06_multispecies_figures.py` | `results/06b_multispecies_figures/fig2_category_heatmaps.png` |
| Figure 4 (per-category candidates) | `analysis/r3_percategory_topcandidates.py` | `results/r3_percategory/fig_percategory_topcandidates.png` |
| Figure 5 (annotation support and category breakdown) | `06_multispecies_figures.py` | `results/06b_multispecies_figures/fig3_annotation_support_and_breakdown.png` |
| Figure 6 (t-SNE) | `06_multispecies_figures.py` | `results/06b_multispecies_figures/fig4_tsne_embedding_space.png` |
| Figure 7 (cross-species categories) | `07_cross_species_compare.py` | `results/07_cross_species_compare/fig_species_comparison.png` |
| Additional file 1 (candidate workbook) | `08_make_supplementary_tables.py` | `results/08_supplementary_tables/PlantDefenseESM_Supplementary_Candidates.xlsx` |
| Additional file 2 (keyword partition) | `analysis/r2_keyword_specificity.py` | `results/r2_keyword_groups.csv` |
| Additional file 3 (threshold sensitivity) | `analysis/r2_param_sensitivity.py` | `results/r2_param_sensitivity_*.csv` |
| Additional file 4 (supplementary heatmaps) | `06_multispecies_figures.py`, `analysis/r3_percategory_topcandidates.py` | `figS_category_heatmaps_full_*.png`, `figS_percategory_full_*.png` |

## Reproducibility

All stochastic steps use a fixed seed (`random_seed: 42`), so a clean run
reproduces the published tables and figures.

**Input proteomes** (downloaded by `00_download_proteome.py`; URLs in that script):

| Species | RefSeq assembly | Proteins (after the ≥30-aa filter) |
| --- | --- | --- |
| *Arabidopsis thaliana* | `GCF_000001735.4` (TAIR10.1) | 48,207 |
| *Oryza sativa* | `GCF_001433935.1` (IRGSP-1.0) | 42,575 |
| *Vitis vinifera* | `GCF_030704535.1` (ASM3070453v1) | 40,632 |

**Fixed parameters**: ESM-2 checkpoint `esm2_t33_650M_UR50D` (33 layers, 1,280-d),
mean-pooled per-protein embedding excluding the BOS token; cosine similarity to
per-category centroids → per-category Z-scores → percentile ranks; tiers strict
≥ 99.5th, moderate ≥ 99.0th, lenient ≥ 97.0th percentile with
`top_n_per_category = 50`; `max_seq_len = 1022`; minimum length 30 aa;
`batch_size = 4`; t-SNE `perplexity = 30`, `n_sample = 5000`.

Small output tables and figures backing the manuscript are committed under
`results/`. The per-species ESM-2 embeddings
(`results/<species>/02_embed_proteins/proteome_embeddings.npz`, 165–209 MB each)
are regenerable and are not tracked.

## Manuscript version

This branch corresponds to the **second revision** of the associated manuscript.
A tagged release and archival DOI will be deposited on acceptance so that the
computational materials cannot change after publication without a record.

## License

MIT. See `LICENSE`.

## Associated publication

> **Behnamian, S., & Boyouk, N.** *Protein language model embeddings for
> proteome-wide prioritization of plant defense gene candidates across species.*
> Globe Institute, University of Copenhagen, Øster Voldgade 5–7, 1350 Copenhagen K,
> Denmark (`sara.behnamian@sund.ku.dk`); Kempten University of Applied Sciences,
> Bahnhofstraße 61, 87435 Kempten, Germany
> (`naghmeh.boyouk@stud.hs-kempten.de`).

**Abstract.** Identifying the full complement of defense genes across plant
proteomes remains challenging, particularly for species with incomplete functional
annotations. Here we present PlantDefenseESM, a computational pipeline that
leverages protein language model embeddings to identify and prioritize candidate
defense proteins at the proteome scale without requiring species-specific training
or curated gene ontology databases. We generated 1,280-dimensional embeddings for
all proteins in the proteomes of *Arabidopsis thaliana* (48,207 proteins),
*Oryza sativa* (42,575) and *Vitis vinifera* (40,632) using ESM-2, a
transformer-based model pre-trained on 250 million protein sequences. Defense
candidates were identified by cosine similarity to category centroids defined by
33 experimentally validated anchor proteins spanning six functional classes. A
multi-tier selection strategy identified 2,807, 2,442 and 2,354 moderate-tier
candidates respectively. Annotation-based enrichment analysis against RefSeq
descriptions showed 3.35–4.22-fold enrichment of defense-annotated proteins among
candidates (Fisher's exact test, *p* < 10⁻¹⁹⁹ in all species), and 55–59% of
candidates lacked defense-associated RefSeq keywords. Cross-species comparison
revealed a conserved category hierarchy with lineage-specific expansions
consistent with known biology. The pipeline requires only a reference proteome as
input and provides a scalable framework for defense candidate prioritization,
though its performance may still depend on the Arabidopsis-heavy anchor set and on
taxonomic biases in the ESM-2 training data.

**Keywords.** protein language model; ESM-2; plant innate immunity; proteome-wide
classification; NBS-LRR; defense candidate prioritization; *Arabidopsis thaliana*;
*Oryza sativa*; *Vitis vinifera*
