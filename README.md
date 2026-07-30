# PlantDefenseESM

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21673071.svg)](https://doi.org/10.5281/zenodo.21673071)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

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
- **ESM-2 is less sensitive than profile-based homology search.** On the curated
  *Arabidopsis* GO defense set, at a candidate count matched exactly to the ESM-2
  moderate tier, iterative profile search (jackhmmer) beat the embeddings on every
  measure: 697 curated proteins recovered versus 482, precision 0.248 versus
  0.172, recall 0.301 versus 0.208, F1 0.272 versus 0.188 and AUPRC 0.183 versus
  0.150. InterPro domain retrieval was more precise still (0.564) though the least
  sensitive (recall 0.144). **Profile-based homology search should remain the
  primary tool; pLM scoring is an adjunct to it, not a replacement.**
- **The contribution is complementarity, not accuracy.** 65 curated defense
  proteins were recovered by ESM-2 alone — by none of the seven homology- and
  domain-based methods benchmarked here. They are not compositionally unusual
  (low-complexity content, repeat content and length are indistinguishable from
  the candidate set as a whole); what distinguishes them is family membership.
  They come from families with no representative among the 33 anchors, including
  respiratory burst oxidase homologues, SYP121/SYP122 syntaxins, copines and
  papain-family cysteine proteases, which is why an anchor-seeded homology search
  cannot reach them. Measured against narrower panels the figure is higher: 72
  against jackhmmer plus domain retrieval, and 70 against the single-sequence
  alignment baseline plus domain retrieval.
- **PSI-BLAST and filter-off searches have been run** (`17_psiblast_filteroff.py`,
  BLAST+ 2.17.0). Disabling low-complexity and composition-based corrections
  (`-seg no -comp_based_stats 0`) changed recovery by two curated proteins in
  blastp and eight in PSI-BLAST out of roughly 700, so filter settings do not
  account for the difference between the two approaches. PSI-BLAST (687 curated
  proteins recovered) agreed closely with jackhmmer (697). Single-pass blastp was
  the most sensitive arm at matched candidate count (705 recovered), although the
  iterated profile methods ranked the proteome better overall (AUPRC 0.194 versus
  0.187), so the sensitivity advantage of profile search is not realised with
  this small and heterogeneous anchor panel. The ESM-2-only count is therefore a
  measured value rather than an upper bound.
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

- **`16_profile_baselines.py`** — ESM-2 against iterative profile search
  (jackhmmer via pyhmmer: each of the 33 anchors as query, three iterations,
  per-iteration inclusion threshold E < 0.001, every protein scored by its best
  full-sequence bitscore to any anchor, candidate count matched to the ESM-2
  moderate tier) and InterPro domain retrieval. Also produces the mutually
  exclusive three-method overlap and characterises the ESM-2-only residual by
  length, truncation status, low-complexity content and repeat content. The
  PSI-BLAST and filter-off arms sketched here were superseded by
  `17_psiblast_filteroff.py`, which runs them against a local BLAST+ install.
  `results/16_profile_baselines/`

- **`17_psiblast_filteroff.py`** — the PSI-BLAST and low-complexity-filter-off
  arms requested by the editor. Requires BLAST+ 2.17.0 or later on `PATH`. Runs
  four searches of the 33 anchors against the proteome — blastp and PSI-BLAST
  (three iterations, per-iteration inclusion threshold E < 0.001), each with
  default settings and with `-seg no -comp_based_stats 0` — scores every protein
  by its best bitscore to any anchor, matches the candidate count to the ESM-2
  moderate tier, and evaluates each arm against the curated GO defense set.
  Optional phmmer arms (`--max --nonull2`) via `--with-phmmer`. Note that the
  output directory is `results/13_psiblast_filteroff/`, which predates the
  renumbering of the script.
  `results/13_psiblast_filteroff/{arm_metrics.csv, arm_scores.csv, recovery_matrix.csv, residual_esm_only.csv, residual_features.csv}`

- **`18_final_residual.py`** — folds the arms above together with the cached
  jackhmmer scores, the Smith–Waterman alignment set and InterPro domain
  retrieval into a single seven-method panel, and reports the ESM-2-only residual
  against all of them, together with the proteins that the wider panel recovers
  and the method that recovered each one.
  `results/18_final_residual/{final_metrics.csv, final_residual.csv, recovered_from_72.csv, summary.md}`

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
pip install pyhmmer                                                    # r2_novel_pfam_support.py and 16_profile_baselines.py
pip install tabulate                                                   # 17_psiblast_filteroff.py summary output
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
python 16_profile_baselines.py

# PSI-BLAST and filter-off arms (needs BLAST+ on PATH; see Installation).
# The two input lists are derived from positive_membership_<species>.csv:
#   curated_go_set.csv = every row; domain_set.csv = rows where domain == True
python 17_psiblast_filteroff.py --species arabidopsis_thaliana --curated results\15_baseline_overlap\curated_go_set.csv --domain-set results\15_baseline_overlap\domain_set.csv --threads 8
python 18_final_residual.py

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

Table 1 is a summary of prior literature and Figure 1 is a schematic overview of
the pipeline; neither is generated by a script.

| Manuscript item | Script | Output |
| --- | --- | --- |
| Table 2 (RefSeq keyword enrichment by tier) | `04_validate_annotations.py` | `results/<species>/04_validate_annotations/enrichment_tests.csv` |
| Table 3 (anchor robustness) | `11_anchor_robustness.py` | `results/11_anchor_robustness/robustness_summary.csv` |
| Table 4 (benchmark vs curated GO set) | `09_benchmark_curated.py` | `results/09_benchmark_curated/benchmark_metrics.csv`, `benchmark_overall.csv` |
| Table 5 (family recall) | `10_benchmark_families.py` | `results/09_benchmark_curated/families_recall.csv` |
| Table 6 (method-overlap breakdown) | `15_baseline_overlap.py` | `results/15_baseline_overlap/overlap_table_arabidopsis_thaliana.csv` |
| Table 7 (windowed embedding) | `13_windowed_embed.py` → `14_truncation_sensitivity.py` | `results/14_truncation_sensitivity/` |
| Table 8 (ESM-2 vs homology search vs domain retrieval) | `16_profile_baselines.py` (jackhmmer) → `17_psiblast_filteroff.py` (blastp, PSI-BLAST, filter-off arms) | `results/16_profile_baselines/`, `results/13_psiblast_filteroff/arm_metrics.csv` |
| Table 9 (three-method mutually exclusive overlap; ESM-2-only residual) | `16_profile_baselines.py`; full-panel residual from `18_final_residual.py` | `results/16_profile_baselines/`, `results/18_final_residual/summary.md` |
| Table 10 (cross-species summary) | `07_cross_species_compare.py` | `results/07_cross_species_compare/cross_species_summary.csv` |
| Figure 2 (Z-score distributions) | `06_multispecies_figures.py` | `results/06b_multispecies_figures/fig1_zscore_distributions.png` |
| Figure 3 (category heatmaps) | `06_multispecies_figures.py` | `results/06b_multispecies_figures/fig2_category_heatmaps.png` |
| Figure 4 (per-category candidates) | `analysis/r3_percategory_topcandidates.py` | `results/r3_percategory/fig_percategory_topcandidates.png` |
| Figure 5 (annotation support and category breakdown) | `06_multispecies_figures.py` | `results/06b_multispecies_figures/fig3_annotation_support_and_breakdown.png` |
| Figure 6 (t-SNE) | `06_multispecies_figures.py` | `results/06b_multispecies_figures/fig4_tsne_embedding_space.png` |
| Figure 7 (cross-species categories) | `07_cross_species_compare.py` | `results/07_cross_species_compare/fig_species_comparison.png` |
| Additional file 1 (candidate workbook) | `08_make_supplementary_tables.py` | `results/08_supplementary_tables/PlantDefenseESM_Supplementary_Candidates.xlsx` |
| Additional file 2 (keyword partition) | `analysis/r2_keyword_specificity.py` | `results/r2_keyword_groups.csv`, submitted as `PlantDefenseESM_Additional_file_2_keyword_partition.csv` |
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

This repository is versioned alongside the manuscript. **Cite the tagged release,
not the `main` branch**, which may move ahead of the published record.

| Release | Manuscript version | Zenodo DOI |
| --- | --- | --- |
| `v1.0.0` | Revision 4, first snapshot | [10.5281/zenodo.21673072](https://doi.org/10.5281/zenodo.21673072) |
| `v1.0.1` | Revision 4, submitted to *BMC Bioinformatics* (submission ID 642bb21d-50eb-413e-8f0d-48481aa415d4) | [10.5281/zenodo.21679034](https://doi.org/10.5281/zenodo.21679034) |
| `v1.0.2` | Revision 5: PSI-BLAST and low-complexity-filter-off benchmark arms (editor comments E4 and E5) | [10.5281/zenodo.21710001](https://doi.org/10.5281/zenodo.21710001) |

The concept DOI [10.5281/zenodo.21673071](https://doi.org/10.5281/zenodo.21673071)
always resolves to the most recent release; the badge above uses it.

The release is archived on Zenodo so that the computational materials cannot
change after publication without a record.

### A note on file naming

Some script and output names predate the terminology revision and are retained so
that existing result paths stay valid: `04_validate_annotations.py` and its
`validated_results.csv` perform annotation-based keyword support (not validation),
`VALIDATION_KEYWORDS` in `shared.py` is the defense-associated keyword list, and
`analysis/r2_novel_domain_support.py` / `r2_novel_pfam_support.py` operate on
keyword-negative candidates. The candidate class labels inside every output file
follow the evidence vocabulary above.

## License

MIT. See `LICENSE`.

## Associated publication

Please cite both the archived software release and the publication.

**Software:**

> Behnamian, S., & Boyouk, N. (2026). *PlantDefenseESM: a training-free pipeline
> for proteome-wide prioritization of plant defense candidates using ESM-2
> embeddings* (v1.0.2) [Software]. Zenodo.
> https://doi.org/10.5281/zenodo.21710001

**Publication:**

> **Behnamian, S.**<sup>1,2</sup>, **& Boyouk, N.**<sup>3</sup> *Protein language
> model embeddings are less sensitive than profile-based homology search for
> proteome-wide prioritization of plant defense genes across species.*
>
> <sup>1</sup> Department of Biology, Lund University, Lund 22362, Sweden
> (`sara.behnamian@biol.lu.se`)
> <sup>2</sup> Globe Institute, University of Copenhagen, Øster Voldgade 5–7,
> 1350 Copenhagen K, Denmark (`sara.behnamian@sund.ku.dk`)
> <sup>3</sup> Kempten University of Applied Sciences, Bahnhofstraße 61,
> 87435 Kempten, Germany (`naghmeh.boyouk@stud.hs-kempten.de`)

**Background.** Identifying the full complement of defense genes across plant
proteomes remains challenging, particularly for species with incomplete functional
annotations. Protein language models (pLMs) are widely proposed as an alternative
to homology-based search, but they have rarely been benchmarked against the
profile-based methods that represent the state of the art in remote homology
detection. Here we present PlantDefenseESM, a training-free pipeline that scores
every protein in a proteome by ESM-2 embedding similarity to category centroids
defined by 33 experimentally validated anchor proteins spanning six defense
classes, and we benchmark it directly against sequence homology search.

**Results.** We generated 1,280-dimensional embeddings for all proteins in the
proteomes of *Arabidopsis thaliana* (48,207 proteins), *Oryza sativa* (42,575) and
*Vitis vinifera* (40,632). A multi-tier selection strategy identified 2,807, 2,442
and 2,354 moderate-tier candidates respectively, enriched 3.35–4.22-fold for
defense-annotated proteins over the proteome background (one-sided Fisher's exact
test, *p* < 10⁻¹⁹⁹ in all species). Benchmarked against a curated Gene Ontology
defense set in *A. thaliana* at a matched candidate count, however, iterative
profile search (jackhmmer) outperformed the embeddings on every measure: precision
0.248 versus 0.172, recall 0.301 versus 0.208, F1 0.272 versus 0.188 and AUPRC
0.183 versus 0.150, while InterPro domain retrieval was more precise still
(0.564). The embeddings are therefore less sensitive than profile-based homology
search, not more. They did, however, recover 65 curated defense proteins that
no homology- or domain-based method tested identified. These were not
compositionally unusual: low-complexity content, repeat content and sequence
length were indistinguishable from the candidate set as a whole. They belonged
instead to families with no representative among the anchors, including
respiratory burst oxidase homologues, SYP121/SYP122 syntaxins, copines and
papain-family cysteine proteases.

**Conclusions.** Profile-based homology search should remain the primary tool for
identifying defense genes in plant proteomes, and pLM embeddings should not be
presented as a replacement for it. Their value here is complementary: because
centroid similarity in embedding space is not confined to the sequence families of
the anchors themselves, the embeddings surface defense-related proteins that an
anchor-seeded homology search cannot reach. We therefore recommend pLM scoring as
an adjunct to, rather than a substitute for, profile-based search.

**Keywords.** protein language model; ESM-2; plant innate immunity; proteome-wide
classification; NBS-LRR; defense candidate prioritization; *Arabidopsis thaliana*;
*Oryza sativa*; *Vitis vinifera*
