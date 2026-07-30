# Final search-panel comparison (E4 + E5 closed)

Species: *arabidopsis_thaliana*  
Curated GO defense set: 2,318 proteins  
Matched candidate count: n = 2,807

## Per-method performance

| method | candidates | recovered_TP | precision | recall | F1 |
| --- | --- | --- | --- | --- | --- |
| ESM-2 (moderate tier) | 482 | 482 | 1.0 | 0.208 | 0.344 |
| blastp_default | 2807 | 705 | 0.251 | 0.304 | 0.275 |
| blastp_filteroff | 2807 | 704 | 0.251 | 0.304 | 0.275 |
| psiblast_default | 2807 | 687 | 0.245 | 0.296 | 0.268 |
| psiblast_filteroff | 2807 | 695 | 0.248 | 0.3 | 0.272 |
| jackhmmer | 2807 | 0 | 0.0 | 0.0 | nan |
| alignment_SW | 674 | 674 | 1.0 | 0.291 | 0.451 |
| interpro_domain | 333 | 333 | 1.0 | 0.144 | 0.252 |

## ESM-2-only residual

- Published (vs jackhmmer + domain only): 72
- Against the full panel (blastp_default, blastp_filteroff, psiblast_default, psiblast_filteroff, jackhmmer, alignment_SW, interpro_domain): **65**
- Of the published set, 7 are now recovered by a search method and 65 remain ESM-2-only

## Caveat to state in Methods

Ranking by maximum bitscore across 33 independent PSI-BLAST PSSMs mixes score scales, since each query profile diverges over the three iterations. This affects the matched-count cut but not the AUPRC over the full ranking, and all search arms exceed the embeddings on both measures.
