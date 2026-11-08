# Publication Figure Index

This order accompanies "AdjScope-Tx: a diagnostic framework for covariate
adjustment in transcriptomic analysis". It differs from the
numerical display order produced by the full analysis runner.

## Main Figures

| Figure | Comparison | Supplementary source |
| --- | --- | --- |
| 1 | Count-to-inference FDR, power and interval coverage | S20.10 |
| 2 | Cage means and residual correlation | S17.1 |
| 3 | Within-cage genotype overlap | S17.1 |
| 4 | General paired coefficient-test calibration | S20.5 |
| 5 | Human coefficient estimates and discovery thresholds | S3.5; clinical_pairs.csv.gz |
| 6 | Multiple-testing families | S3.6 |
| 7 | Seven mouse analysis arms | S12.7-S12.8 and gene-level seven-arm SE summaries |
| 8 | Precision components and design information | S14 |
| 9 | Design information and conditional screening | S21.1-S21.2, S21.4, S21.7 |
| 10 | Paired-test coverage, FDR, Bonferroni and extreme tails | S23.1-S23.2 |

## Supplementary Figures

| Figure | Comparison | Supplementary source |
| --- | --- | --- |
| S1 | Human donor quality | S5 |
| S2 | Astrocyte coefficient estimates | S8 |
| S3 | Mouse panel consistency and focal genes | S12-S13 |
| S4 | Cage ICC and precision | S15 |
| S5 | Cage allocation, rejection and power | S17.1 |
| S6 | Coefficient and precision changes across human contexts | S4 |
| S7 | Design-only precision approximation | S16 |
| S8 | Joint F and paired tests | S20.1-S20.4 |
| S9 | Diagnostic workflow | Conceptual diagram; no numerical input |
| S10 | Selected human coefficient/SE components | S4.4 |
| S11 | Donor-conditioned calibration | S20.1-S20.4 |
| S12 | Cancellation and paired-test power | S20.5 |
| S13 | Analytical and simulated cage-mean bias | S21.1, S21.3, S21.8 |
| S14 | Full count paired-null coverage and both BH families | S22.1 |
| S15 | Six-arm tail controls with Monte Carlo intervals | S23.1 |

`source_index.csv` identifies every included data file, row count and checksum.
The source-file names in that index document provenance; the compact files under
`data/` contain all numerical inputs needed for these figures.

For Fig. 7, `mouse_se_ratios.csv.gz` contains 36 specification-by-arm median SE
ratios derived from the gene-level seven-arm table. It is included alongside the
42-row arm summary and 30-row matched comparison summary.

Rendering these fixed inputs recreates the publication figures. Numerical
reproduction from public expression data is a separate analysis entry, documented
in the repository README.
