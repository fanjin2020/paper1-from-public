# Results and Completion

## Output Tree

| Directory | Content |
| --- | --- |
| core/raw | Verified public downloads |
| core/prepared | Newly generated discovery, human and mouse inputs |
| core/results | Main coefficients, paired changes, families, precision, variance and calibration |
| matched/results | Eight human variants, seven mouse arms, count-pipeline summaries and failure inventory |
| benchmark/final_results | General paired calibration and five-model cage comparison |
| supplementary | Nine supplementary groups, per-stage inputs/results and completion receipts |
| figures | Eight main numerical displays, PDF/PNG and CSV |
| reporting/S01 ... S20 | Full-precision supplementary CSV exports |
| reporting/block_sources.csv | 105 workbook-block references to computed source tables |
| reporting/source_index.csv | Export size, row count and SHA256 |
| reporting/figures | Twelve supplementary PDF/PNG displays and their numerical sources |

These are source tables, not an Excel/Word typesetting engine. A single table may supply multiple display blocks. Full families are exported, including null results and missing test slots; a selected figure need not display every source row. Numerical identity checks use the estimates from the same experiment.

## Key Results

- core/results/human_results.csv.gz: coefficients, SE, df, p and reserved-family q.
- core/results/human_paired_changes.csv.gz: covariance-aware and HC3 difference tests.
- core/results/human_precision.csv.gz: exact per-gene two-factor SE decomposition.
- core/results/mouse_precision.csv.gz: design, posterior-SD and reweighting factors.
- core/results/mouse_variance.csv.gz: cage variance, boundaries and availability.
- matched/results/count_pipeline_summary.csv: scenario/model FDR, power and coverage.
- matched/results/count_pipeline_unavailable_genes.csv: full missing gene-level records.
- matched/results/count_pipeline_diagnostic_inventory.csv: dataset/model denominator counts.

Gene rows and dataset-model counts have different units. Match their explicit keys; do not infer the number of missing evaluations from the number of unique genes.

## Completion States

INTEGRATED_ANALYSES_COMPLETE means all requested main and supplementary processes and exports finished. RESULT_EXPORTS_COMPLETE identifies the final reporting export. Neither status asserts that the results equal every manuscript number. PRIMARY_ANALYSES_COMPLETE is deliberately narrower and is used only with --primary-only.

Implementation coverage in inputs/paper_coverage.json is separate from execution. Source integration, unit tests, successful server execution and numerical agreement with a manuscript are four different checks.

Unavailable individual tests and prespecified non-estimability remain explicit. Technical errors stop a stage. Whole failed empirical arms cannot count as successful matched-model tasks. Missing p values occupy reserved multiplicity slots but remain missing in exported inference.

## Figures

Main figures show human estimates, testing families, seven mouse arms, precision decomposition, general paired calibration, count validation, cage means/correlation and within-cage allocation.

Supplementary Figures S1-S12 show donor quality, astrocyte estimates, mouse panel consistency, cage ICC, layout rejection, coefficient/precision changes, design approximation, joint/paired tests, the diagnostic workflow, selected coefficient changes, donor-conditioned calibration and general-grid cancellation/power. PDF is the vector output; PNG is a preview.

## Supplementary Simulation Inputs

The completed `supplementary/simulations/attempt_*` directory contains
`simulation_animal_order.csv` and `simulation_input_order.json`, recording the
animal-to-random-draw mapping and its hashes. These accompany `scenarios.csv`,
`layout_units.csv`, `layout_registry.csv`, R versions and per-repetition results.
Existing completed experiments without the two additional order files remain
interpretable using their prepared mouse metadata and layout inputs; collection
does not invent new execution records for them.

## Collection

collect_results.py copies the completed reporting/ and figures/ directories with completion and environment records. Supplementary numerical summaries and compact per-repetition result tables are retained where required for verification. It excludes raw counts, prepared expression matrices and cell-level annotations. It creates a normal directory, never a ZIP.

Do not publish runtime results or software environments as part of the minimal source checkout. No remote archive is created automatically.
