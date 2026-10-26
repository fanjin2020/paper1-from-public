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
| reporting/block_sources.csv | 112 workbook-block references to computed source tables |
| reporting/Supplementary_Tables.xlsx | 20 worksheets and 112 complete sections, with the S01 catalog |
| reporting/workbook_manifest.json | Workbook SHA256, source mapping, counts and header rows |
| reporting/source_index.csv | Export size, row count and SHA256 |
| reporting/figures | Twelve supplementary PDF/PNG displays and their numerical sources |

The workbook is assembled from the indexed CSVs. A single source can supply multiple blocks. Full families are exported, including null results and missing test slots; a selected figure need not display every source row. Word layout is not regenerated.

S6.6 contains context-specific human sex counts, S10.2 contains CAMERA results, and S20.13-S20.17 contain per-dataset count metrics, paired metrics, unavailable positions, generator parameters and complete simulation jobs. These seven blocks extend the 105-block layout in archived v2.1.0. The new exporter requires exactly 112 block identifiers; collection can also read completed legacy 105-block runs for re-export.

S1.2 uses public entry points and input-version definitions. Workbook protocol displays omit machine-specific directories and identify generated input filenames. CSV/JSON protocol records and download receipts remain complete. Protocol job arrays are displayed as a reference to the losslessly expanded S20.17 instead of exceeding Excel's 32,767-character cell limit. Unexpected overlong cells stop the export rather than being truncated.

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
