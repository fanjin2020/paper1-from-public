# Release Validation

[English run guide](../README_EN.md) | [Chinese version](VALIDATION.md)

## Code Identity

Analysis version: **1.2.2**. Original GitHub package date: **9 October 2026**. English documentation added on **10 October 2026**.

The package's 32 execution, configuration, and test files were checked byte-for-byte against the executed copies saved with the completed server analysis. Their SHA-256 values matched. Documentation and packaging changes did not alter fitting, filtering, correction, random seeds, or optimizer rules.

`release_manifest.json` records packaged file checksums. It describes file identity, not a guarantee that a future run on different inputs or software will produce identical results. The manifest excludes its own checksum.

## Completed Server Analysis

The checked full analysis finished with `status: completed` and version 1.2.2.

| Check | Result |
|---|---:|
| Result integrity | Passed |
| Mouse expression numerical fit failures | 0 |
| Variance-component fit failures | 0 |
| Four-model simulation numerical failures | 0 |
| dream warning records | 1 |
| Unresolved dream warnings | 0 |
| Maximum relative error in the exact-SE check | 0 |
| Structurally non-estimable simulation positions | 4,000 |

The 4,000 non-estimable positions arise from the inability to separate genotype and cage fixed effects in pure-genotype cages. They are distinct from numerical failure. The recorded near-boundary dream warning was retained in the diagnostics and resolved through agreement between two optimizers.

## Packaging Checks

The original release preparation reran Python unit tests and checked ZIP integrity and file checksums. It did not download the counts or repeat the full R analysis during packaging. Full-analysis evidence comes from the completed server run; local unit tests provide a separate check.

The English-documentation update checks document links, command-option names, file presence, and unchanged original code checksums. It does not constitute another execution of the Linux analysis.

Every actual run still executes R dependency checks, numerical regression tests, and a small official dream API smoke test before real-data analysis. Matching a release version does not replace checking the local environment.

## Interpreting Validation Outputs

Use `status.json` and `results/result_audit.json` for the requested task's completion and diagnostic state. Read any failure or warning records alongside the corresponding denominators.

Historical manuscript anchors are documented in [PAPER_MAPPING_EN.md](PAPER_MAPPING_EN.md). Their comparisons are separate from numerical fit validation. Current results should be matched to the intended manuscript version, model, gene universe, and testing family.
