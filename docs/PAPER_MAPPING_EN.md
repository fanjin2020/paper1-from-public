# Manuscript and Output Mapping

[English run guide](../README_EN.md) | [Chinese version](PAPER_MAPPING.md)

This mapping refers to the Paper 1 manuscript assembled on 9 October 2026. Code figure numbers identify analysis modules; manuscript figure numbers identify the publication layout.

## Core Result Sources

All paths below are relative to the work directory's `results/` folder.

| Manuscript analysis | Computed outputs |
|---|---|
| Human clinical coefficients, SEs, and threshold transitions | `human_results.csv.gz`, restricted to `clinical_status` |
| Paired coefficient differences and HC3 sensitivity | `human_paired_changes.csv.gz` |
| Testing families and denominator choices | `correction_scope.csv`, `human_correction_scope_summary.csv` |
| Five mouse modelling strategies | `mouse_results.csv.gz`, `mouse_background_fits.csv.gz` |
| Cage-fixed precision decomposition | `mouse_precision.csv.gz` |
| Cage variance, ICC, and boundary fits | `mouse_variance.csv.gz`, `mouse_dream_diagnostics.csv` |
| Four-model layout simulation | `simulation_summary.csv`, `simulation_failures.csv` |
| Conditional paired-test calibration | `calibration/calibration_summary.csv` |

See [OUTPUTS_EN.md](OUTPUTS_EN.md) for row definitions and statistical fields.

## Figure Mapping

| Generated figure identifier | Content | Main manuscript correspondence |
|---|---|---|
| `Figure1_human_effect_precision` | Clinical coefficients, q values, and paired ratios | Figure 1; some quantities also inform Figure 7 |
| `Figure2_paired_changes` | Selected clinical coefficient differences and intervals | Human results accompanying Figure 1 and paired-test supplementary tables |
| `Figure3_mouse_precision` | Six-specification SE factor distributions | Figure 4 |
| `Figure4_cage_variance` | ICC and model availability | Figure 3 and supplementary Figure S4 |
| `Figure5_design_simulation` | Four-layout performance | Figures 5 and 6 and supplementary Figure S5 |
| `Figure6_paired_test_calibration` | Calibration, SE/SD, and power | Figure 8 and the paired-calibration supplementary display |
| `Figure7_mouse_testing_families` | Mouse discovery counts under explicit families | Figure 3 and related source tables |
| `Figure8_SE_implementation` | Direct/contrast SE extraction and moderation settings | SE implementation diagnostic data |

Exact generated filenames are defined in `pipeline/figures.py` and listed in the output figure index. Figures use the current calculated results, but do not automatically recreate Word layouts or all supplementary displays. Do not overwrite submission figures merely because their numbers match.

## Checked Result Landmarks

The following values identify the scope of the checked run. They are not inputs to gene selection or fitting:

| Quantity | Checked value |
|---|---|
| Human planned family | 17,064 positions per arm |
| Human testable positions | 11,232 per arm |
| Testable clinical-status positions | 3,744 |
| Clinical discoveries under the joint family | 40 base, zero quality-adjusted |
| Median absolute coefficient ratio within those 40 | Approximately 0.928 |
| Median SE ratio within those 40 | Approximately 1.019 |
| Common mouse background | 13,312 genes |
| Six-specification precision records | 79,872 |
| Background-gene median fixed/base SE ratios | 1.495-1.745 across specifications |
| Two-weight-strategy variance-component records | 159,744 |
| Unresolved dream warnings in the checked completed run | Zero |

The 1.495-1.745 range is for the full gene background, not the 138-gene testable candidate panel. Keep gene universe, model, and testing family explicit when comparing summaries.

## Historical Anchor Table

`pipeline/audit.py` still writes `manuscript_comparison.csv` using earlier manuscript anchors. Its six mouse SE references are not the current manuscript values. The unchanged source preserves identity with the completed analysis; these reference constants are used only for an additional comparison and do not enter model fitting or multiple-testing calculations.

A failed historical-anchor comparison therefore requires reading the reference version. It does not, by itself, mean that the current manuscript is inconsistent or that the scientific run failed.

## Analyses Outside This Entry Point

Other cohorts' exploratory analyses, additional broad-cell-class mixture experiments, historical simulation grids, and the full publication table/figure formatting have separate workflows. This repository entry point is not a complete rerun of every supplementary study.

Public files and ortholog annotations can change at the source. Download receipts record URL, time, size, and SHA-256. Complete verified inputs can be reused to refit models without downloading them again.
