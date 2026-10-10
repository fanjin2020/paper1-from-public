# Outputs and Result Interpretation

[English run guide](../README_EN.md) | [Chinese version](OUTPUTS.md)

## Main Result Tables

Paths in this section are relative to the work directory's `results/` folder.

| File | One row represents | Principal information |
|---|---|---|
| `panel.csv` | One selected frontal gene | `gene_symbol`, `beta`, `p`, `q`, `n` |
| `ortholog_panel.csv` | One human panel endpoint | Human symbol and `mouse_gene` mapping |
| `human_results.csv.gz` | Gene x subclass x contrast x arm | `beta`, `SE`, `df`, `p`, `q_joint`, intervals, status |
| `human_paired_changes.csv.gz` | Quality-adjustment difference for one endpoint | `delta`, `SE_delta`, `p`, `q_joint`, `SE_HC3`, `p_HC3`, `q_HC3_joint` |
| `human_designs.csv` | Subclass x contrast | Sample count, design dimensions, estimability |
| `human_gene_membership.csv` | Subclass x panel gene | Input presence, symbol matching, filter retention |
| `correction_scope.csv` | Declared family x denominator strategy | Planned/available counts, denominator, discoveries |
| `human_correction_scope_summary.csv` | Arm x scope x denominator strategy | Clinical discoveries and retention of the original clinical set |
| `mouse_results.csv.gz` | Panel endpoint x specification x model | Effects, intervals, `q_within_model`, `q_base_fixed_joint`, `q_three_model_joint`, `q_dream_joint`, `rho` |
| `mouse_precision.csv.gz` | Background gene x specification | `SE_ratio`, `design_factor`, `posterior_factor`, `reweight_factor` |
| `mouse_variance.csv.gz` | Background gene x specification x weight strategy | Cage variance, residual scale, ICC, pair correlation, boundaries, convergence |
| `mouse_SE_implementation_audit.csv.gz` | Background gene x specification x principal model | Direct/contrast/exact-QR SE checks and ordinary/robust moderation comparisons |
| `mouse_background_fits.csv.gz` | Background gene x specification x model | Full-background effects, intervals, availability |
| `mouse_model_status.csv` | Specification x model | Planned, available, failed |
| `mouse_fit_failures.csv` | Unavailable background fit | Gene, specification, model, status |
| `mouse_model_messages.csv` | Model warning or error | Specification, model, original message |
| `mouse_dream_diagnostics.csv` | Localized gene warning or unlocalized warning | Gene, specification, model, final status, message |
| `mouse_dream_retry_attempts.csv` | Optimizer retry | Raw estimates, warnings, acceptance checks |
| `mouse_dream_background_impact.csv.gz` | Background gene in an affected model | Moderated results with and without restored raw-fit rows |
| `mouse_designs.csv` | Mouse specification | Animals, cages, mixed cages, comparison information |
| `simulation_summary.csv` | Layout x correlation x shift x target effect x model | Bias, empirical SD, estimated SE, coverage, rejection, width, availability |
| `simulation_failures.csv` | Numerical failure or structural non-estimability | Scenario, replicate, and reason |
| `calibration/calibration_summary.csv` | Paired-test calibration condition | Performance measures and quality-covariate F-test rejection |
| `result_audit.json` | Overall result audit | Denominator checks, failures, historical-anchor comparisons |
| `manuscript_comparison.csv` | Historical manuscript anchor | Old value, calculated value, tolerance, comparison result |

`simulation_replicates.csv.gz` retains replicate-level outputs locally but is excluded from the review ZIP. Downloaded objects and prepared count matrices also remain outside the compact review package.

Compressed CSVs can be read directly:

```python
import pandas as pd

human = pd.read_csv("work_public/results/human_results.csv.gz")
```

## Availability and Failures

A header-only `mouse_dream_diagnostics.csv` means that no warning requiring localization was captured. A warning row without a gene identifier means that localization remains unresolved; it does not mean that the warning can be ignored.

Unavailable fits remain in planned denominators and are not negative findings. In the current simulation output:

- `numerical_failed` counts numerical failures.
- `not_estimable` counts structurally unidentified contrasts.
- `failure_rate` is the numerical failure fraction.
- `not_estimable_rate` is the structural non-estimability fraction.
- `unavailable_rate` includes both causes of unavailability.
- `coverage` and `rejection` use available fits.
- `coverage_all_planned` and `rejection_all_planned` use all planned repeats.

Older manuscript source tables used `failure_rate` for overall unavailability. Compare its meaning with the current `unavailable_rate`, not with the current numerical `failure_rate`.

Likewise, a correlation-parameter boundary rate and a dream singular-fit fraction describe different diagnostics. An unavailable or inapplicable diagnostic should not be silently recoded as zero.

## Figures

Each analysis figure has a same-named source CSV, PDF, and PNG. The PNG output uses 600 dpi. Larger supporting tables remain compressed under `results/`.

| Code figure | Content |
|---|---|
| `Figure1_human_effect_precision` | Clinical coefficients, q transitions, coefficient ratios, and SE ratios |
| `Figure2_paired_changes` | Paired differences and 95% intervals for the base-significant clinical records |
| `Figure3_mouse_precision` | SE decomposition distributions across six mouse specifications |
| `Figure4_cage_variance` | ICC under shared/native weights, availability, failures, and boundaries |
| `Figure5_design_simulation` | Bias, rejection, width, coverage, and failure across layouts and models |
| `Figure6_paired_test_calibration` | Estimated SE versus empirical SD, null rejection, and two-SD power |
| `Figure7_mouse_testing_families` | Mouse discoveries under method-specific and joint corrections |
| `Figure8_SE_implementation` | Contrast parameterization and empirical Bayes contributions to SE differences |

These identifiers are analysis-module figure numbers. They do not automatically reproduce the manuscript figure sequence or layout. See [PAPER_MAPPING_EN.md](PAPER_MAPPING_EN.md).

## Comparing a Run with the Manuscript

Check inputs and units before significance counts:

1. Compare panel size and gene identities.
2. Check SEA-AD donor eligibility, subclasses, complete-case populations, and quality-variable definitions.
3. Check mouse inclusion/exclusion, animal IDs, cage and slide labels, and cortical annotations.
4. Compare gene backgrounds, TMM, voom, empirical Bayes settings, and testing-family definitions.
5. Match records by gene, context, contrast, and model before comparing coefficients, SEs, intervals, and q values.

`manuscript_comparison.csv` contains earlier manuscript anchors, including obsolete mouse SE references. Consult the current manuscript version and [PAPER_MAPPING_EN.md](PAPER_MAPPING_EN.md) when interpreting those comparisons. Shared non-significance is not evidence of numerical reproduction, and the BH scope must not be changed simply to recover a desired discovery count.

## Scope of This Entry Point

The workflow covers the core public-data analyses. Other exploratory cohorts, all historical supplementary simulations, and the exact layout of every manuscript figure are not reproduced by this single entry point. Read the manuscript's supplementary methods for those separate analyses.
