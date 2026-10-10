# Statistical Methods and Definitions

[English run guide](../README_EN.md) | [Chinese version](METHODS.md)

## Analysis Questions

The human analysis compares clinical/pathological coefficients before and after adding quality covariates, using the same donors, response, and voom weights. The mouse analysis compares no-cage, cage-fixed, common-correlation, and dream models and decomposes their precision differences.

Historical discovery counts are comparison values, not fitting targets. Differences from a manuscript version should be traced to inputs and methods rather than removed by changing thresholds.

## Human Contrasts

| Contrast identifier | Target | Base covariates | Added terms |
|---|---|---|---|
| `clinical_status` | Dementia versus No dementia | Age, sex | PMI, RIN, pH |
| `clinical_status_adjusted_for_braak` | The same clinical comparison | Age, sex, Braak | PMI, RIN, pH |
| `braak_same_clinical_population` | Change per Braak stage | Age, sex | PMI, RIN, pH |

Within each subclass, all comparisons use a common complete-case population. The clinical base design supplies edgeR expression filtering, followed by TMM normalization and voom mean-variance modelling. The three contrasts and both model arms share the resulting response and weights.

Constant age or sex columns are omitted to avoid redundant columns; other rank deficiencies are not handled by silently dropping covariates. At least 100 background genes are required for residual-variance estimation. Filtered-out panel genes remain unavailable planned positions.

The primary limma empirical Bayes analysis uses `robust=FALSE` and `trend=FALSE`. Robust moderation is evaluated separately in the labelled SE implementation diagnostic.

## Paired Coefficient Differences

For common diagonal weights W and nested designs X0 and X1, define the coefficient maps and the difference between their target rows:

```text
Ak = (Xk' W Xk)^(-1) Xk' W
d  = A1[j, ] - A0[j, ]
delta = d y
```

The model-based SE is:

```text
SE(delta) = sqrt(s1^2 * sum_i(d_i^2 / w_i))
```

Here `s1^2` is the full-model unmoderated weighted residual mean square. The estimates share the same response, so adding two independent-estimate variances is inappropriate.

The HC3 alternative uses full-model residuals `e_i` and weighted leverage `h_i`:

```text
SE_HC3(delta) = sqrt(sum_i(d_i^2 * e_i^2 / (1 - h_i)^2))
```

Both tests use the full-model residual t degrees of freedom, approximately for HC3. Degenerate contrasts and leverage equal to one do not supply an ordinary test. Inference is conditional on the fitted weights and does not propagate weight estimation or panel selection. These are difference tests, not equivalence tests.

## Mouse Models

Whole-section and cortical analyses share the background passing `filterByExpr` in both scopes. The six specifications combine two scopes with linear age, categorical age, or categorical age plus slide. Sex is included; genotype is coded WT = 0 and 5X = 1.

| Model | Implementation |
|---|---|
| `baseline` | Voom/limma without a cage term |
| `cage_fixed` | Cage fixed effects, with voom weights re-estimated for that design |
| `cage_correlation` | Initial voom and common-correlation estimation, reweighted voom using that correlation, a second correlation estimate, then block fitting |
| `dream_shared` | Official dream random cage intercept using baseline response and voom weights |
| `dream_native` | Official dream after `voomWithDreamWeights` |

An additional cage-fixed fit retaining baseline weights is used for precision decomposition. It does not represent an additional biological replicate.

## Three-Factor SE Decomposition

Let B denote the baseline model, S the cage-fixed model with baseline weights, and F the cage-fixed model with re-estimated weights.

```text
Design factor     = unscaled_SE(S) / unscaled_SE(B)
Posterior factor  = sqrt(posterior_variance(S) / posterior_variance(B))
Reweighting factor = moderated_SE(F) / moderated_SE(S)

moderated_SE(F) / moderated_SE(B)
    = Design factor * Posterior factor * Reweighting factor
```

The identity is checked per gene. Separately summarized medians need not obey the same product. The unweighted design-only prediction is supplied as an approximation, not as an independent validation of the fitted model.

### Variance and boundary diagnostics

`fitVarPartModel` extracts cage variance under shared and native weights. Conditional ICC is cage variance divided by cage variance plus mean observation-specific residual variance. Mean within-cage pair correlation is also reported and has a distinct definition.

Near-zero variance, singularity, and nonconvergence are separate outcomes. A converged zero-variance estimate can remain usable. Each gene retains both planned weight-strategy positions.

Warnings in the original dream fit are localized with the same response, weights, formula, and official API. `bobyqa` and `nloptwrap` must agree under the recorded numerical tolerances and have no unresolved warnings. Accepted raw-fit rows are merged before common-background moderation. Unresolved or unlocalized warnings remain reported. A separate successful variance-extraction fit does not establish convergence of the original dream fit. See [DREAM_REPAIR_EN.md](DREAM_REPAIR_EN.md).

## SE Implementation Checks

For gene-specific weights W, the exact target unscaled SE is:

```text
sqrt([(X' W X)^(-1)]jj)
```

The implementation uses weighted QR rather than an explicit matrix inverse and checks the result against limma's direct genotype-coefficient `stdev.unscaled`.

Using the same response and weights, a group-means parameterization with `contrasts.fit` extracts 5X minus WT. Direct versus transformed contrasts and ordinary versus robust moderation form four diagnostic combinations. Equal coefficients with different SEs may reflect covariance approximations or moderation settings rather than different expression effects.

The limma documentation describes approximate unscaled contrast SEs for nonorthogonal designs with gene-specific precision weights. The primary analysis therefore uses the direct coefficient verified by weighted QR. Additional checks compute `sqrt(c' (X' W X)^(-1) c)` and refit after official `contrastAsCoef`; both must agree with the exact direct-coefficient result. Centred and uncentred age parameterizations are retained for diagnosing the approximate calculation. Exact estimates should be invariant to an age shift.

Implementation references: [limma manual](https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf) and [dream documentation](https://bioconductor.org/packages/release/bioc/vignettes/variancePartition/inst/doc/dream.html).

## Multiple-Testing Families

| Output | Planned family |
|---|---|
| Discovery panel | All testable genes within each discovery region |
| Human associations | All subclasses x three contrasts x panel genes, separately for each model arm |
| Paired differences | All subclasses x three contrasts x panel genes; model-based and HC3 tests corrected separately |
| `q_base_fixed_joint` | Two methods x six specifications x panel positions |
| `q_three_model_joint` | Baseline, cage-fixed, and common-correlation methods x six specifications x panel positions |
| `q_dream_joint` | Two dream weight strategies x six specifications x panel positions |
| `q_within_model` | Six specifications x panel positions, separately for each mouse method |

Unavailable p values temporarily contribute one to the planned BH denominator; exported p/q values remain missing. They are not counted as negative findings.

`correction_scope.csv` evaluates partitioned families and planned versus available-only denominators. `human_correction_scope_summary.csv` reports clinical discovery counts and retention of the original clinical set. These are descriptive comparisons of declared families. `human_gene_membership.csv` records symbol matching and expression filtering; symbols are case-normalized without deleting genes to recover historical test counts.

## Cage-Layout Simulations

The default grid uses the observed 74-animal linear-age layout plus three synthetic layouts with 72 animals, six animals per cage, and mixed-cage fractions 0, 0.5, and 1. The observed design retains age and sex; synthetic designs contain the intercept, genotype, and applicable cage terms. SVD removes redundant nuisance columns without deleting genotype to manufacture estimability.

The grid crosses cage correlation 0 or 0.2, genotype-related cage mean shift 0 or 0.5, and target effect 0 or 0.3, producing 32 scenarios. Each replicate generates 64 Gaussian genes. Nontarget effects alternate between 0 and 0.3; the first gene uses the scenario's target effect. Generating age and sex effects are zero, although the observed-layout fits retain those covariates.

All 64 genes contribute to common-correlation estimation; four-model performance is evaluated for the first gene. Pure-genotype cages produce `not_estimable` fixed-effect contrasts, separately from numerical `failed` fits.

There are 500 replicates per scenario and 64,000 planned model records. The simulation uses unmoderated inference. Scenario seeds and stage records make the generated data independent of task order and resumption. Observed animals are ordered by SAMPLE strings, source cage labels are sorted before integer encoding, and synthetic cages use integer ordering. Exported unit lists and matrices preserve these correspondences.

## Paired-Test Calibration

Calibration uses the current complete-case donor designs, two fixed weight profiles, three error laws, and three mean conditions. Model-based and HC3 tests are evaluated separately. Mean conditions are zero quality effects, nonzero quality effects with zero target difference, and a target difference of two analytical sampling SDs under the matched reference model. The default is 2,500 replicates.

This evaluates conditional pointwise inference under the specified designs and weights. It does not calibrate the complete limma/voom count pipeline or BH-FDR procedure.

## Performance Denominators

Outputs report bias, empirical sampling SD, mean and RMS estimated SE, coverage, interval width, rejection, and failure rates. Conditional coverage uses available intervals; a separate column reports covered-and-available outcomes over all planned repeats. Failures remain in the planned denominator and are not interpreted as absence of false positives.
