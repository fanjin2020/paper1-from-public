# Publication package

The manuscript contains ten main figures, fifteen supplementary figures and
Supplementary Tables S1-S23. The core analysis and the design-geometry benchmark
have separate experiment identities and output directories.

## Fast verification and figure reproduction

Run commands from the repository root. Replace the example input and output
paths under `$HOME/adjscope` with your own locations. The publication snapshot
is identified by `PUBLICATION_VERSION.json`; engine protocol versions remain
separate from the publication bundle version.

```bash
conda activate paper1
python verify_publication.py
python paper1_calibration_report.py
python paper1_plot_figures.py
```

These commands use supplied compact numerical extracts. They do not download
study data or refit models. The figure launcher prints the new output directory.

## Completed design benchmark

To summarize the completed server experiment without fitting again:

```bash
python paper1_design_fix.py \
  --previous-results "$HOME/adjscope/p1_design_completed"
```

This writes a separate, checksum-verified compact report. The original result
directory remains read-only. Exact binomial intervals are used under the global
null; the nonnull false discovery proportion uses an explicitly labeled Student-t
approximation and a conservative bound when all FDP values are zero.
Historical per-gene diagnostics that were not recorded are not
reconstructed from summaries.

## New design experiment

Only run this section when a new full simulation is required. It can take days.
The prepared inputs must come from the public-input reproduction documented in
the root README.

```bash
python paper1_design_fix.py --run-new \
  --output "$HOME/adjscope/p1_design_publication_36c" \
  --profile full --workers 36 \
  --data-work "$HOME/adjscope/paper1_complete_from_public_v1/core"
```

The named output must be new or empty. Each count dataset is checkpointed
separately. After an interruption, resume this exact experiment with:

```bash
python paper1_design_fix.py --resume --workers 36 \
  --output "$HOME/adjscope/p1_design_publication_36c"
```

Source, settings and software are checked before reusing completed work. Version
0.1.0 runs require their original source to resume; the report-only command above
can read their completed exports. Do not change a protocol version to bypass
these checks.

## Interpretation

The within-group information ratio uses exposure, group membership and nuisance
covariates. It measures estimability and the homoscedastic precision cost. The
conditional screening rule also requires the generating cage-mean amplitude,
cage standard deviation and nonnull fraction. It abstains when these are unknown
or unsupported. A 7.5% Monte Carlo screening limit is not nominal 5% FDR control.

The count benchmark uses 50 datasets per design-scenario combination. Under the
global null, one event out of 50 has a one-sided 95% upper bound of 9.14%.
Small event counts therefore produce wide uncertainty intervals.

## Files and provenance

- `publication_figures/`: all plotting code and checksum-verified inputs.
- `publication_tables/`: compact numerical sections and deterministic CSV export.
- `design_extension/`: design diagnostics, simulations, reporting and tests.
- `publication_calibration/`: completed paired-calibration inputs and read-only reports.
- `paired_count_calibration/` and `paired_tail_diagnostics/`: original simulation engines.
- `config.json`: unchanged core reproduction protocol version 2.1.0.
- `CITATION.cff`: citation metadata for the complete publication source snapshot.
- `LICENSE`: MIT license for source code and accompanying documentation.

The full pipeline generates complete gene- and dataset-level exports. The
submission tables retain summaries and explicitly defined discovery-associated
extracts; no multiple-testing correction is recomputed on a selected extract.

Public repository uploads and DOI deposits are separate author actions. The
software does not publish files or modify an existing archive.
