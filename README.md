# AdjScope-Tx

Public-input reproduction for "AdjScope-Tx: a diagnostic framework for covariate adjustment in transcriptomic analysis".

## License

Source code and accompanying documentation are provided under the [MIT License](LICENSE).
Public source datasets retain the terms specified by their original repositories;
this software license does not replace those data-use terms.

## Verify the submission materials

```bash
python verify_publication.py
python paper1_calibration_report.py
python paper1_plot_figures.py
```

These commands verify supplied numerical inputs, produce calibration summaries,
and redraw the publication figures without downloading data or refitting models.
See [publication instructions](docs/PUBLICATION.md) and
[publication snapshot metadata](PUBLICATION_VERSION.json).

## Count-to-Paired-Contrast Calibration

The separate [tail and weight-control extension](paired_tail_diagnostics/README.md)
evaluates six scenarios selected after inspection of the broader count grid,
using fresh seeds, independent-weight and fixed-weight controls, a Gaussian
control, extreme-p counts and Bonferroni family-wise error rates.
Run `bash run_paired_tail_diagnostics.sh --profile full --workers 36` for
3,000 new diagnostic datasets. It preserves the previous calibration code and
results; it does not select a replacement test automatically.

The new [count calibration extension](paired_count_calibration/README.md) tests
the complete count -> filtering/TMM -> estimated voom weights -> two model fits
-> paired coefficient-difference pipeline. It uses the existing paired kernel,
not the genotype-effect tests from the design benchmark. Run
`bash run_paired_count_calibration.sh --profile full --workers 36` in an existing
R environment. It creates a separate, receipt-verified output and a compact
`review.zip`. The manuscript includes this 24,000-dataset calibration in
Table S22 and the 3,000 targeted diagnostic replicates in Table S23.

## Design-Geometry Benchmark

The separate [AdjScope-Tx design extension](design_extension/README.md) adds
outcome-free within-cage information diagnostics, conditional analytic bias
checks, a five-model count benchmark and held-out design evaluation. Run
`bash run_design_extension.sh` for a 36-worker reduced pilot. It uses prepared
public count inputs from a completed reproduction, recalibrates and refits the
extension, and writes to a new directory. It never imports fitted results.
The manuscript includes 48,000 Gaussian and 4,800 count datasets from this
extension, alongside the observed-layout benchmark of 1,600 count datasets.
The extension writes separately from the core analysis.
Design information alone does not determine FDR; the extension evaluates
explicit generating assumptions and reports unsupported selections.

The included extension is version 0.2.0. Completed legacy exports can be
summarized without model fitting using `paper1_design_fix.py`; only experiments
created by the same source version may be resumed. See
[publication and execution instructions](docs/PUBLICATION.md) before choosing
between plotting, report export, and a new simulation.

## Publication Figures

The standalone [publication figure package](publication_figures/README.md) redraws
ten main and fifteen supplementary figures from compact, checksum-verified
source-table extracts included in this repository:

```bash
python publication_figures/make_figures.py
```

PDF, editable SVG and PNG outputs are written to `publication_figures/output/`.
Plotting dependencies are listed separately in `publication_figures/requirements.txt`.
On the server, `python paper1_plot_figures.py` selects an existing compatible
environment and writes to a new short-named directory beside the checkout,
without installing packages or running analyses.
See the [figure index](publication_figures/FIGURE_INDEX.md) for the manuscript order
and numerical sources.

The [compact publication tables](publication_tables/README.md) contain the
numerical sections supplied as Supplementary Tables S1-S23. Run
`python verify_publication.py` to verify these extracts, figure inputs and the
design-held-out result counts without fitting models.

This entry recreates manuscript artwork from fixed data. It does not download
study inputs or refit models. The public-input analysis below does not read these
extracts: it regenerates derived data and results. Its numerical result displays
retain their own figure order, distinct from the publication figure order.

## Run

Run commands from the repository root. Paths below use `$HOME/adjscope` as an
example output parent; replace them with your own input and output locations.

Use Linux and an existing compatible Python 3.12 / R 4.5.3 environment. The launcher does not install or upgrade software.

```bash
conda activate paper1
python reproduce_from_public.py --check-env
python reproduce_from_public.py --workers 60 \
  --output "$HOME/adjscope/paper1_complete_from_public_v1"
```

A fresh output must be empty. All analyses run by default: panel discovery, primary human/mouse models, matched sensitivities, count simulations, external cohorts, cognition, cell annotation/interactions, mixtures, allocation simulations, result exports and figures.

The convenience entry is `bash run_full_reproduction.sh`. It uses 60 workers,
chooses a new fixed-name output beside the checkout, and offers existing adjacent
paper1 download directories as verified raw caches. Use `--check-env` to inspect
the environment and release configuration without running analyses, or `--resume` to continue that exact experiment.
PAPER1_OUTPUT and PAPER1_WORKERS override its defaults.

Both `--check-env` and the analysis runner validate `config.json` against the code release before creating an output directory. Synchronize source files and configuration together; do not bypass a mismatch by changing an unrelated configuration's version label.

Only public download objects may be reused across experiments. Add `--raw-cache /path/to/downloads_or_work` to reuse objects with matching URL, size and SHA256 receipts. Prepared matrices, inferred cell labels, fitted results and simulation batches are never imported. Missing public inputs are downloaded. The starting point is deposited counts and microarray expression, not FASTQ alignment or CEL processing.

## Software

Exact central Python pins are in `requirements.txt`; R checks are in `inputs/environment.json`. Supplementary requirements are in `requirements-supplementary.txt`. Cell annotation requires compatible scikit-learn, threadpoolctl and scrublet. Their installed versions are recorded and frozen for each experiment; an exact classifier build is documented by the run's software record.

Use repeated `--r-library /existing/library` to add installed R libraries. Use `--python /existing/python` to select Python. Existing optional Python packages can be exposed with repeated `--python-package-path /existing/site-packages`; the selected interpreter's own packages retain precedence. Imports and central numerical versions must pass before data processing. Software reuse is separate from analysis-result reuse.

The launcher checks the active interpreter, virtual environments in the repository and its `paper1-from-public` subdirectory, and existing source environments under the sibling `paper1_public_evidence` directory. Python symlinks retain their virtual-environment invocation path; they are not resolved to the underlying Conda executable. A `PYTHON NOT SELECTED` message does not stop checks of the remaining candidates. To select an existing environment explicitly, use `--python "$PWD/.raw-reproduction-venv/bin/python"`. Missing compatible environments cause a stop, never automatic installation or relaxed version checks.

Sixty workers apply to independent simulation tasks. Memory-intensive cell annotation and preparation are intentionally serial or limited. R/BLAS threads are restricted inside workers. No GPU is used. Complete count-to-inference validation can take days; hardware memory alone does not predict duration.

## Monitor and Resume

```bash
python run.py --status --output "$HOME/adjscope/paper1_complete_from_public_v1"
python reproduce_from_public.py --workers 60 --resume \
  --output "$HOME/adjscope/paper1_complete_from_public_v1"
```

Repeat the original cache/environment options when resuming. Only checksum-verified completed work from this same experiment is reused. Unfinished supplementary modules restart in separate attempt directories; completed simulation scenarios have their own receipts. Changing code, scientific settings or software requires a fresh output. Worker counts may change. Never remove a writer lock to force concurrent writers.

`python run.py --plan` lists the execution graph. `python run.py --coverage` lists implemented scientific modules. Coverage means code integration, not that this machine has rerun all public data or that all manuscript numbers have matched.

## Results

To regenerate publication summaries from a completed run without fitting again,
use `python paper1_reexport_results.py --run /path/to/completed/run`.
It creates a new verified `p1_review_YYYYMMDD_HHMMSS` sibling directory and leaves
the experiment unchanged. The collector uses the same short default name. Keep the
transfer parent short and avoid duplicate nested folder names on Windows.
See [re-export instructions](docs/REEXPORT.md) for environment and provenance details.

- `core/`: newly prepared data, primary fits and donor-conditioned calibration.
- `matched/`: eight human sensitivity variants, seven mouse arms and 1,600 count-to-inference datasets.
- `benchmark/`: general paired calibration and five-model cage benchmark.
- `supplementary/`: nine public-input supplementary analysis groups.
- `figures/`: eight main result displays with numerical sources.
- `reporting/S01` through `S20`: full-precision CSV source tables.
- `reporting/block_sources.csv`: source mapping for all 112 supplementary workbook blocks.
- `reporting/Supplementary_Tables.xlsx`: 20 worksheets, 112 blocks and an S01 catalog with header-row locations.
- `reporting/workbook_manifest.json`: workbook checksum, block sources and record counts.
- `reporting/figures/`: 12 supplementary numerical displays.

Repeated table/figure views share numerical sources. The workbook includes S6.6 (context-specific sex counts), S10.2 (CAMERA) and S20.13-S20.17 (dataset-level count results, paired results, unavailable positions, generator parameters and complete simulation jobs). Its data come from the completed experiment, not a reference workbook. The analysis entry does not regenerate Word layout or pixel-identical publication artwork. Use the separate publication figure entry for manuscript artwork. Scientific results can be compared by their explicit gene, model, sample and testing-family keys.

S1.2 gives public data entry points and input versions. Original download receipts remain in the run records. Protocol filenames such as `mouse_units.csv` identify generated inputs. Machine-specific directories are omitted only from the workbook display; complete protocol records remain in CSV and JSON. The complete job list is expanded in S20.17 instead of being truncated to Excel's cell-length limit.

The archived v2.1.0 export map contains 105 blocks. This source update uses reporting version `publication_exports_20261026_v3` and exports 112 blocks directly. The re-exporter accepts completed 105-block runs without refitting them. No public tag or archived DOI is changed by editing this checkout.

S18 simulations retain the order of the prepared animal metadata. Their draw-position
table and input/order hashes are recorded with the simulation outputs. Use the same
input order, seeds, settings and software when comparing individual realizations;
changing only row order can change finite Monte Carlo summaries. See [Methods](docs/METHODS.md).

Collect the completed result tables without raw data or ZIP files:

```bash
python collect_results.py \
  --run "$HOME/adjscope/paper1_complete_from_public_v1" \
  --output "$HOME/adjscope/paper1_complete_results_for_review"
```

## Code Map

| Entry | Purpose |
| --- | --- |
| reproduce_from_public.py | Existing-environment selection; no package installation |
| run.py | One experiment, execution graph and explicit continuation |
| pipeline/core.py | Public inputs, aggregation, primary inference |
| pipeline/matched.py | Matched sensitivities and complete count simulations |
| pipeline/benchmark.py | General paired and five-model cage benchmarks |
| pipeline/supplementary/driver.py | Age, quality, external cohorts, cognition, annotation, interactions, mixtures and layouts |
| pipeline/reporting.py | Derived summaries and 20-section numeric export |
| pipeline/workbook.py | Streaming 112-block XLSX assembly; no extra spreadsheet package |
| pipeline/figures.py; pipeline/supplementary_figures.py | Data-driven PDF/PNG displays |
| publication_figures/make_figures.py | Standalone publication artwork from included compact data; no fitting |
| pipeline/cache.py; pipeline/lifecycle.py | Download-only reuse and experiment identity |
| analysis/ | R preparation, fitting and inference |
| tests/ | Synthetic numerical, file-boundary and workflow tests |

```bash
python -m unittest discover -s tests -v
```

Tests and rendered synthetic fixtures are not scientific replications. The small,
curated extracts in `publication_figures/data/` are included for figure recreation.
Do not add raw data, full fitted outputs, environments or runtime logs to the source
repository. No upload or DOI creation is performed by these scripts.

[中文说明](README_CN.md) | [Data](docs/DATA.md) | [Methods](docs/METHODS.md) | [Outputs](docs/OUTPUTS.md)
