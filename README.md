# Model Adjustment in Alzheimer Disease Transcriptomes

Public-input reproduction for "Distinguishing coefficient changes from precision costs in brain transcriptome analysis", prepared for Journal of Neuroscience Methods.

## Run

Use Linux and an existing compatible Python 3.12 / R 4.5.3 environment. The launcher does not install or upgrade software.

```bash
conda activate paper1
cd /home/research/git/repository/paper1-from-publicV1
python reproduce_from_public.py --check-env
python reproduce_from_public.py --workers 60 \
  --output /home/research/git/repository/paper1_complete_from_public_v1
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
python run.py --status --output /home/research/git/repository/paper1_complete_from_public_v1
python reproduce_from_public.py --workers 60 --resume \
  --output /home/research/git/repository/paper1_complete_from_public_v1
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

Repeated table/figure views share numerical sources. The workbook includes S6.6 (context-specific sex counts), S10.2 (CAMERA) and S20.13-S20.17 (dataset-level count results, paired results, unavailable positions, generator parameters and complete simulation jobs). Its data come from the completed experiment, not a reference workbook. Word layout and pixel-identical artwork are not regenerated. Scientific results can be compared by their explicit gene, model, sample and testing-family keys.

S1.2 gives public data entry points and input versions. Original download receipts remain in the run records. Protocol filenames such as `mouse_units.csv` identify generated inputs. Machine-specific directories are omitted only from the workbook display; complete protocol records remain in CSV and JSON. The complete job list is expanded in S20.17 instead of being truncated to Excel's cell-length limit.

The archived v2.1.0 export map contains 105 blocks. This source update uses reporting version `publication_exports_20261026_v3` and exports 112 blocks directly. The re-exporter accepts completed 105-block runs without refitting them. No public tag or archived DOI is changed by editing this checkout.

S18 simulations retain the order of the prepared animal metadata. Their draw-position
table and input/order hashes are recorded with the simulation outputs. Use the same
input order, seeds, settings and software when comparing individual realizations;
changing only row order can change finite Monte Carlo summaries. See [Methods](docs/METHODS.md).

Collect the completed result tables without raw data or ZIP files:

```bash
python collect_results.py \
  --run /home/research/git/repository/paper1_complete_from_public_v1 \
  --output /home/research/git/repository/paper1_complete_results_for_review
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
| pipeline/cache.py; pipeline/lifecycle.py | Download-only reuse and experiment identity |
| analysis/ | R preparation, fitting and inference |
| tests/ | Synthetic numerical, file-boundary and workflow tests |

```bash
python -m unittest discover -s tests -v
```

Tests and rendered synthetic fixtures are not scientific replications. Do not upload raw data, fitted outputs, environments or runtime logs into the source repository. No upload or DOI creation is performed by these scripts.

[中文说明](README_CN.md) | [Data](docs/DATA.md) | [Methods](docs/METHODS.md) | [Outputs](docs/OUTPUTS.md)
