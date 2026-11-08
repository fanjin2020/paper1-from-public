# Manuscript figures

Recreate ten main figures and fifteen supplementary figures from the included
compact source-table extracts. No study downloads, model fitting or simulation
are required. All paths are relative to this directory.

## Run

For a server with existing study environments, run from the repository root:

```bash
python paper1_plot_figures.py
```

This launcher checks installed environments without downloading packages, reports
the selected font, and writes to a new `p1_fig_YYYYMMDD_HHMMSS_*` directory beside
the checkout. Use `--python /existing/python` to select an interpreter explicitly,
or `--check-env` to check without rendering. It never refits models or downloads
study data.

Use Python 3.10 or later with the packages in `requirements.txt`:

```bash
python -m pip install -r requirements.txt
python make_figures.py
```

Existing environments with the listed dependencies can run the second command
directly. Install dependencies in a dedicated environment when needed.

The default output is `output/` inside this package, regardless of the current
working directory. To choose a different folder or subset:

```bash
python make_figures.py --output ./output --figures 1 5 S4
python make_figures.py --tiff
```

From the repository root, use `python publication_figures/make_figures.py`.
The analysis launcher does not call this script or read its data. To regenerate
analyses from public inputs, use the repository's `reproduce_from_public.py` entry.
See [FIGURE_INDEX.md](FIGURE_INDEX.md) for the publication order and source mapping.

## Outputs

- `Fig1.pdf` to `Fig10.pdf` and `FigS1.pdf` to `FigS15.pdf`: publication figures.
- `editable/*.svg`: vector figures with editable text.
- `previews/*.png`: 300-dpi previews.
- `tiff/*.tiff`: optional 600-dpi losslessly compressed raster exports.
- `layout_checks.json`: page-boundary checks for text labels.

Figures are 174 mm wide. Arial is used when installed; DejaVu Sans is the fallback.
PDF fonts are embedded as TrueType, and figure numbering matches the manuscript.
The model palette and marker shapes are defined at the top of `make_figures.py`.

## Data and statistical definitions

`source_index.csv` maps each data file to its supplementary table and figure,
records its row count, and supplies a SHA256 checksum. All input checksums are
verified before plotting. Gene-level boxplots are reconstructed from exact
summary statistics; density plots retain the gene-level coordinate pairs.

The `file` column identifies an included file relative to this package. The
`source_file` column identifies its source export for provenance; that original
export is not needed to run these figures. All required plotting inputs are
included under `data/`.

Intervals in simulation figures retain their original definitions: Wilson
intervals for binomial rates, dataset-level Monte Carlo intervals for the count
pipeline, and 1.96 Monte Carlo standard errors for coefficient bias. Gene IQRs
and variation across fixed designs are descriptive summaries. Unestimable
coefficients remain missing. The original numerical values are not changed.

Figures 9 and S13 use the completed design-geometry benchmark. For this study,
global-null FDR intervals are exact binomial Clopper-Pearson intervals;
nonnull FDR intervals are bounded Student-t approximations for the mean
dataset-level false discovery proportion, with conservative bounds for all-zero
FDP samples. Figure 9 uses two-sided 95%
intervals. The conditional screening rule uses a separate one-sided 95%
Monte Carlo upper limit of 7.5%, with nominal BH level 5%. These are distinct
quantities. The held-out rule makes 16 selections among 32 combinations,
with eight structural exclusions and eight unsupported selections.

The critical human q-bin counts and headline count-pipeline FDR values are
checked against their source data by `verify_inputs()`.

Figures 10, S14 and S15 use the completed paired count and tail-control studies.
Figure 10 displays all six targeted scenarios per arm; its scenario ranges are
not confidence intervals. Figure S15 supplies pointwise Wilson Monte Carlo
intervals from 500 datasets per scenario. Planned and retained BH families are
distinct; Gaussian controls use independent responses on the same design and
retained gene indices. The report-only entry `paper1_calibration_report.py`
exports the corresponding summaries without modifying tests or selecting an
estimator. Ordinary interval coverage is not treated as evidence of small-p
tail calibration.

## Visual references

The presentation uses a shared model palette, redundant marker shapes,
light gridlines, direct labels where useful, and panel layouts matched to the
statistical comparison. All graphs are drawn from the included study data.
No published figure or artwork is reproduced.

Open-access articles in *Interdisciplinary Sciences: Computational Life Sciences*
consulted for figure presentation:

- [Causal Transformer for Learning Embeddings from Structured Medical History Records and Multi-Source Data Integration for Complex Disease Risk Prediction](https://link.springer.com/article/10.1007/s12539-025-00749-9)
- [A Framework Integrating Spiking Cortical Circuit Modeling and Simulation-Based Inference to Probe Biomarkers of Cortical Dysfunction in Alzheimer's Disease](https://link.springer.com/article/10.1007/s12539-026-00817-8)
- [ResNeXt-Based Rescoring Model for Proteoform Characterization in Top-Down Mass Spectra](https://link.springer.com/article/10.1007/s12539-025-00701-x)

The article's [submission guidelines](https://link.springer.com/journal/12539/submission-guidelines)
remain the authority for upload requirements. Supply the vector PDF files;
the PNG files are previews, not high-resolution line-art substitutes.
