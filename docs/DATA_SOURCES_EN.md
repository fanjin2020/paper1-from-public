# Data Sources and Processing

[English run guide](../README_EN.md) | [Chinese version](DATA_SOURCES.md)

## SEA-AD: Human Postmortem Single-Nucleus Counts

- [Official data description](https://brain-map.org/consortia/sea-ad/our-data)
- [Source S3 bucket](https://sea-ad-single-cell-profiling.s3.amazonaws.com/)

The program lists `PFC/RNAseq/donors_objects/` and selects the 80 objects ending in `_SEAAD_DFC_RNAseq_final-nuclei.2026-06-22.h5ad`. The release date is recorded in `config.json`; a later release is not selected automatically. The inventory records object names, sizes, and ETags.

The donor workbook is `sea-ad_cohort_donor_metadata_072524.xlsx`, linked from the official data page. Its full URL is stored in `URLS` in `pipeline/sources.py`. Required fields include Donor ID, Cognitive Status, PMI, RIN, and Brain pH.

### Processing

1. Check that each H5AD donor identifier agrees with its filename, the region is DFC, and cell barcodes are unique.
2. Retain nuclei marked by the source as `Used in analysis`, excluding `Neurotypical reference`. Parse Braak stage, age, and sex from their designated fields.
3. Read raw sparse counts from `layers/UMIs`. Normalized expression or `X` cannot replace this layer.
4. Sum counts by donor and source Subclass. Sum rows with duplicate gene symbols. A gene absent from the source is not assigned a fabricated zero. Inconsistent gene spaces across donors stop preparation for review.
5. Join quality measurements by Donor ID and check clinical status against the H5AD metadata.
6. Restrict clinical comparisons to Dementia and No dementia with complete quality covariates. Require at least 20 recovered nuclei per donor-subclass unit, at least five donors per clinical group, and at least five residual degrees of freedom. Model comparisons within a subclass use the same eligible donors.

PMI, RIN, and pH have different measurement contexts. In particular, workbook RIN is not a separate quality measurement for every sequenced library. Interpretation should retain the source definitions. Cell annotations are those supplied by the source.

## GSE233208: Mouse Spatial Expression and Animal Identity

- [GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233208)
- [Source publication](https://doi.org/10.1038/s41588-024-01961-x)

The workflow downloads `GSE233208_5XFAD_seurat_processed_annotated.rds.gz` and six sample metadata files: pilot, 8samples, and four date-specific sample tables. Download URLs are enumerated in the code.

Animal identity is reconstructed from `Supp. Table 2` in `SupplementaryTables.xlsx`, supplied in the publication's supplementary archive. GEO metadata establish sample eligibility; the supplementary table supplies animal, cage, and ear-notch identifiers. Tables are joined exactly by `SAMPLE`, then checked for age, sex, and genotype agreement.

### Processing

- Use the object's `SAMPLE` field rather than incomplete `SampleID` values. Animal identifiers are not inferred from spot barcodes or filenames.
- Preserve inclusion/exclusion records. Pilot samples with missing sex are not silently added.
- Stop if an animal maps to multiple SAMPLE values so that repeated-section relationships can be reviewed.
- Require unique cage/ear-notch combinations. Spatial spots are aggregated rather than treated as independent animals.
- Read the Spatial assay counts layer, including the SeuratObject 5 `layer` interface.
- Aggregate all eligible spots for the whole-section scope. The cortical scope includes only `ctx-deep-layers` and `ctx-upper-layers`, excluding olfactory cortex.
- Require at least 100 cortical spots. Both scopes and all models use the same eligible animals.
- Calculate animal counts, cage counts, and mixed-genotype cage counts from the retrieved metadata rather than selecting animals from a historical 74-animal list.

The public RDS may have two gzip layers. Python removes the extra outer layer, and R reads the remaining compressed object. Preparation stops if the object cannot be read or the Spatial counts layer is absent; processed expression is not substituted.

## GSE118553: Candidate-Panel Definition

[GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118553)

The workflow uses the GEO series matrix and official GPL10558 platform annotation. The matrix contains author-processed microarray expression, not raw probe intensities. The pipeline does not add an unrequested log transformation.

1. Read individual, tissue, disease state, gender, and age for each GSM.
2. Retain frontal, entorhinal, and temporal cortex records. Collapse technical replicates within donor-region combinations by the median at the probe level.
3. Accept uniquely annotated, valid gene symbols; collapse multiple probes per gene by their median.
4. Fit AD versus control separately in frontal and entorhinal cortex, adjusting for centred age and sex with female as the reference. Require at least ten complete cases per group.
5. Apply BH correction over all valid genes within each region. Frontal genes with q < 0.05 define the candidate panel, independently of SEA-AD and mouse results.

The manuscript's panel contains 237 genes. The pipeline reconstructs membership rather than shipping those names as a fixed input. A different panel size requires checking platform annotation, source expression, and model definitions before claiming reproduction of that panel.

## MGI: Human-Mouse Orthologs

[Official homology report](https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt)

Using DB Class Key, retain homology classes with one human symbol and one mouse symbol and require globally unique one-to-one mappings. Missing or ambiguous mappings remain unavailable; their planned panel positions are preserved.

MGI updates its report. Download time, URL, and SHA-256 are recorded locally so annotation-version differences can be traced.

## Download Integrity and Reuse

Downloads retain TLS verification, support resumption protected by ETag or Last-Modified, and check response ranges and complete lengths. Completed downloads receive SHA-256 records. Reuse requires agreement among the receipt, URL, length, and checksum.

A locally computed checksum detects subsequent changes to that local file. It is not a publisher-provided authentication value when the publisher supplies no checksum. Conflicting source metadata trigger review rather than guessed replacements.

Use the [English run guide](../README_EN.md#2-reuse-existing-data-without-downloading) to reuse a complete work directory without acquiring the data again.
