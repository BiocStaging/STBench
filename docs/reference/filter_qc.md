# Filter a SpatialExperiment using QC thresholds

Filter spatial locations using log library size and genes using log mean
expression.

## Usage

``` r
filter_qc(
  spe,
  library_cutoff = NULL,
  gene_cutoff = NULL,
  recalculate = TRUE,
  assay_name = "counts",
  mito_pattern = "^MT-"
)
```

## Arguments

- spe:

  A SpatialExperiment object containing QC metrics.

- library_cutoff:

  Optional minimum log library size.

- gene_cutoff:

  Optional minimum log mean gene expression.

- recalculate:

  Logical. If TRUE, QC metrics are recalculated after filtering. Default
  is TRUE.

- assay_name:

  Name of assay containing raw counts. Default is \`"counts"\`.

- mito_pattern:

  Regular expression identifying mitochondrial genes. Default is
  \`"^MT-"\`.

## Value

A filtered SpatialExperiment object.

## Details

QC metrics must first be calculated using \`add_qcs()\`.
