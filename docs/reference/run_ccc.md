# Run CellChat across normalization methods

Runs CellChat cell-cell communication analysis across multiple
normalized versions of the same SpatialExperiment object.

## Usage

``` r
run_ccc(
  spelist,
  celltype_col,
  database,
  spatial_factors,
  interaction_range,
  contact_range,
  assay_name = "logcounts",
  nboot = 100,
  min_cells = 10,
  trim = 0.1,
  scale_distance = 0.1,
  variable_both = FALSE,
  verbose = TRUE
)
```

## Arguments

- spelist:

  Named list of SpatialExperiment objects. Names should correspond to
  normalization methods.

- celltype_col:

  Character. Column in colData(spe) containing cell-type labels.

- database:

  CellChat database object, e.g. CellChatDB.mouse.

- spatial_factors:

  Data frame containing CellChat spatial factors, typically ratio and
  tol.

- interaction_range:

  Numeric. Maximum spatial interaction distance.

- contact_range:

  Numeric. Contact-dependent interaction distance.

- assay_name:

  Character. Assay containing normalized expression. Default is
  "logcounts".

- nboot:

  Integer. Number of bootstrap iterations.

- min_cells:

  Integer. Minimum number of spots/cells per group.

- trim:

  Numeric. Trimming proportion for truncated mean.

- scale_distance:

  Numeric. Distance scaling used by CellChat.

- variable_both:

  Logical. Passed to identifyOverExpressedInteractions().

- verbose:

  Logical. Print progress messages.

## Value

Named list containing CellChat results for each method.
