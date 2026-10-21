# Plot CCC signaling pathway strength

Plot CCC signaling pathway strength

## Usage

``` r
plot_ccc_pathways(
  ccc_result,
  pathways = NULL,
  top_n = 10,
  method_order = NULL,
  method_cols = NULL,
  log_scale = FALSE
)
```

## Arguments

- ccc_result:

  Output from run_ccc().

- pathways:

  Optional character vector of pathways.

- top_n:

  Number of strongest pathways when pathways = NULL.

- method_order:

  Optional method order.

- method_cols:

  Optional named vector of method colors.

- log_scale:

  Logical. Use log10 y-axis.

## Value

A ggplot object.
