# Plot CCC ligand-receptor interactions

Plot CCC ligand-receptor interactions

## Usage

``` r
plot_ccc_lr(
  ccc_result,
  pathways = NULL,
  interactions = NULL,
  top_n = 15,
  method_order = NULL,
  method_cols = NULL,
  log_scale = FALSE
)
```

## Arguments

- ccc_result:

  Output from run_ccc().

- pathways:

  Optional pathways to include.

- interactions:

  Optional LR interactions to include.

- top_n:

  Number of strongest LR interactions.

- method_order:

  Optional method order.

- method_cols:

  Optional named vector of method colors.

- log_scale:

  Logical. Use log10 y-axis.

## Value

A ggplot object.
