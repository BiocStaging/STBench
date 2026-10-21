# Plot CellChat circle networks

Produces CellChat circle plots for one or more pathways across
normalization methods.

## Usage

``` r
plot_ccc_circle(ccc_result, pathways, methods = NULL, method_order = NULL)
```

## Arguments

- ccc_result:

  Output from run_ccc().

- pathways:

  Character vector of pathways.

- methods:

  Optional normalization methods.

- method_order:

  Optional method order.

## Value

Invisibly returns NULL.

## Details

Rows correspond to pathways and columns correspond to normalization
methods.
