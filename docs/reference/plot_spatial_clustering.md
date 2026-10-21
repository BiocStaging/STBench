# Plot spatial clustering results

Visualise predicted spatial clusters and ground-truth annotations across
normalisation methods.

## Usage

``` r
plot_spatial_clustering(
  norm_list,
  cluster_col = "banksy_cluster",
  methods = names(norm_list),
  method_order = NULL,
  point_size = 0.5,
  palette = NULL,
  ncol = 3,
  legend_point_size = 4,
  theme_name = "classic",
  ground_truth_col = NULL,
  show = c("predicted", "both", "truth"),
  truth_palette = NULL
)
```

## Arguments

- norm_list:

  Named list of SpatialExperiment objects or nested lists of
  SpatialExperiment objects.

- cluster_col:

  Column containing predicted clusters.

- methods:

  Optional normalisation methods to plot.

- method_order:

  Optional method order.

- point_size:

  Spatial point size.

- palette:

  Optional named colours for predicted clusters.

- ncol:

  Number of columns for spatial facets.

- legend_point_size:

  Point size in legends.

- theme_name:

  Plot theme.

- ground_truth_col:

  Ground-truth annotation column.

- show:

  One of "both", "predicted", or "truth".

- truth_palette:

  Optional named colours for ground truth.

## Value

List containing data, palettes, individual plots, and combined plot.

## Details

When show = "both", the ground truth is displayed once per spatial
slice, followed by predicted cluster maps across all selected
normalisation methods.

Predicted clusters share one colour scale and legend. Ground truth uses
a separate palette and legend.
