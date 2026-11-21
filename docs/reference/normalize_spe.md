# Apply multiple normalisation methods

Apply multiple normalisation methods

## Usage

``` r
normalize_spe(
  spe,
  methods = c("none", "scran", "sct", "spanorm", "giotto"),
  params = list()
)
```

## Arguments

- spe:

  A SpatialExperiment object.

- methods:

  Character vector containing normalisation methods.

- params:

  Named list containing method-specific parameters.

## Value

A named list of normalised objects.

## Examples

``` r
data(spe_example)

norm <- normalize_spe(
  spe_example,
  methods = c("None", "SpaNorm"),
  params = list(
    spanorm = list(
      sample_p = 0.15,
      df_tps = 6,
      adj_method = "logpac"
    )
  )
)
#> norm_none completed in 0.001 mins
#> (1/2) Fitting SpaNorm model
#> Warning: more singular values/vectors requested than available
#> Warning: You're computing too large a percentage of total singular values, use a standard svd instead.
#> Warning: did not converge--results might be invalid!; try increasing work or maxit
#> Warning: encountered non-positive size factor estimates
#> Warning: encountered non-positive size factor estimates
#> Warning: encountered non-positive size factor estimates
#> Warning: encountered non-positive size factor estimates
#> Warning: encountered non-positive size factor estimates
#> 944 cells/spots sampled to fit model
#> iter:  1, estimating gene-wise dispersion
#> iter:  1, log-likelihood: -22005.930127
#> iter:  1, fitting NB model
#> iter:  1, iter:  1, log-likelihood: -22005.930127
#> iter:  1, iter:  2, log-likelihood: -19566.704287
#> iter:  1, iter:  3, log-likelihood: -18244.866630
#> iter:  1, iter:  4, log-likelihood: -17469.350498
#> iter:  1, iter:  5, log-likelihood: -17126.983282
#> iter:  1, iter:  6, log-likelihood: -16972.268055
#> iter:  1, iter:  7, log-likelihood: -16921.251322
#> iter:  1, iter:  8, log-likelihood: -16879.850748
#> iter:  1, iter:  9, log-likelihood: -16860.633883
#> iter:  1, iter: 10, log-likelihood: -16852.821303
#> iter:  1, iter: 11, log-likelihood: -16852.733699 (converged)
#> iter:  2, estimating gene-wise dispersion
#> iter:  2, log-likelihood: -16670.454257
#> iter:  2, fitting NB model
#> iter:  2, iter:  1, log-likelihood: -16670.454257
#> iter:  2, iter:  2, log-likelihood: -16623.611272
#> iter:  2, iter:  2, log-likelihood: -16623.611272
#> iter:  2, iter:  2, log-likelihood: -16623.611272
#> iter:  2, iter:  3, log-likelihood: -16623.611272 (converged)
#> iter:  3, estimating gene-wise dispersion
#> iter:  3, log-likelihood: -16621.106109
#> iter:  3, fitting NB model
#> iter:  3, iter:  1, log-likelihood: -16621.106109
#> iter:  3, iter:  1, log-likelihood: -16621.106109
#> iter:  3, iter:  1, log-likelihood: -16621.106109
#> iter:  3, iter:  2, log-likelihood: -16621.106109
#> iter:  3, iter:  2, log-likelihood: -16621.106109
#> iter:  3, iter:  2, log-likelihood: -16621.106109
#> iter:  3, iter:  3, log-likelihood: -16621.106109 (converged)
#> iter:  4, estimating gene-wise dispersion
#> iter:  4, log-likelihood: -16621.106109
#> iter:  4, fitting NB model
#> iter:  4, iter:  1, log-likelihood: -16621.106109
#> iter:  4, iter:  1, log-likelihood: -16621.106109
#> iter:  4, iter:  1, log-likelihood: -16621.106109
#> iter:  4, iter:  2, log-likelihood: -16621.106109
#> iter:  4, iter:  2, log-likelihood: -16621.106109
#> iter:  4, iter:  2, log-likelihood: -16621.106109
#> iter:  4, iter:  3, log-likelihood: -16621.106109 (converged)
#> iter:  5, log-likelihood: -16621.106109 (converged)
#> (2/2) Normalising data
#> norm_spanorm completed in 0.083 mins
```
