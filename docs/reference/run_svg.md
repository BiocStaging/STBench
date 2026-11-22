# Run spatially variable gene detection

Run spatially variable gene detection

## Usage

``` r
run_svg(norm_list, method = "nnSVG", ...)
```

## Arguments

- norm_list:

  Named list of normalised SpatialExperiment objects.

- method:

  SVG detection method.

- ...:

  Additional parameters passed to the SVG method.

## Value

Named list containing SVG results.

## Examples

``` r
data(norm_example)

svg_results <- run_svg(
  norm_list = norm_example,
  method = "nnSVG",
  n_threads = 1,
  verbose = FALSE
)

names(svg_results)
#> [1] "None"    "SpaNorm"
head(as.data.frame(
  SummarizedExperiment::rowData(svg_results[[1]])
))
#>             sigma.sq     tau.sq       phi     loglik runtime       mean
#> HSPB7   0.0004910903 0.02677370  8.539111  361.31556   0.072 0.03561908
#> KCNC3   0.0588907271 0.17492278 27.562535 -630.84074   0.115 0.24149342
#> COL25A1 0.0075259833 0.07753537 18.737159 -170.74820   0.106 0.09589966
#> MAL     0.1803866545 0.24361150 33.665686 -889.39771   0.099 0.36704269
#> CPXM2   0.0002814358 0.05030879 54.400979   68.97371   0.074 0.05788464
#> TM4SF4  0.0001563186 0.01502439 57.799690  637.14608   0.079 0.01878617
#>                var     spcov     prop_sv  loglik_lm    LR_stat rank
#> HSPB7   0.02729862 0.6221541 0.018011884  360.65616  1.3188143   26
#> KCNC3   0.23457624 1.0048891 0.251870514 -654.58969 47.4978949    6
#> COL25A1 0.08512410 0.9046166 0.088477115 -176.13722 10.7780297   14
#> MAL     0.42515362 1.1571392 0.425442070 -935.27388 91.7523443    2
#> CPXM2   0.05062784 0.2898187 0.005563047   69.11804 -0.2886487   48
#> TM4SF4  0.01519245 0.6655291 0.010297185  637.26597 -0.2397655   45
#>                 pval         padj
#> HSPB7   5.171579e-01 9.945343e-01
#> KCNC3   4.852474e-11 4.043728e-10
#> COL25A1 4.566470e-03 1.630882e-02
#> MAL     0.000000e+00 0.000000e+00
#> CPXM2   1.000000e+00 1.000000e+00
#> TM4SF4  1.000000e+00 1.000000e+00
```
