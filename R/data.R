#' Example spatial transcriptomics dataset
#'
#' A subset of the publicly available STARmap PLUS mouse brain
#' spatial transcriptomics dataset, obtained using the
#' `STARmapPLUS_mouseBrain()` function from the Bioconductor
#' package `STexampleData`.
#'
#' The dataset was generated using STARmap PLUS, an imaging-based
#' spatial transcriptomics technology.
#'
#' A random subset of 100 genes and 200 cells was selected
#' using seed 123. Genes and cells with zero total counts
#' were removed.
#'
#' The dataset contains spatial coordinates, raw gene counts,
#' and original cell-type and tissue-region annotations.
#'
#' Two additional annotation columns were created:
#' `cell_type`, containing molecular cell-type annotations,
#' and `layer_guess`, containing molecular tissue-region
#' annotations.
#'
#' @name spe_example
#' @docType data
#' @format A SpatialExperiment object.
#' @source Bioconductor STexampleData package:
#'   \url{https://bioconductor.org/packages/STexampleData/}
NULL


#' Example normalised spatial transcriptomics datasets
#'
#' Normalisation results generated from `spe_example`
#' using the STBench `normalize_spe()` function.
#'
#' Two normalisation methods are included:
#' `None`, representing the unnormalised baseline,
#' and `SpaNorm`, representing spatially aware normalisation.
#'
#' The original cell-type and tissue-region annotations
#' are retained for downstream benchmarking.
#'
#' @name norm_example
#' @docType data
#' @format A named list of SpatialExperiment objects.
#' @source Derived from the STARmap PLUS mouse brain dataset
#'   available through the Bioconductor STexampleData package.
NULL


#' Example spatially variable gene results
#'
#' Precomputed spatially variable gene detection results
#' generated from `norm_example` using the STBench
#' `run_svg()` function with the nnSVG method.
#'
#' Results are provided for the normalisation methods
#' `None` and `SpaNorm`.
#'
#' This dataset allows users to evaluate and visualise
#' spatially variable gene statistics and rankings
#' without repeating the SVG detection procedure.
#'
#' @name svg_example
#' @docType data
#' @format A named list of SpatialExperiment objects
#'   containing nnSVG results.
#' @source Derived from the STARmap PLUS mouse brain dataset
#'   available through the Bioconductor STexampleData package.
NULL


#' Example normalised data for cell-cell communication
#'
#' A larger subset of the publicly available STARmap PLUS
#' mouse brain dataset, prepared for demonstrating
#' cell-cell communication and Gene Ontology analysis.
#'
#' An initial random subset of 300 genes and 1000 cells
#' was selected using seed 123. Genes and cells with
#' zero total counts were removed.
#'
#' Gene symbols were harmonised with the mouse CellChat
#' database where unambiguous matches were available.
#'
#' The subset was normalised using `normalize_spe()`
#' with the `None` and `SpaNorm` methods.
#'
#' Original cell-type and tissue-region annotations
#' are retained for downstream analyses.
#'
#' @name norm_ccc_example
#' @docType data
#' @format A named list of SpatialExperiment objects.
#' @source Bioconductor STexampleData package:
#'   \url{https://bioconductor.org/packages/STexampleData/}
NULL
