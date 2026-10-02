# ============================================================
# Individual normalisation methods
# ============================================================


#' No normalisation
#'
#' Applies log1p transformation to raw counts.
#'
#' @param spe A SpatialExperiment or SingleCellExperiment object.
#'
#' @return A modified object containing logcounts.
#' @export
norm_none <- function(spe) {

  start_time <- Sys.time()

  SingleCellExperiment::logcounts(spe) <-
    log1p(SingleCellExperiment::counts(spe))

  end_time <- Sys.time()

  message(
    "norm_none completed in ",
    round(
      difftime(
        end_time,
        start_time,
        units = "mins"
      ),
      3
    ),
    " mins"
  )

  spe
}


#' Scran normalisation
#'
#' @param spe A SpatialExperiment or SingleCellExperiment object.
#'
#' @return Normalised object.
#' @export
norm_scran <- function(spe) {

  start_time <- Sys.time()

  cl <- scran::quickCluster(spe)

  spe <- scran::computeSumFactors(
    spe,
    clusters = cl
  )

  SingleCellExperiment::sizeFactors(spe) <-
    pmax(
      1e-8,
      SingleCellExperiment::sizeFactors(spe)
    )

  spe <- scater::logNormCounts(spe)

  end_time <- Sys.time()

  message(
    "norm_scran completed in ",
    round(
      difftime(
        end_time,
        start_time,
        units = "mins"
      ),
      3
    ),
    " mins"
  )

  spe
}


#' SCTransform normalisation
#'
#' @param spe A SpatialExperiment or SingleCellExperiment object.
#' @param seed Random seed.
#' @param vst_flavor SCTransform vst flavour.
#' @param verbosity Verbosity passed to sctransform.
#' @param ... Additional arguments passed to sctransform::vst().
#'
#' @return Normalised object.
#' @export
norm_sct <- function(
    spe,
    seed = 123,
    vst_flavor = "v2",
    verbosity = 1,
    ...
) {

  start_time <- Sys.time()

  set.seed(seed)

  vst_res <- sctransform::vst(
    SingleCellExperiment::counts(spe),
    vst.flavor = vst_flavor,
    verbosity = verbosity,
    batch_var = NULL,
    cell_attr = as.data.frame(
      SummarizedExperiment::colData(spe)
    ),
    ...
  )

  spe <- spe[
    rownames(vst_res$y),
    ,
    drop = FALSE
  ]

  SingleCellExperiment::logcounts(spe) <-
    vst_res$y

  end_time <- Sys.time()

  message(
    "norm_sct completed in ",
    round(
      difftime(
        end_time,
        start_time,
        units = "mins"
      ),
      3
    ),
    " mins"
  )

  spe
}


#' Giotto normalisation
#'
#' @param spe A SpatialExperiment object.
#' @param ... Additional arguments passed to Giotto normalisation.
#'
#' @return Normalised object.
#' @export
norm_giotto <- function(
    spe,
    ...
) {

  start_time <- Sys.time()

  giotto_obj <- Giotto::createGiottoObject(
    raw_exprs = SingleCellExperiment::counts(spe),
    spatial_locs =
      SpatialExperiment::spatialCoords(spe)
  )

  giotto_obj <- Giotto::normalizeGiotto(
    giotto_obj,
    ...
  )

  SingleCellExperiment::logcounts(spe) <-
    giotto_obj@norm_scaled_expr

  end_time <- Sys.time()

  message(
    "norm_giotto completed in ",
    round(
      difftime(
        end_time,
        start_time,
        units = "mins"
      ),
      3
    ),
    " mins"
  )

  spe
}


#' SpaNorm normalisation
#'
#' @param spe A SpatialExperiment object.
#' @param sample_p Sampling proportion. If NULL, 3000/ncol(spe)
#'   is used.
#' @param df_tps Degrees of freedom for thin-plate spline.
#' @param adj_method Adjustment method.
#' @param seed Random seed.
#' @param ... Additional arguments passed to SpaNorm.
#'
#' @return Normalised object.
#' @export
norm_spanorm <- function(
    spe,
    sample_p = NULL,
    df_tps = 6,
    adj_method = "logpac",
    seed = 123,
    ...
) {

  start_time <- Sys.time()

  if (!is.null(SingleCellExperiment::sizeFactors(spe))) {
    SingleCellExperiment::sizeFactors(spe) <- NULL
  }

  if (is.null(sample_p)) {
    sample_p <- min(
      1,
      3000 / ncol(spe)
    )
  }

  set.seed(seed)

  spe <- SpaNorm::SpaNorm(
    spe,
    sample.p = sample_p,
    df.tps = df_tps,
    adj.method = adj_method,
    ...
  )

  end_time <- Sys.time()

  message(
    "norm_spanorm completed in ",
    round(
      difftime(
        end_time,
        start_time,
        units = "mins"
      ),
      3
    ),
    " mins"
  )

  spe
}
