
# Normalisation runner



#' Apply multiple normalisation methods
#'
#' @param spe A SpatialExperiment object.
#' @param methods Character vector containing normalisation methods.
#' @param params Named list containing method-specific parameters.
#'
#' @return A named list of normalised objects.
#'
#' @examples
#' \dontrun{
#' norm <- normalize_spe(
#'   spe,
#'   methods = c("none", "scran")
#' )
#' }
#'
#' @export
normalize_spe <- function(
    spe,
    methods = c(
      "none",
      "scran",
      "sct",
      "spanorm",
      "giotto"
    ),
    params = list()
) {

  norm_functions <- list(
    None = norm_none,
    scran = norm_scran,
    sctransform = norm_sct,
    SpaNorm = norm_spanorm,
    Giotto = norm_giotto
  )

  unknown_methods <- setdiff(
    methods,
    names(norm_functions)
  )

  if (length(unknown_methods) > 0) {

    stop(
      "Unknown normalisation method(s): ",
      paste(
        unknown_methods,
        collapse = ", "
      )
    )
  }

  out <- lapply(
    methods,
    function(method) {

      fun <- norm_functions[[method]]

      method_params <- params[[method]]

      if (is.null(method_params)) {
        method_params <- list()
      }

      do.call(
        fun,
        c(
          list(spe = spe),
          method_params
        )
      )
    }
  )

  names(out) <- methods

  out
}
