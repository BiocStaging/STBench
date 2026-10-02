# ============================================================
# SVG detection
# ============================================================


#' Run spatially variable gene detection
#'
#' @param norm_list Named list of normalised SpatialExperiment objects.
#' @param method SVG detection method.
#' @param ... Additional parameters passed to the SVG method.
#'
#' @return Named list containing SVG results.
#'
#' @export
run_svg <- function(
    norm_list,
    method = "nnSVG",
    ...
) {

  .check_named_list(norm_list)

  if (method != "nnSVG") {

    stop(
      "Currently supported SVG method: nnSVG"
    )
  }

  out <- lapply(
    norm_list,
    function(spe) {

      nnSVG::nnSVG(
        spe,
        ...
      )
    }
  )

  names(out) <- names(norm_list)

  out
}
