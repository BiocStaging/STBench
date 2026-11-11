
# Internal utility functions



.extract_rowdata <- function(
    spe,
    gene_id = "rownames"
) {

  rd <- as.data.frame(
    SummarizedExperiment::rowData(spe)
  )

  if (gene_id == "rownames") {

    rd$gene <- rownames(rd)

  } else {

    if (!gene_id %in% colnames(rd)) {

      stop(
        "Gene identifier column not found: ",
        gene_id
      )
    }

    rd$gene <- rd[[gene_id]]
  }

  rd
}


.check_named_list <- function(x) {

  if (!is.list(x)) {
    stop("Input must be a list.")
  }

  if (is.null(names(x))) {
    stop("Input list must be named.")
  }

  invisible(TRUE)
}


.resolve_methods <- function(
    x,
    method_order = NULL
) {

  .check_named_list(x)

  available_methods <- names(x)[
    !vapply(
      x,
      is.null,
      logical(1)
    )
  ]

  if (is.null(method_order)) {
    return(available_methods)
  }

  missing_methods <- setdiff(
    method_order,
    available_methods
  )

  if (length(missing_methods) > 0) {

    warning(
      "Methods not available and skipped: ",
      paste(
        missing_methods,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  intersect(
    method_order,
    available_methods
  )
}


.scale_zero_one <- function(x) {

  xmin <- min(
    x,
    na.rm = TRUE
  )

  xmax <- max(
    x,
    na.rm = TRUE
  )

  if (
    !is.finite(xmin) ||
    !is.finite(xmax) ||
    xmax == xmin
  ) {

    return(
      rep(
        0,
        length(x)
      )
    )
  }

  (x - xmin) /
    (xmax - xmin)
}
