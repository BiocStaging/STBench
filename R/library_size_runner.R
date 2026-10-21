
# ============================================================
# Normalised library-size comparison
# ============================================================


#' Compare normalised library sizes across methods
#'
#' Calculate and visualise total normalised expression
#' across spatial locations for different normalisation methods.
#'
#' Log-normalised expression is back-transformed using expm1(),
#' summed across genes, and transformed again using log1p().
#'
#' Produces a boxplot and spatial maps.
#'
#' By default, sctransform and Giotto are excluded because
#' their transformed expression values may not represent
#' log1p-transformed non-negative expression.
#'
#' @param norm Named list of normalised SpatialExperiment objects.
#' @param methods Character vector of methods to compare.
#'   If NULL, all methods except sctransform and Giotto are used.
#' @param assay_name Name of the log-normalised assay.
#'   Default is "logcounts".
#' @param point_size Point size for spatial maps.
#' @param theme_name Plot theme.
#'
#' @return A list containing comparison data, boxplot,
#'   spatial plot, and combined plot.
#'
#' @export
plot_library_size_comparison <- function(
    norm,
    methods = NULL,
    assay_name = "logcounts",
    point_size = 0.5,
    theme_name = "classic"
) {

  # ----------------------------------------------------------
  # Select methods
  # ----------------------------------------------------------

  excluded <- c(
    "sct",
    "sctransform",
    "giotto"
  )

  if (
    is.null(names(norm)) ||
    anyNA(names(norm)) ||
    any(!nzchar(names(norm))) ||
    anyDuplicated(names(norm))
  ) {
    stop("norm must have unique, non-empty method names.")
  }

  if (is.null(methods)) {

    methods <- names(norm)[
      !tolower(names(norm)) %in% excluded
    ]

  }

  unknown <- setdiff(
    methods,
    names(norm)
  )

  if (length(unknown) > 0) {

    stop(
      "Unknown methods: ",
      paste(unknown, collapse = ", ")
    )
  }

  if (length(methods) == 0) {
    stop("No methods selected.")
  }

  methods <- unique(methods)


  # ----------------------------------------------------------
  # Calculate normalised library sizes
  # ----------------------------------------------------------

  comparison_data <- lapply(
    methods,
    function(method) {

      spe <- norm[[method]]

      if (
        !assay_name %in%
        SummarizedExperiment::assayNames(spe)
      ) {

        stop(
          "Assay '",
          assay_name,
          "' not found for method: ",
          method
        )
      }


      # Extract log-normalised expression

      log_mat <- SummarizedExperiment::assay(
        spe,
        assay_name
      )


      # Back-transform the log-normalised expression

      norm_mat <- expm1(log_mat)


      # Calculate the total normalised expression per spot

      expression_total <- MatrixGenerics::colSums2(
        norm_mat
      )


      # Apply log1p to the normalised expression totals

      log_total <- log1p(
        expression_total
      )


      # Extract spatial coordinates

      coords <- SpatialExperiment::spatialCoords(spe)

      if (ncol(coords) < 2) {

        stop(
          "Two spatial coordinates are required for: ",
          method
        )
      }


      data.frame(
        method = method,
        spot = seq_len(ncol(spe)),
        x = coords[, 1],
        y = coords[, 2],
        expression_total = expression_total,
        log_expression_total = log_total
      )
    }
  )


  # ----------------------------------------------------------
  # Combine data
  # ----------------------------------------------------------

  df <- dplyr::bind_rows(
    comparison_data
  )

  df$method <- factor(
    df$method,
    levels = methods
  )


  # ----------------------------------------------------------
  # Boxplot comparison
  # ----------------------------------------------------------

  p_boxplot <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$method,
      y = .data$log_expression_total,
      fill = .data$method
    )
  ) +
    ggplot2::geom_boxplot(
      width = 0.65,
      alpha = 0.8,
      colour = "black",
      outlier.size = 0.5
    ) +
    ggplot2::labs(
      x = "Normalisation method",
      y = "Log total normalised expression"
    ) +
    .get_plot_theme(theme_name) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 30,
        hjust = 1
      )
    )


  # ----------------------------------------------------------
  # Spatial maps
  # ----------------------------------------------------------

  p_spatial <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$x,
      y = .data$y,
      colour = .data$log_expression_total
    )
  ) +
    ggplot2::geom_point(
      size = point_size
    ) +
    ggplot2::facet_wrap(
      ~ method
    ) +
    ggplot2::coord_fixed() +
    viridis::scale_colour_viridis(
      option = "D",
      name = "Log \nexpression"
    ) +
    ggplot2::labs(
      x = " X",
      y = " Y"
    ) +
    .get_plot_theme(theme_name) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold"
      )
    )


  # ----------------------------------------------------------
  # Combine plots
  # ----------------------------------------------------------

  p_combined <- patchwork::wrap_plots(
    p_boxplot,
    p_spatial,
    ncol = 1
  )


  # ----------------------------------------------------------
  # Return
  # ----------------------------------------------------------

  return(
    list(
      data = df,
      boxplot = p_boxplot,
      spatial = p_spatial,
      plot = p_combined
    )
  )
}
