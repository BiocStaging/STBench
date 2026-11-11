

# Spatial expression plots



#' Plot spatial expression across normalisation methods
#'
#' Visualise spatial expression of selected genes across
#' multiple normalisation methods.
#'
#' Expression can be scaled between 0 and 1 either across
#' all selected methods for each gene or independently
#' within each gene and normalisation method.
#'
#' When scale_expr = FALSE, original assay values are used.
#' Each gene receives a separate plot with a common colour
#' scale across the selected methods.
#'
#' @param norm_list Named list of normalised SpatialExperiment
#'   objects.
#' @param genes Character vector of genes to plot.
#' @param methods Optional character vector of normalisation
#'   methods. If NULL, all available methods are used.
#' @param assay_name Name of the expression assay.
#' @param method_order Optional ordering of selected methods.
#' @param scale_expr Logical. Whether to scale expression.
#' @param scale_by Scaling approach: "gene" or "gene_method".
#' @param point_size Spatial point size.
#' @param colour_option Viridis colour option.
#' @param ncol_genes Number of columns for unscaled plots.
#' @param theme_name Plot theme.
#'
#' @return A list containing plotting data, individual
#'   gene plots when unscaled, and the combined plot.
#'
#' @export
plot_spatial_expression <- function(
    norm_list,
    genes,
    assay_name = "logcounts",
    method_order = NULL,
    scale_expr = TRUE,
    point_size = 0.5,
    colour_option = "magma",
    ncol_genes = 1,
    theme_name = "classic",
    methods = NULL,
    scale_by = c("gene", "gene_method")
) {

  # ----------------------------------------------------------
  # Resolve selected methods
  # ----------------------------------------------------------

  scale_by <- match.arg(scale_by)

  if (is.null(methods)) {
    methods <- names(norm_list)
  }

  if (length(methods) == 0) {
    stop("No normalisation methods were selected.")
  }

  unknown_methods <- setdiff(
    methods,
    names(norm_list)
  )

  if (length(unknown_methods) > 0) {
    stop(
      "Unknown normalisation methods: ",
      paste(unknown_methods, collapse = ", ")
    )
  }

  methods <- unique(methods)

  if (!is.null(method_order)) {

    unknown_order <- setdiff(
      method_order,
      names(norm_list)
    )

    if (length(unknown_order) > 0) {
      stop(
        "Unknown methods in method_order: ",
        paste(unknown_order, collapse = ", ")
      )
    }

    methods <- c(
      intersect(method_order, methods),
      setdiff(methods, method_order)
    )
  }


  # ----------------------------------------------------------
  # Extract expression and coordinates
  # ----------------------------------------------------------

  df_list <- list()

  for (method in methods) {

    spe <- norm_list[[method]]


    # --------------------------------------------------------
    # Check assay
    # --------------------------------------------------------

    if (
      !assay_name %in%
      SummarizedExperiment::assayNames(spe)
    ) {

      warning(
        "Assay '",
        assay_name,
        "' not found for method: ",
        method,
        call. = FALSE
      )

      next
    }


    # --------------------------------------------------------
    # Genes available in this object
    # --------------------------------------------------------

    genes_use <- intersect(
      genes,
      rownames(spe)
    )

    if (length(genes_use) == 0) {
      next
    }


    # --------------------------------------------------------
    # Spatial coordinates
    # --------------------------------------------------------

    coords <- as.data.frame(
      SpatialExperiment::spatialCoords(spe)
    )

    if (ncol(coords) < 2) {
      stop(
        "spatialCoords(spe) must contain at least two columns."
      )
    }


    # --------------------------------------------------------
    # Expression matrix
    # --------------------------------------------------------

    expr_mat <- SummarizedExperiment::assay(
      spe,
      assay_name
    )


    # --------------------------------------------------------
    # Extract each gene
    # --------------------------------------------------------

    for (gene in genes_use) {

      temp <- data.frame(
        x = coords[, 1],
        y = coords[, 2],
        expression = as.numeric(
          expr_mat[
            gene,
            ,
            drop = TRUE
          ]
        ),
        gene = gene,
        method = method,
        stringsAsFactors = FALSE
      )

      df_list[[
        paste(method, gene, sep = "_")
      ]] <- temp
    }
  }


  # ----------------------------------------------------------
  # Check extracted data
  # ----------------------------------------------------------

  if (length(df_list) == 0) {
    stop(
      "None of the requested genes were available ",
      "in the selected methods."
    )
  }


  # ----------------------------------------------------------
  # Combine
  # ----------------------------------------------------------

  df <- dplyr::bind_rows(
    df_list
  )

  df$method <- factor(
    df$method,
    levels = methods
  )

  df$gene <- factor(
    df$gene,
    levels = genes
  )


  # ==========================================================
  # SCALED EXPRESSION
  # ==========================================================

  if (scale_expr) {

    # --------------------------------------------------------
    # Scale each gene across all selected methods
    # --------------------------------------------------------

    if (scale_by == "gene") {

      df_scaled <- df |>
        dplyr::group_by(
          .data$gene
        ) |>
        dplyr::mutate(
          expression_scaled = .scale_zero_one(
            .data$expression
          )
        ) |>
        dplyr::ungroup()

    }


    # --------------------------------------------------------
    # Scale each gene independently within each method
    # --------------------------------------------------------

    if (scale_by == "gene_method") {

      df_scaled <- df |>
        dplyr::group_by(
          .data$gene,
          .data$method
        ) |>
        dplyr::mutate(
          expression_scaled = .scale_zero_one(
            .data$expression
          )
        ) |>
        dplyr::ungroup()

    }


    # --------------------------------------------------------
    # Faceted spatial expression maps
    # --------------------------------------------------------

    p <- ggplot2::ggplot(
      df_scaled,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        colour = .data$expression_scaled
      )
    ) +
      ggplot2::geom_point(
        size = point_size
      ) +
      ggplot2::facet_grid(
        method ~ gene
      ) +
      ggplot2::coord_fixed() +
      viridis::scale_colour_viridis(
        option = colour_option,
        limits = c(0, 1),
        name = "Scaled\nexpression"
      ) +
      ggplot2::labs(
        x = "x",
        y = "y"
      ) +
      .get_plot_theme(
        theme_name
      ) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(
          size = 9
        ),
        axis.text = ggplot2::element_text(
          size = 8
        ),
        legend.position = "right"
      )


    return(
      list(
        data = df_scaled,
        plot = p
      )
    )
  }


  # ==========================================================
  # ORIGINAL EXPRESSION
  # ==========================================================

  gene_plots <- list()

  for (gene_use in genes) {

    gene_df <- df[
      as.character(df$gene) == gene_use,
      ,
      drop = FALSE
    ]

    if (nrow(gene_df) == 0) {
      next
    }


    # --------------------------------------------------------
    # Gene-specific colour scale across selected methods
    # --------------------------------------------------------

    gene_limits <- range(
      gene_df$expression,
      na.rm = TRUE
    )


    # --------------------------------------------------------
    # Gene-specific plot
    # --------------------------------------------------------

    p_gene <- ggplot2::ggplot(
      gene_df,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        colour = .data$expression
      )
    ) +
      ggplot2::geom_point(
        size = point_size
      ) +
      ggplot2::facet_wrap(
        ~ method,
        nrow = 1
      ) +
      ggplot2::coord_fixed() +
      viridis::scale_colour_viridis(
        option = colour_option,
        limits = gene_limits,
        name = "Expression"
      ) +
      ggplot2::labs(
        title = gene_use,
        x = "x",
        y = "y"
      ) +
      .get_plot_theme(
        theme_name
      ) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(
          size = 9
        ),
        axis.text = ggplot2::element_text(
          size = 8
        ),
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold"
        ),
        legend.position = "right"
      )

    gene_plots[[gene_use]] <- p_gene
  }


  # ----------------------------------------------------------
  # Combine separate gene plots
  # ----------------------------------------------------------

  combined_plot <- patchwork::wrap_plots(
    gene_plots,
    ncol = ncol_genes
  )


  return(
    list(
      data = df,
      plots = gene_plots,
      plot = combined_plot
    )
  )
}
