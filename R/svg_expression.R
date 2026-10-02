# ============================================================
# Spatial expression plots
# ============================================================


#' Plot spatial expression across normalisation methods
#'
#' Plot spatial expression of selected genes across multiple
#' normalisation methods.
#'
#' When `scale_expr = TRUE`, expression is scaled from 0 to 1
#' separately for each gene across all normalisation methods.
#' The resulting maps are displayed using a facet grid with
#' methods in rows and genes in columns.
#'
#' When `scale_expr = FALSE`, the original expression values
#' from the selected assay are retained. A separate plot is
#' created for each gene, with normalisation methods shown as
#' facets. Therefore, each gene has its own colour scale.
#'
#' @param norm_list Named list of normalised SpatialExperiment
#'   objects.
#' @param genes Character vector containing genes to plot.
#' @param assay_name Name of assay containing expression values.
#'   Default is `"logcounts"`.
#' @param method_order Optional character vector specifying the
#'   order of normalisation methods.
#' @param scale_expr Logical. If TRUE, expression is scaled
#'   between 0 and 1 within each gene across all methods.
#' @param point_size Point size for spatial locations.
#' @param colour_option Viridis colour option.
#' @param ncol_genes Number of columns used when combining
#'   gene-specific unscaled plots.
#' @param theme_name Plot theme. One of `"classic"`, `"bw"`,
#'   `"minimal"`, or `"gray"`.
#'
#' @return A list containing the plotting data and plot.
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
    theme_name = "classic"
) {

  # ----------------------------------------------------------
  # Resolve methods
  # ----------------------------------------------------------

  methods <- .resolve_methods(
    norm_list,
    method_order
  )

  if (length(methods) == 0) {
    stop("No valid normalisation methods were found.")
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
        paste(
          method,
          gene,
          sep = "_"
        )
      ]] <- temp
    }
  }


  # ----------------------------------------------------------
  # Check extracted data
  # ----------------------------------------------------------

  if (length(df_list) == 0) {

    stop(
      "None of the requested genes were available in the supplied objects."
    )
  }


  # ----------------------------------------------------------
  # Combine
  # ----------------------------------------------------------

  df <- do.call(
    rbind,
    df_list
  )

  rownames(df) <- NULL


  # ----------------------------------------------------------
  # Factor ordering
  # ----------------------------------------------------------

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
  #
  # Scale each gene across ALL methods.
  #
  # Therefore:
  # PLP1 across none/scran/SCT/etc uses one 0-1 scale.
  # MBP gets its own 0-1 scaling.
  #
  # This allows comparison of the same gene across methods.
  # ==========================================================

  if (scale_expr) {

    df_scaled <- df |>
      dplyr::group_by(
        .data$gene
      ) |>
      dplyr::mutate(

        expression_scaled =
          .scale_zero_one(
            .data$expression
          )

      ) |>
      dplyr::ungroup()


    # --------------------------------------------------------
    # Facet grid:
    #
    # columns = genes
    # rows    = methods
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
        limits = c(
          0,
          1
        ),
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

        panel.grid =
          ggplot2::element_blank(),

        strip.background =
          ggplot2::element_blank(),

        strip.text =
          ggplot2::element_text(
            size = 9
          ),

        axis.text =
          ggplot2::element_text(
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
  #
  # Original assay values are retained.
  #
  # Each gene gets a SEPARATE plot.
  # Methods are wrapped within each gene.
  #
  # Therefore:
  # - PLP1 has one colour scale across methods
  # - MBP has another colour scale across methods
  # - GFAP has another colour scale across methods
  #
  # Methods remain directly comparable within each gene.
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
    # Expression limits for this gene across all methods
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

        panel.grid =
          ggplot2::element_blank(),

        strip.background =
          ggplot2::element_blank(),

        strip.text =
          ggplot2::element_text(
            size = 9
          ),

        axis.text =
          ggplot2::element_text(
            size = 8
          ),

        plot.title =
          ggplot2::element_text(
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
