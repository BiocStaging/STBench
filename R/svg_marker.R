# ============================================================
# SVG marker gene analysis
# ============================================================


#' Extract SVG statistics for selected genes
#'
#' @param svg_norm Named list of SVG results.
#' @param genes Character vector of genes.
#' @param stat SVG statistic column.
#' @param gene_id Gene identifier.
#' @param method_order Optional method order.
#'
#' @return Data frame.
#' @export
extract_svg_stats <- function(
    svg_norm,
    genes,
    stat = "LR_stat",
    gene_id = "rownames",
    method_order = NULL
) {

  methods <- .resolve_methods(
    svg_norm,
    method_order
  )

  out <- lapply(
    methods,
    function(method) {

      spe <- svg_norm[[method]]

      rd <- .extract_rowdata(
        spe,
        gene_id
      )

      if (!stat %in% colnames(rd)) {

        stop(
          "Statistic '",
          stat,
          "' not found for method ",
          method
        )
      }

      rd <- rd[
        rd$gene %in% genes,
        c(
          "gene",
          stat
        ),
        drop = FALSE
      ]

      names(rd)[2] <- "svg_stat"

      rd$method <- method

      rd
    }
  )

  df <- do.call(
    rbind,
    out
  )

  rownames(df) <- NULL

  df$method <- factor(
    df$method,
    levels = methods
  )

  df$gene <- factor(
    df$gene,
    levels = genes
  )

  df
}


#' Plot SVG statistic for selected genes
#'
#' @param svg_norm Named list of SVG results.
#' @param genes Genes to plot.
#' @param stat SVG statistic.
#' @param method_order Optional method order.
#' @param method_colors Optional named vector of method colours.
#' @param transform Transformation: log1p, log, or none.
#' @param point_size Point size.
#' @param gene_id Gene identifier.
#' @param theme_name Plot theme.
#'
#' @return A list containing data and plot.
#' @export
plot_svg_stat <- function(
    svg_norm,
    genes,
    stat = "LR_stat",
    method_order = NULL,
    method_colors = NULL,
    transform = c(
      "log1p",
      "log",
      "none"
    ),
    point_size = 3,
    gene_id = "rownames",
    theme_name = "classic"
) {

  transform <- match.arg(transform)

  df <- extract_svg_stats(
    svg_norm = svg_norm,
    genes = genes,
    stat = stat,
    gene_id = gene_id,
    method_order = method_order
  )

  if (transform == "log1p") {

    df$value <- log1p(
      df$svg_stat
    )

    y_lab <- paste0(
      "log1p(",
      stat,
      ")"
    )

  } else if (transform == "log") {

    df <- df[
      df$svg_stat > 0,
      ,
      drop = FALSE
    ]

    df$value <- log(
      df$svg_stat
    )

    y_lab <- paste0(
      "log(",
      stat,
      ")"
    )

  } else {

    df$value <- df$svg_stat
    y_lab <- stat
  }

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = method,
      y = value,
      colour = method
    )
  ) +
    ggplot2::geom_point(
      size = point_size
    ) +
    ggplot2::facet_wrap(
      ~ gene,
      scales = "free_y"
    ) +
    ggplot2::labs(
      x = "Normalisation method",
      y = y_lab
    ) +
    .get_plot_theme(
      theme_name
    ) +
    ggplot2::theme(
      axis.text.x =
        ggplot2::element_text(
          angle = 45,
          hjust = 1
        ),
      strip.background =
        ggplot2::element_blank(),
      legend.position = "none"
    )

  if (!is.null(method_colors)) {

    p <- p +
      ggplot2::scale_colour_manual(
        values = method_colors
      )
  }

  list(
    data = df,
    plot = p
  )
}


#' Plot SVG ranks for selected genes
#'
#' @param svg_norm Named list of SVG results.
#' @param genes Genes to plot.
#' @param method_order Optional method order.
#' @param method_colors Optional method colours.
#' @param point_size Point size.
#' @param gene_id Gene identifier.
#' @param theme_name Plot theme.
#'
#' @return List containing data and plot.
#' @export
plot_svg_rank <- function(
    svg_norm,
    genes,
    method_order = NULL,
    method_colors = NULL,
    point_size = 3,
    gene_id = "rownames",
    theme_name = "classic"
) {

  df <- extract_svg_stats(
    svg_norm = svg_norm,
    genes = genes,
    stat = "rank",
    gene_id = gene_id,
    method_order = method_order
  )

  names(df)[
    names(df) == "svg_stat"
  ] <- "rank"

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = method,
      y = rank,
      colour = method
    )
  ) +
    ggplot2::geom_point(
      size = point_size
    ) +
    ggplot2::facet_wrap(
      ~ gene
    ) +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(
      x = "Normalisation method",
      y = "SVG rank"
    ) +
    .get_plot_theme(
      theme_name
    ) +
    ggplot2::theme(
      axis.text.x =
        ggplot2::element_text(
          angle = 45,
          hjust = 1
        ),
      strip.background =
        ggplot2::element_blank(),
      legend.position = "none"
    )

  if (!is.null(method_colors)) {

    p <- p +
      ggplot2::scale_colour_manual(
        values = method_colors
      )
  }

  list(
    data = df,
    plot = p
  )
}
