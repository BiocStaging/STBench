# ============================================================
# SVG gene-set analyses
# ============================================================


#' Extract SVG statistics for a gene set
#'
#' @param svg_norm Named SVG results.
#' @param genes Character vector containing genes.
#' @param stat SVG statistic.
#' @param gene_id Gene identifier.
#' @param method_order Optional method order.
#'
#' @return Data frame.
#' @export
extract_svg_gene_set <- function(
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
          "' not found for method: ",
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

  df
}


#' Plot SVG statistic distribution for a gene set
#'
#' @param svg_norm Named SVG results.
#' @param genes Gene set.
#' @param stat SVG statistic.
#' @param method_order Optional method order.
#' @param method_colors Optional named colours.
#' @param gene_id Gene identifier.
#' @param transform Transformation.
#' @param plot_type boxplot or violin.
#' @param add_jitter Add gene-level points.
#' @param point_size Point size.
#' @param alpha Point transparency.
#' @param theme_name Plot theme.
#'
#' @return List containing data and plot.
#'
#' @export
plot_svg_gene_set_distribution <- function(
    svg_norm,
    genes,
    stat = "LR_stat",
    method_order = NULL,
    method_colors = NULL,
    gene_id = "rownames",
    transform = c(
      "log1p",
      "log",
      "none"
    ),
    plot_type = c(
      "boxplot",
      "violin"
    ),
    add_jitter = TRUE,
    point_size = 1.5,
    alpha = 0.6,
    theme_name = "classic"
) {

  transform <- match.arg(
    transform
  )

  plot_type <- match.arg(
    plot_type
  )

  methods <- .resolve_methods(
    svg_norm,
    method_order
  )

  df <- extract_svg_gene_set(
    svg_norm = svg_norm,
    genes = genes,
    stat = stat,
    gene_id = gene_id,
    method_order = methods
  )

  if (nrow(df) == 0) {

    stop(
      "None of the requested genes were found."
    )
  }

  df$method <- factor(
    df$method,
    levels = methods
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
      fill = method,
      colour = method
    )
  )

  if (plot_type == "boxplot") {

    p <- p +
      ggplot2::geom_boxplot(
        alpha = 0.65,
        outlier.shape = NA
      )

  } else {

    p <- p +
      ggplot2::geom_violin(
        alpha = 0.65,
        trim = FALSE
      )
  }

  if (add_jitter) {

    p <- p +
      ggplot2::geom_jitter(
        width = 0.15,
        size = point_size,
        alpha = alpha
      )
  }

  p <- p +
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
      legend.position = "none"
    )

  if (!is.null(method_colors)) {

    p <- p +
      ggplot2::scale_fill_manual(
        values = method_colors
      ) +
      ggplot2::scale_colour_manual(
        values = method_colors
      )
  }

  list(
    data = df,
    plot = p
  )
}


#' Proportion of a gene set among top-ranked SVGs
#'
#' @param svg_norm Named SVG results.
#' @param genes Gene set, such as housekeeping genes.
#' @param top_n Vector containing top-N values.
#' @param method_order Optional method order.
#' @param method_colors Optional named colours.
#' @param gene_id Gene identifier.
#' @param theme_name Plot theme.
#'
#' @return List containing data and plot.
#'
#' @export
plot_svg_topn_gene_set <- function(
    svg_norm,
    genes,
    top_n = c(
      50,
      100,
      200,
      300
    ),
    method_order = NULL,
    method_colors = NULL,
    gene_id = "rownames",
    theme_name = "classic"
) {

  methods <- .resolve_methods(
    svg_norm,
    method_order
  )

  top_n <- sort(
    unique(
      as.integer(
        top_n
      )
    )
  )

  top_n <- top_n[
    top_n > 0
  ]

  df <- lapply(
    methods,
    function(method) {

      spe <- svg_norm[[method]]

      rd <- .extract_rowdata(
        spe,
        gene_id
      )

      if (!"rank" %in% colnames(rd)) {

        stop(
          "'rank' not found for method: ",
          method
        )
      }

      rd <- rd[
        !is.na(rd$rank) &
          !is.na(rd$gene),
        ,
        drop = FALSE
      ]

      rd <- rd[
        order(rd$rank),
        ,
        drop = FALSE
      ]

      n_use <- top_n[
        top_n <= nrow(rd)
      ]

      if (length(n_use) == 0) {
        return(NULL)
      }

      do.call(
        rbind,
        lapply(
          n_use,
          function(n) {

            top_genes <- rd$gene[
              seq_len(n)
            ]

            data.frame(
              method = method,
              top_n = n,
              n_gene_set =
                sum(
                  top_genes %in% genes
                ),
              proportion =
                mean(
                  top_genes %in% genes
                )
            )
          }
        )
      )
    }
  )

  df <- Filter(
    Negate(is.null),
    df
  )

  if (length(df) == 0) {

    stop(
      "No valid top-N values were available."
    )
  }

  df <- do.call(
    rbind,
    df
  )

  rownames(df) <- NULL

  df$method <- factor(
    df$method,
    levels = methods
  )

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = top_n,
      y = proportion,
      colour = method,
      group = method
    )
  ) +
    ggplot2::geom_line(
      linewidth = 1
    ) +
    ggplot2::geom_point(
      size = 2
    ) +
    ggplot2::labs(
      x = "Top N SVGs",
      y = "Gene-set proportion",
      colour = "Method"
    ) +
    .get_plot_theme(
      theme_name
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
