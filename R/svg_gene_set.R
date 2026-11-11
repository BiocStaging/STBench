
# SVG gene-set analyses



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
      x = NULL,#"Normalisation method",
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




# SVG gene-set proportion among top-ranked genes



#' Proportion of a gene set among top-ranked SVGs
#'
#' Calculate and visualise the proportion of genes belonging
#' to a user-defined gene set among the top-ranked spatially
#' variable genes (SVGs) across normalisation methods.
#'
#' The gene set may contain housekeeping genes, negative-control
#' genes, marker genes, or any other user-defined gene set.
#'
#' @param svg_norm Named list of SVG results.
#' @param genes Character vector containing genes of interest.
#' @param top_n Numeric vector of top-N SVG thresholds.
#' @param method_order Optional normalisation method order.
#' @param method_colors Optional named vector of method colours.
#' @param gene_id Gene identifier used in SVG results.
#' @param theme_name Plot theme.
#' @param gene_set_name Name of the gene set.
#' @param title Optional plot title.
#' @param subtitle Optional plot subtitle.
#' @param alpha Transparency of points and lines.
#' @param point_size Point size.
#' @param line_width Line width.
#' @param line_type Line type.
#' @param show_points Logical. Display points at each top-N value.
#' @param show_legend Logical. Display the method legend.
#' @param percentage Logical. Display y-axis as percentages.
#'
#' @return List containing data and plot.
#'
#' @export
plot_svg_topn_gene_set <- function(
    svg_norm,
    genes,
    top_n = c(50, 100, 200, 300),
    method_order = NULL,
    method_colors = NULL,
    gene_id = "rownames",
    theme_name = "classic",
    gene_set_name = "Gene set",
    title = NULL,
    subtitle = NULL,
    alpha = 1,
    point_size = 2,
    line_width = 1,
    line_type = "solid",
    show_points = TRUE,
    show_legend = TRUE,
    percentage = FALSE
) {

  # ----------------------------------------------------------
  # Resolve methods
  # ----------------------------------------------------------

  methods <- .resolve_methods(
    svg_norm,
    method_order
  )


  # ----------------------------------------------------------
  # Validate gene set and top-N values
  # ----------------------------------------------------------

  genes <- unique(
    as.character(genes)
  )

  genes <- genes[
    !is.na(genes) & nzchar(genes)
  ]

  if (length(genes) == 0) {
    stop("The supplied gene set is empty.")
  }

  if (
    !is.numeric(top_n) ||
    anyNA(top_n) ||
    any(!is.finite(top_n)) ||
    any(top_n <= 0) ||
    any(top_n != floor(top_n))
  ) {
    stop("top_n must contain positive integers.")
  }

  top_n <- sort(
    unique(as.integer(top_n))
  )

  if (length(top_n) == 0) {
    stop("At least one top-N value is required.")
  }


  # ----------------------------------------------------------
  # Extract top-ranked SVGs
  # ----------------------------------------------------------

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
          is.finite(rd$rank) &
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

            n_gene_set <- sum(
              top_genes %in% genes
            )

            data.frame(
              method = method,
              gene_set = gene_set_name,
              top_n = n,
              n_gene_set = n_gene_set,
              proportion = n_gene_set / n,
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )


  # ----------------------------------------------------------
  # Combine data
  # ----------------------------------------------------------

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


  # ----------------------------------------------------------
  # Axis labels
  # ----------------------------------------------------------

  y_label <- paste0(
    gene_set_name,
    " proportion"
  )

  if (is.null(title)) {

    title <- paste0(
      gene_set_name,
      " among top-ranked SVGs"
    )
  }


  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$top_n,
      y = .data$proportion,
      colour = .data$method,
      group = .data$method
    )
  ) +
    ggplot2::geom_line(
      linewidth = line_width,
      linetype = line_type,
      alpha = alpha
    )


  if (show_points) {

    p <- p +
      ggplot2::geom_point(
        size = point_size,
        alpha = alpha
      )
  }


  p <- p +
    ggplot2::scale_x_continuous(
      breaks = top_n
    ) +
    ggplot2::labs(
      title = NULL,#title,
      subtitle = subtitle,
      x = "Top N SVGs",
      y = y_label,
      colour = "Method"
    ) +
    .get_plot_theme(
      theme_name
    ) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      legend.position = if (show_legend) {
        "right"
      } else {
        "none"
      }
    )


  # ----------------------------------------------------------
  # Percentage axis
  # ----------------------------------------------------------

  if (percentage) {

    p <- p +
      ggplot2::scale_y_continuous(
        labels = scales::label_percent(
          accuracy = 1
        ),
        expand = ggplot2::expansion(
          mult = c(0.02, 0.08)
        )
      )
  }


  # ----------------------------------------------------------
  # Method colours
  # ----------------------------------------------------------

  if (!is.null(method_colors)) {

    p <- p +
      ggplot2::scale_colour_manual(
        values = method_colors
      )
  }


  # ----------------------------------------------------------
  # Return
  # ----------------------------------------------------------

  return(
    list(
      data = df,
      plot = p
    )
  )
}





# SVG comparison across multiple gene sets



#' Compare SVG statistics across multiple gene sets
#'
#' Plot distributions of spatial variability statistics
#' for multiple user-defined gene sets across normalisation
#' methods.
#'
#' Gene sets may include housekeeping genes, negative-control
#' genes, biological marker genes, or any custom gene sets.
#'
#' Each gene set is displayed in a separate facet.
#'
#' @param svg_norm Named list of SVG results.
#' @param gene_sets Named list of character vectors containing
#'   genes for each gene set.
#' @param stat SVG statistic. Default is "LR_stat".
#' @param method_order Optional normalisation method order.
#' @param method_colors Optional named vector of method colours.
#' @param gene_id Gene identifier.
#' @param transform Transformation: "log1p", "log", or "none".
#' @param plot_type Plot type: "boxplot" or "violin".
#' @param scales Facet scales: "fixed", "free_y", or "free".
#' @param facet_ncol Number of facet columns.
#' @param add_jitter Logical. Display individual gene points.
#' @param point_size Size of jittered points.
#' @param alpha Transparency of boxplots or violins.
#' @param jitter_alpha Transparency of jittered points.
#' @param box_width Width of boxplots.
#' @param title Optional plot title.
#' @param subtitle Optional plot subtitle.
#' @param x_label X-axis label.
#' @param y_label Optional y-axis label.
#' @param show_legend Logical. Display method legend.
#' @param theme_name Plot theme.
#'
#' @return List containing plotting data and ggplot.
#'
#' @export
plot_svg_gene_sets <- function(
    svg_norm,
    gene_sets,
    stat = "LR_stat",
    method_order = NULL,
    method_colors = NULL,
    gene_id = "rownames",
    transform = c("log1p", "log", "none"),
    plot_type = c("boxplot", "violin"),
    scales = c("free_y", "fixed", "free"),
    facet_ncol = 2,
    add_jitter = TRUE,
    point_size = 1,
    alpha = 0.65,
    jitter_alpha = 0.5,
    box_width = 0.65,
    title = NULL,
    subtitle = NULL,
    x_label = "Normalisation method",
    y_label = NULL,
    show_legend = FALSE,
    theme_name = "classic"
) {

  # ----------------------------------------------------------
  # Validate input
  # ----------------------------------------------------------

  transform <- match.arg(transform)

  plot_type <- match.arg(plot_type)

  scales <- match.arg(scales)

  if (
    !is.list(gene_sets) ||
    length(gene_sets) == 0 ||
    is.null(names(gene_sets)) ||
    anyNA(names(gene_sets)) ||
    any(!nzchar(names(gene_sets))) ||
    anyDuplicated(names(gene_sets))
  ) {
    stop(
      "gene_sets must be a named list of gene vectors."
    )
  }

  if (!all(vapply(
    gene_sets,
    is.character,
    logical(1)
  ))) {
    stop("Each gene set must be a character vector.")
  }

  methods <- .resolve_methods(
    svg_norm,
    method_order
  )


  # ----------------------------------------------------------
  # Extract SVG statistics
  # ----------------------------------------------------------

  df_list <- lapply(
    names(gene_sets),
    function(set_name) {

      genes <- unique(
        gene_sets[[set_name]]
      )

      genes <- genes[
        !is.na(genes) & nzchar(genes)
      ]

      if (length(genes) == 0) {
        return(NULL)
      }

      df <- extract_svg_gene_set(
        svg_norm = svg_norm,
        genes = genes,
        stat = stat,
        gene_id = gene_id,
        method_order = methods
      )

      if (nrow(df) == 0) {
        return(NULL)
      }

      df$gene_set <- set_name

      df
    }
  )

  df_list <- Filter(
    Negate(is.null),
    df_list
  )

  if (length(df_list) == 0) {
    stop(
      "None of the requested genes were found."
    )
  }

  df <- dplyr::bind_rows(df_list)


  # ----------------------------------------------------------
  # Apply transformations
  # ----------------------------------------------------------

  if (transform == "log1p") {

    if (any(df$svg_stat < -1, na.rm = TRUE)) {
      stop("log1p requires SVG statistics >= -1.")
    }

    df$value <- log1p(
      df$svg_stat
    )

    default_y_label <- paste0(
      "log1p(",
      stat,
      ")"
    )

  } else if (transform == "log") {

    df <- df[
      !is.na(df$svg_stat) &
        df$svg_stat > 0,
      ,
      drop = FALSE
    ]

    df$value <- log(
      df$svg_stat
    )

    default_y_label <- paste0(
      "log(",
      stat,
      ")"
    )

  } else {

    df$value <- df$svg_stat

    default_y_label <- stat
  }

  df <- df[
    is.finite(df$value),
    ,
    drop = FALSE
  ]

  if (nrow(df) == 0) {
    stop("No finite SVG statistics available.")
  }


  # ----------------------------------------------------------
  # Factor ordering
  # ----------------------------------------------------------

  df$method <- factor(
    df$method,
    levels = methods
  )

  df$gene_set <- factor(
    df$gene_set,
    levels = names(gene_sets)
  )

  if (is.null(y_label)) {
    y_label <- default_y_label
  }


  # ----------------------------------------------------------
  # Initialise plot
  # ----------------------------------------------------------

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$method,
      y = .data$value,
      fill = .data$method,
      colour = .data$method
    )
  )


  # ----------------------------------------------------------
  # Boxplot or violin
  # ----------------------------------------------------------

  if (plot_type == "boxplot") {

    p <- p +
      ggplot2::geom_boxplot(
        width = box_width,
        alpha = alpha,
        outlier.shape = NA,
        colour = "black"
      )

  } else {

    p <- p +
      ggplot2::geom_violin(
        alpha = alpha,
        trim = FALSE
      )
  }


  # ----------------------------------------------------------
  # Individual gene points
  # ----------------------------------------------------------

  if (add_jitter) {

    p <- p +
      ggplot2::geom_jitter(
        width = 0.15,
        height = 0,
        size = point_size,
        alpha = jitter_alpha
      )
  }


  # ----------------------------------------------------------
  # Facets and labels
  # ----------------------------------------------------------

  p <- p +
    ggplot2::facet_wrap(
      ~ gene_set,
      ncol = facet_ncol,
      scales = scales
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = y_label,
      fill = "Normalisation method",
      colour = "Normalisation method"
    ) +
    .get_plot_theme(
      theme_name
    ) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      legend.position = if (show_legend) {
        "right"
      } else {
        "none"
      }
    )


  # ----------------------------------------------------------
  # Method colours
  # ----------------------------------------------------------

  if (!is.null(method_colors)) {

    p <- p +
      ggplot2::scale_fill_manual(
        values = method_colors
      ) +
      ggplot2::scale_colour_manual(
        values = method_colors
      )
  }


  # ----------------------------------------------------------
  # Return
  # ----------------------------------------------------------

  return(
    list(
      data = df,
      plot = p
    )
  )
}
