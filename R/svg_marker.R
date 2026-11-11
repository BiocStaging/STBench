
# SVG marker gene analysis



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

    # y_lab <- paste0(
    #   "log(",
    #   stat,
    #   ")"
    # )
    y_lab <- "log (LR statistic)"

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



#' Plot SVG Rank Against SVG Statistic
#'
#' Visualises the relationship between SVG statistics and
#' SVG ranks for selected genes across normalisation methods.
#'
#' Each point represents a normalisation method, with
#' genes displayed in separate facets.
#'
#' @param svg_norm Named list of SVG results.
#' @param genes Genes to plot.
#' @param stat SVG statistic. Default is "LR_stat".
#' @param method_order Optional method order.
#' @param method_colors Optional named vector of method colours.
#' @param transform Transformation: log1p, log, or none.
#' @param scales Facet scales: fixed, free, free_x, or free_y.
#' @param ncol Number of facet columns.
#' @param nrow Number of facet rows.
#' @param point_size Point size.
#' @param gene_id Gene identifier.
#' @param theme_name Plot theme.
#'
#' @return A list containing data and plot.
#'
#' @export
plot_svg_rank_vs_stat <- function(
    svg_norm,
    genes,
    stat = "LR_stat",
    method_order = NULL,
    method_colors = NULL,
    transform = c("log1p", "log", "none"),
    scales = c("fixed", "free", "free_x", "free_y"),
    ncol = NULL,
    nrow = NULL,
    point_size = 3,
    gene_id = "rownames",
    theme_name = "classic"
) {

  transform <- match.arg(transform)
  scales <- match.arg(scales)

  # Extract SVG statistics

  df_stat <- extract_svg_stats(
    svg_norm = svg_norm,
    genes = genes,
    stat = stat,
    gene_id = gene_id,
    method_order = method_order
  )

  # Extract SVG ranks

  df_rank <- extract_svg_stats(
    svg_norm = svg_norm,
    genes = genes,
    stat = "rank",
    gene_id = gene_id,
    method_order = method_order
  )

  # Combine statistics and ranks

  df <- dplyr::inner_join(
    df_stat,
    df_rank[, c("method", "gene", "svg_stat")],
    by = c("method", "gene"),
    suffix = c("_stat", "_rank")
  )

  names(df)[names(df) == "svg_stat_rank"] <- "rank"

  # Apply transformation

  if (transform == "log1p") {

    if (any(df$svg_stat_stat < -1, na.rm = TRUE)) {
      stop("log1p requires SVG statistics >= -1.")
    }

    df$value <- log1p(df$svg_stat_stat)

    x_lab <- paste0("log1p(", stat, ")")

  } else if (transform == "log") {

    df <- df[
      !is.na(df$svg_stat_stat) &
        df$svg_stat_stat > 0,
      ,
      drop = FALSE
    ]

    df$value <- log(df$svg_stat_stat)

    x_lab <- paste0("log(", stat, ")")

  } else {

    df$value <- df$svg_stat_stat
    x_lab <- stat
  }

  # Remove missing or non-finite values

  df <- df[
    is.finite(df$value) &
      is.finite(df$rank),
    ,
    drop = FALSE
  ]

  # Plot

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$value,
      y = .data$rank,
      colour = .data$method
    )
  ) +

    ggplot2::geom_point(
      size = point_size
    ) +

    ggplot2::facet_wrap(
      ~ gene,
      scales = scales,
      ncol = ncol,
      nrow = nrow
    ) +

    ggplot2::scale_y_reverse() +

    ggplot2::labs(
      x = x_lab,
      y = "SVG rank",
      colour = "Normalisation method"
    ) +

    .get_plot_theme(theme_name) +

    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      panel.grid = ggplot2::element_blank()
    )

  # Optional method colours

  if (!is.null(method_colors)) {

    p <- p +
      ggplot2::scale_colour_manual(
        values = method_colors
      )
  }

  # Return

  list(
    data = df,
    plot = p
  )
}





#' Plot SVG Rank Correlation Across Normalisation Methods
#'
#' Calculates pairwise Spearman rank correlations of SVG
#' rankings across normalisation methods.
#'
#' Only genes shared by all selected methods are used.
#' The heatmap displays the lower triangle of the
#' correlation matrix, including the diagonal.
#'
#' @param svg_norm Named list of SVG results.
#' @param methods Optional methods to compare.
#' @param gene_id Gene identifier. Default is "rownames".
#' @param method_order Optional order of normalisation methods.
#' @param method_colors Three colours for low, middle, and high correlations.
#' @param show_values Display correlation values. Default TRUE.
#' @param digits Number of decimal places.
#' @param theme_name Plot theme.
#'
#' @return A list containing the correlation matrix,
#'   heatmap data, shared gene ranks, and ggplot object.
#'
#' @export
plot_svg_rank_correlation <- function(
    svg_norm,
    methods = NULL,
    gene_id = "rownames",
    method_order = NULL,
    method_colors = c(
      "#B2182B",
      "#F7F7F7",
      "#2166AC"
    ),
    show_values = TRUE,
    digits = 2,
    theme_name = "classic"
) {

  # Validate method names

  if (
    is.null(names(svg_norm)) ||
    anyNA(names(svg_norm)) ||
    any(!nzchar(names(svg_norm))) ||
    anyDuplicated(names(svg_norm))
  ) {
    stop(
      "svg_norm must have unique, non-empty method names."
    )
  }

  # Select methods

  if (is.null(methods)) {
    methods <- names(svg_norm)
  }

  methods <- unique(methods)

  if (!is.null(method_order)) {

    if (
      length(method_order) != length(methods) ||
      !setequal(method_order, methods)
    ) {
      stop(
        "method_order must contain all selected methods."
      )
    }

    methods <- method_order
  }

  if (
    length(methods) < 2L ||
    any(!methods %in% names(svg_norm))
  ) {
    stop(
      "At least two valid normalisation methods are required."
    )
  }

  # Extract SVG ranks

  rank_list <- lapply(
    methods,
    function(method) {

      result <- svg_norm[[method]]

      if (gene_id == "rownames") {
        genes <- rownames(result)
      } else {
        genes <- result[[gene_id]]
      }

      if (
        is.null(genes) ||
        anyNA(genes) ||
        anyDuplicated(genes)
      ) {
        stop(
          "Valid unique gene identifiers are required for: ",
          method
        )
      }

      df <- extract_svg_stats(
        svg_norm = svg_norm[method],
        genes = genes,
        stat = "rank",
        gene_id = gene_id
      )

      data.frame(
        gene = as.character(df$gene),
        rank = as.numeric(df$svg_stat)
      )
    }
  )

  names(rank_list) <- methods

  # Merge rankings using shared genes

  rank_wide <- rank_list[[1]]
  names(rank_wide)[2] <- methods[1]

  for (i in seq_along(methods)[-1]) {

    temp <- rank_list[[i]]
    names(temp)[2] <- methods[i]

    rank_wide <- merge(
      rank_wide,
      temp,
      by = "gene",
      all = FALSE,
      sort = FALSE
    )
  }

  # Keep genes with valid ranks across all methods

  rank_matrix <- as.matrix(
    rank_wide[, methods, drop = FALSE]
  )

  valid <- rowSums(!is.finite(rank_matrix)) == 0L

  rank_wide <- rank_wide[
    valid,
    ,
    drop = FALSE
  ]

  if (nrow(rank_wide) < 2L) {
    stop(
      "At least two genes with valid ranks are required."
    )
  }

  # Spearman rank correlation

  correlation_matrix <- stats::cor(
    rank_wide[, methods, drop = FALSE],
    method = "spearman"
  )

  # Convert correlation matrix to long format

  correlation_data <- as.data.frame(
    as.table(correlation_matrix)
  )

  names(correlation_data) <- c(
    "method_x",
    "method_y",
    "correlation"
  )

  # Keep lower triangle including diagonal

  correlation_data$row_id <- match(
    correlation_data$method_y,
    methods
  )

  correlation_data$col_id <- match(
    correlation_data$method_x,
    methods
  )

  correlation_data <- correlation_data[
    correlation_data$row_id >= correlation_data$col_id,
    ,
    drop = FALSE
  ]

  # Set method order

  correlation_data$method_x <- factor(
    correlation_data$method_x,
    levels = methods
  )

  correlation_data$method_y <- factor(
    correlation_data$method_y,
    levels = rev(methods)
  )

  # Plot heatmap

  p <- ggplot2::ggplot(
    correlation_data,
    ggplot2::aes(
      x = .data$method_x,
      y = .data$method_y,
      fill = .data$correlation
    )
  ) +

    ggplot2::geom_tile(
      colour = "white",
      linewidth = 0.5,alpha=0.8
    ) +

    ggplot2::scale_fill_viridis_c(
      name = "Spearman\ncorrelation"
    ) +

    ggplot2::coord_fixed() +

    ggplot2::labs(
      x = NULL,
      y = NULL
    ) +

    .get_plot_theme(theme_name) +

    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      axis.ticks = ggplot2::element_blank()
    )

  # Add correlation coefficients

  if (show_values) {

    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(
          label = sprintf(
            paste0("%.", digits, "f"),
            .data$correlation
          )
        ),
        size = 3.5
      )
  }

  # Return

  list(
    correlation_matrix = correlation_matrix,
    data = correlation_data,
    ranks = rank_wide,
    plot = p
  )
}


#' Compare Marker Gene Expression Between Cell Groups
#'
#' Performs Welch's two-sample t-test for selected marker
#' genes, comparing expression in the specified cell group
#' against all other cells across normalisation methods.
#'
#' @param norm Named list of normalised SpatialExperiment objects.
#' @param genes Character vector of marker genes.
#' @param marker_labels Expected cell type or region for each gene.
#'   Must have the same length as genes.
#' @param region_col Column name in colData containing cell
#'   type or spatial region annotations.
#' @param assay_name Expression assay. Default is "logcounts".
#' @param method_order Optional order of normalisation methods.
#' @param method_colors Optional named vector of method colours.
#' @param point_size Point size. Default is 3.
#' @param scales Facet scales: fixed, free, free_x, or free_y.
#' @param ncol Number of facet columns.
#' @param nrow Number of facet rows.
#' @param theme_name Plot theme.
#'
#' @return A list containing t-test results and a ggplot object.
#'
#' @export
plot_marker_ttest <- function(
    norm,
    genes,
    marker_labels,
    region_col = "cell_type",
    assay_name = "logcounts",
    method_order = NULL,
    method_colors = NULL,
    point_size = 3,
    scales = "free_y",
    ncol = NULL,
    nrow = NULL,
    theme_name = "classic"
) {

  scales <- match.arg(
    scales,
    c("fixed", "free", "free_x", "free_y")
  )

  # Validate inputs

  if (length(genes) != length(marker_labels)) {
    stop("genes and marker_labels must have the same length.")
  }

  if (length(genes) == 0L) {
    stop("At least one marker gene is required.")
  }

  if (
    is.null(names(norm)) ||
    anyNA(names(norm)) ||
    any(!nzchar(names(norm))) ||
    anyDuplicated(names(norm))
  ) {
    stop("norm must have unique, non-empty method names.")
  }

  methods <- names(norm)

  if (!is.null(method_order)) {

    if (!setequal(method_order, methods) ||
        length(method_order) != length(methods)) {
      stop("method_order must contain all methods exactly once.")
    }

    methods <- method_order
  }

  # Perform t-tests

  results <- lapply(methods, function(method) {

    spe <- norm[[method]]

    if (!region_col %in% colnames(
      SummarizedExperiment::colData(spe)
    )) {
      stop(
        "Column '", region_col,
        "' not found for method: ", method
      )
    }

    if (!assay_name %in%
        SummarizedExperiment::assayNames(spe)) {
      stop(
        "Assay '", assay_name,
        "' not found for method: ", method
      )
    }

    region <- as.character(
      SummarizedExperiment::colData(spe)[[region_col]]
    )

    expr_mat <- SummarizedExperiment::assay(
      spe,
      assay_name
    )

    dplyr::bind_rows(
      lapply(seq_along(genes), function(i) {

        gene <- genes[i]
        marker <- marker_labels[i]

        if (!gene %in% rownames(expr_mat)) {
          return(NULL)
        }

        expression <- as.numeric(expr_mat[gene, ])

        valid <- is.finite(expression) & !is.na(region)

        group_yes <- expression[
          valid & region == marker
        ]

        group_no <- expression[
          valid & region != marker
        ]

        n_yes <- length(group_yes)
        n_no <- length(group_no)

        if (n_yes < 2 || n_no < 2) {
          return(NULL)
        }

        tt <- tryCatch(
          stats::t.test(group_yes, group_no),
          error = function(e) NULL
        )

        if (is.null(tt)) {
          return(NULL)
        }

        mean_yes <- mean(group_yes)
        mean_no <- mean(group_no)

        data.frame(
          method = method,
          gene = gene,
          marker_label = marker,
          n_yes = n_yes,
          n_no = n_no,
          mean_yes = mean_yes,
          mean_no = mean_no,
          mean_diff = mean_yes - mean_no,
          t_stat = unname(tt$statistic),
          abs_t = abs(unname(tt$statistic)),
          p_value = tt$p.value,
          correct_direction = mean_yes > mean_no
        )
      })
    )
  })

  df <- dplyr::bind_rows(results)

  if (nrow(df) == 0L) {
    stop("No valid marker t-tests could be performed.")
  }

  # Adjust p-values across all performed tests

  df$p_adj <- stats::p.adjust(
    df$p_value,
    method = "BH"
  )

  # Set factor order

  df$method <- factor(
    df$method,
    levels = methods
  )

  df$gene <- factor(
    df$gene,
    levels = unique(genes)
  )

  # Plot absolute t-statistics

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$method,
      y = .data$abs_t,
      colour = .data$method
    )
  ) +

    ggplot2::geom_point(
      size = point_size
    ) +

    ggplot2::facet_wrap(
      ~ gene,
      scales = scales,
      ncol = ncol,
      nrow = nrow
    ) +

    ggplot2::labs(
      x =NULL,# "Method",
      y = "|t-statistic|"
    ) +

    .get_plot_theme(theme_name) +

    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      legend.position = "none",
      panel.grid = ggplot2::element_blank()
    )

  if (!is.null(method_colors)) {
    p <- p +
      ggplot2::scale_colour_manual(
        values = method_colors
      )
  }

  return(
    list(
      data = df,
      plot = p
    )
  )
}

