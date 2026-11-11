# BANKSY clustering evaluation utilities for an R package

#' Calculate BANKSY clustering evaluation metrics
#'
#' Calculates ARI, NMI, homogeneity, completeness, and spatial average
#' silhouette width (ASW) for each normalisation method and spatial slice.
#' Overall is the unweighted mean of NMI, HOM, and COM; it is a user-defined
#' composite agreement score rather than a standard accuracy measure.
#'
#' ASW uses Euclidean distances between the first two spatial coordinates,
#' not expression-space distances. For large slices, ASW is estimated from
#' a random subset of cells to limit memory requirements.
#'
#' @param norm_list Named list of SpatialExperiment objects, or a named list
#'   whose elements are named lists of SpatialExperiment objects (slices).
#' @param ground_truth_col Name of the colData column containing reference labels.
#' @param cluster_col Name of the colData column containing predicted clusters.
#' @param methods Optional vector of normalisation methods to evaluate.
#' @param nmi_variant NMI normalisation; one of "sum", "max", "min", "sqrt",
#'   or "joint". "sum" uses arithmetic-mean entropy normalization.
#' @param asw_max_cells Maximum number of cells used to calculate spatial ASW
#'   per slice. Larger slices are randomly subsampled.
#' @param seed Integer random seed for reproducible ASW subsampling.
#'
#' @return Data frame with method, slice, n, ARI, NMI, HOM, COM, ASW,
#'   Overall, n_asw, and ASW_sampled.
#' @export
calculate_clustering_metrics <- function(
    norm_list,
    ground_truth_col,
    cluster_col = "banksy_cluster",
    methods = names(norm_list),
    nmi_variant = "sum",
    asw_max_cells = 3000L,
    seed = 123L
) {
  if (!is.list(norm_list) || is.null(names(norm_list)) ||
      anyNA(names(norm_list)) || any(!nzchar(names(norm_list))) ||
      anyDuplicated(names(norm_list))) {
    stop("'norm_list' must be a list with unique, nonempty method names.")
  }
  if (!is.character(ground_truth_col) || length(ground_truth_col) != 1L ||
      is.na(ground_truth_col) || !nzchar(ground_truth_col)) {
    stop("'ground_truth_col' must be a nonempty column name.")
  }
  if (!is.character(cluster_col) || length(cluster_col) != 1L ||
      is.na(cluster_col) || !nzchar(cluster_col)) {
    stop("'cluster_col' must be a nonempty column name.")
  }
  nmi_variant <- match.arg(nmi_variant, c("sum", "max", "min", "sqrt", "joint"))
  if (length(asw_max_cells) != 1L || is.na(asw_max_cells) ||
      !is.finite(asw_max_cells) || asw_max_cells < 3L) {
    stop("'asw_max_cells' must be at least 3.")
  }
  if (!is.character(methods) || anyNA(methods) ||
      anyDuplicated(methods) || !length(methods) ||
      any(!methods %in% names(norm_list))) {
    stop("'methods' must contain distinct, available method names.")
  }
  results <- list()
  index <- 1L

  for (method in methods) {
    method_data <- norm_list[[method]]
    if (methods::is(method_data, "SpatialExperiment")) {
      slices <- list(slice1 = method_data)
    } else if (is.list(method_data) && length(method_data) &&
               !is.null(names(method_data)) &&
               !anyNA(names(method_data)) &&
               all(nzchar(names(method_data))) &&
               !anyDuplicated(names(method_data))) {
      slices <- method_data
    } else {
      stop("Invalid or unnamed slices in method: ", method)
    }

    for (slice_name in names(slices)) {
      spe <- slices[[slice_name]]
      if (!methods::is(spe, "SpatialExperiment")) {
        stop("Expected SpatialExperiment: ", method, " / ", slice_name)
      }
      cd <- SummarizedExperiment::colData(spe)
      missing_cols <- setdiff(c(ground_truth_col, cluster_col), colnames(cd))
      if (length(missing_cols)) {
        stop("Missing colData columns for ", method, " / ", slice_name,
             ": ", paste(missing_cols, collapse = ", "))
      }
      truth <- as.character(cd[[ground_truth_col]])
      predicted <- as.character(cd[[cluster_col]])
      coords <- as.matrix(SpatialExperiment::spatialCoords(spe))
      if (nrow(coords) != length(truth) || ncol(coords) < 2L) {
        stop("Invalid spatial coordinates for ", method, " / ", slice_name)
      }
      coords <- coords[, 1:2, drop = FALSE]
      storage.mode(coords) <- "double"
      keep <- !is.na(truth) & nzchar(truth) &
        !is.na(predicted) & nzchar(predicted) &
        stats::complete.cases(coords) & rowSums(is.finite(coords)) == 2L
      truth <- truth[keep]
      predicted <- predicted[keep]
      coords <- coords[keep, , drop = FALSE]
      n <- length(truth)
      if (!n) {
        warning("No complete observations: ", method, " / ", slice_name,
                call. = FALSE)
        next
      }

      ari <- mclust::adjustedRandIndex(truth, predicted)
      nmi <- aricode::NMI(factor(truth), factor(predicted), variant = nmi_variant)
      hom <- clevr::homogeneity(truth, predicted)
      com <- clevr::completeness(truth, predicted)
      overall <- mean(c(nmi, hom, com))

      n_asw <- n
      sampled <- n > asw_max_cells
      if (sampled) {
        set.seed(seed)
        ids <- sample.int(n, asw_max_cells)
        asw_coords <- coords[ids, , drop = FALSE]
        asw_clusters <- predicted[ids]
        n_asw <- length(ids)
      } else {
        asw_coords <- coords
        asw_clusters <- predicted
      }
      k <- length(unique(asw_clusters))
      if (k < 2L || k >= n_asw) {
        asw <- NA_real_
        warning("ASW undefined (needs 2 <= number of clusters < n): ",
                method, " / ", slice_name, call. = FALSE)
      } else {
        spatial_dist <- stats::dist(asw_coords, method = "euclidean")
        sil <- cluster::silhouette(as.integer(factor(asw_clusters)), spatial_dist)
        asw <- mean(sil[, "sil_width"])
      }
      results[[index]] <- data.frame(
        method = method, slice = slice_name, n = n,
        ARI = ari, NMI = nmi, HOM = hom, COM = com,
        ASW = asw, Overall = overall,
        n_asw = n_asw, ASW_sampled = sampled,
        stringsAsFactors = FALSE
      )
      index <- index + 1L
    }
  }
  if (!length(results)) stop("No metrics could be calculated.")
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}


#' Plot clustering metrics across normalisation methods
#'
#' Each facet represents one clustering metric.
#' Points are coloured by normalisation method.
#'
#' @param metrics_df Data frame from calculate_clustering_metrics().
#' @param metrics Metrics to plot.
#' @param point_size Point size.
#' @param point_alpha Point opacity.
#' @param ncol Number of facet columns.
#' @param scales Facet scales.
#' @param palette Optional named colours for normalisation methods.
#' @param theme_name One of "classic", "bw", "minimal", "gray".
#' @param rotate_x Angle of method labels.
#'
#' @return List containing data and ggplot object (`plot`).
#' @export
plot_clustering_metrics <- function(
    metrics_df,
    metrics = c("ARI", "NMI", "HOM", "COM", "ASW", "Overall"),
    point_size = 4,
    point_alpha = 0.9,
    ncol = 3,
    scales = "free_y",
    palette = NULL,
    theme_name = "classic",
    rotate_x = 45
) {

  # Check input
  required <- c("method", metrics)

  if (!all(required %in% names(metrics_df))) {
    stop("Required columns missing from metrics_df.")
  }

  scales <- match.arg(
    scales,
    c("free_y", "fixed", "free", "free_x")
  )

  theme_name <- match.arg(
    theme_name,
    c("classic", "bw", "minimal", "gray")
  )

  # Prepare plotting data
  df <- do.call(rbind, lapply(metrics, function(metric) {

    data.frame(
      method = as.character(metrics_df$method),
      metric = metric,
      score = as.numeric(metrics_df[[metric]]),
      stringsAsFactors = FALSE
    )

  }))

  df$method <- factor(
    df$method,
    levels = unique(as.character(metrics_df$method))
  )

  df$metric <- factor(
    df$metric,
    levels = metrics
  )

  # Plot
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$method,
      y = .data$score,
      colour = .data$method
    )
  ) +

    ggplot2::geom_point(
      size = point_size,
      alpha = point_alpha
    ) +

    ggplot2::facet_wrap(
      ~metric,
      ncol = ncol,
      scales = scales
    ) +

    ggplot2::labs(
      x = NULL,#"Normalisation method",
      y = "Score",
      colour = "Method"
    ) +

    .get_plot_theme(theme_name) +

    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = rotate_x,
        hjust = 1
      ),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    )

  # Optional colour palette
  if (!is.null(palette)) {
    p <- p +
      ggplot2::scale_colour_manual(
        values = palette
      )
  }

  # Return
  list(
    data = df,
    plot = p
  )
}


#' Rank normalisation methods by overall clustering agreement
#'
#' Methods are ranked using the Overall score, calculated
#' as the mean of NMI, HOM and COM.
#'
#' @param metrics_df Data frame from calculate_clustering_metrics().
#'
#' @return Data frame containing Rank, method and Overall.
#' @export
rank_normalization_methods <- function(
    metrics_df
) {

  # Check input
  required <- c("method", "Overall")

  if (!all(required %in% names(metrics_df))) {
    stop("Required columns missing from metrics_df.")
  }

  if (anyDuplicated(metrics_df$method)) {
    stop("Expected one row per normalisation method.")
  }

  if (anyNA(metrics_df$Overall)) {
    stop("Overall scores cannot be missing.")
  }

  # Rank methods
  out <- metrics_df[
    order(-metrics_df$Overall),
    c("method", "Overall"),
    drop = FALSE
  ]

  out$Rank <- seq_len(nrow(out))

  out <- out[, c("Rank", "method", "Overall")]

  rownames(out) <- NULL

  # Return
  out
}


#' Plot normalisation method ranking
#'
#' Display overall clustering agreement scores across
#' normalisation methods.
#'
#' @param ranking_df Data frame from rank_normalization_methods().
#' @param bar_width Bar width.
#' @param palette Optional named colours for normalisation methods.
#' @param label_size Numeric label size.
#' @param show_values Whether to show score labels.
#' @param theme_name One of "classic", "bw", "minimal", "gray".
#'
#' @return List containing data and ggplot object (`plot`).
#' @export
plot_normalization_ranking <- function(
    ranking_df,
    bar_width = 0.7,
    palette = NULL,
    label_size = 3.5,
    show_values = TRUE,
    theme_name = "classic"
) {

  # Check input
  required <- c("method", "Overall")

  if (!all(required %in% names(ranking_df))) {
    stop("Required columns missing from ranking_df.")
  }

  # Prepare plotting data
  df <- ranking_df[
    order(-ranking_df$Overall),
    ,
    drop = FALSE
  ]

  df$method <- factor(
    df$method,
    levels = rev(as.character(df$method))
  )

  # Plot
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$Overall,
      y = .data$method,
      fill = .data$method
    )
  ) +

    ggplot2::geom_col(
      width = bar_width
    ) +

    ggplot2::coord_cartesian(
      xlim = c(0, 1),
      clip = "off"
    ) +

    ggplot2::labs(
      x = "Overall clustering agreement score",
      y = NULL,#"Normalisation method",
      fill = "Method"
    ) +

    .get_plot_theme(theme_name) +

    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "none"
    )

  # Optional colour palette
  if (!is.null(palette)) {
    p <- p +
      ggplot2::scale_fill_manual(
        values = palette
      )
  }

  # Optional score labels
  if (show_values) {

    p <- p +

      ggplot2::geom_text(
        ggplot2::aes(
          label = sprintf("%.3f", .data$Overall)
        ),
        hjust = -0.15,
        size = label_size
      ) +

      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(
          mult = c(0, 0.12)
        )
      )
  }

  # Return
  list(
    data = df,
    plot = p
  )
}
