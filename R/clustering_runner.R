# Spatial clustering plot


#' Plot spatial clustering results
#'
#' Visualise predicted spatial clusters and ground-truth
#' annotations across normalisation methods.
#'
#' When show = "both", the ground truth is displayed once
#' per spatial slice, followed by predicted cluster maps
#' across all selected normalisation methods.
#'
#' Predicted clusters share one colour scale and legend.
#' Ground truth uses a separate palette and legend.
#'
#' @param norm_list Named list of SpatialExperiment objects
#'   or nested lists of SpatialExperiment objects.
#' @param cluster_col Column containing predicted clusters.
#' @param methods Optional normalisation methods to plot.
#' @param method_order Optional method order.
#' @param ground_truth_col Ground-truth annotation column.
#' @param show One of "both", "predicted", or "truth".
#' @param point_size Spatial point size.
#' @param palette Optional named colours for predicted clusters.
#' @param truth_palette Optional named colours for ground truth.
#' @param ncol Number of columns for spatial facets.
#' @param legend_point_size Point size in legends.
#' @param theme_name Plot theme.
#'
#' @return List containing data, palettes, individual plots,
#'   and combined plot.
#'
#' @examples
#' data(norm_example)
#'
#' banksy_norm <- run_banksy_individual(
#'   norm_list = norm_example,
#'   k_geom = 15,
#'   lambda = 0.8,
#'   resolution = 0.2,
#'   seed = 123
#' )
#'
#' p <- plot_spatial_clustering(
#'   norm_list = banksy_norm,
#'   cluster_col = "banksy_cluster",
#'   ground_truth_col = "layer_guess",
#'   show = "both",
#'   point_size = 0.5,
#'   ncol = 2,
#'   theme_name = "classic"
#' )
#'
#' p$plot
#' @export
plot_spatial_clustering <- function(
  norm_list,
  cluster_col = "banksy_cluster",
  methods = names(norm_list),
  method_order = NULL,
  point_size = 0.5,
  palette = NULL,
  ncol = 3,
  legend_point_size = 4,
  theme_name = "classic",
  ground_truth_col = NULL,
  show = c("predicted", "both", "truth"),
  truth_palette = NULL
) {
  # Check input


  .check_named_list(norm_list)

  show <- match.arg(show)

  if (
    show %in% c("both", "truth") &&
      is.null(ground_truth_col)
  ) {
    stop(
      "ground_truth_col must be supplied when ",
      "show = 'both' or 'truth'."
    )
  }


  # Resolve methods


  missing_methods <- setdiff(
    methods,
    names(norm_list)
  )

  if (length(missing_methods) > 0) {
    stop(
      "Unknown normalisation methods: ",
      paste(missing_methods, collapse = ", ")
    )
  }

  methods <- unique(methods)

  if (!is.null(method_order)) {
    methods <- c(
      intersect(method_order, methods),
      setdiff(methods, method_order)
    )
  }

  if (length(methods) == 0) {
    stop("No normalisation methods selected.")
  }


  # Extract spatial annotations


  df_list <- list()

  for (method in methods) {
    method_data <- norm_list[[method]]


    # Single SPE or multiple slices


    if (
      methods::is(
        method_data,
        "SpatialExperiment"
      )
    ) {
      spe_list <- list(
        slice1 = method_data
      )
    } else if (is.list(method_data)) {
      spe_list <- method_data

      if (
        is.null(names(spe_list)) ||
          anyNA(names(spe_list)) ||
          any(!nzchar(names(spe_list)))
      ) {
        stop(
          "Slice lists must be named for method: ",
          method
        )
      }
    } else {
      stop(
        "Invalid input for method: ",
        method
      )
    }


    # Process each slice


    for (slice_name in names(spe_list)) {
      spe <- spe_list[[slice_name]]

      if (
        !methods::is(
          spe,
          "SpatialExperiment"
        )
      ) {
        stop(
          "Invalid SpatialExperiment for: ",
          method,
          " / ",
          slice_name
        )
      }

      cd <- SummarizedExperiment::colData(spe)


      # Check required columns


      if (
        show %in% c("both", "predicted") &&
          !cluster_col %in% colnames(cd)
      ) {
        stop(
          "Cluster column '",
          cluster_col,
          "' not found for: ",
          method
        )
      }

      if (
        show %in% c("both", "truth") &&
          !ground_truth_col %in% colnames(cd)
      ) {
        stop(
          "Ground-truth column '",
          ground_truth_col,
          "' not found for: ",
          method
        )
      }


      # Spatial coordinates


      coords <- as.data.frame(
        SpatialExperiment::spatialCoords(spe)
      )

      if (
        ncol(coords) < 2 ||
          nrow(coords) != ncol(spe)
      ) {
        stop(
          "Invalid spatial coordinates for: ",
          method,
          " / ",
          slice_name
        )
      }


      # Construct plotting data


      predicted <- if (
        cluster_col %in% colnames(cd)
      ) {
        as.character(cd[[cluster_col]])
      } else {
        rep(NA_character_, ncol(spe))
      }

      truth <- if (
        !is.null(ground_truth_col) &&
          ground_truth_col %in% colnames(cd)
      ) {
        as.character(cd[[ground_truth_col]])
      } else {
        rep(NA_character_, ncol(spe))
      }

      df_list[[
        paste(method, slice_name, sep = "_")
      ]] <- data.frame(
        x = coords[, 1],
        y = coords[, 2],
        method = method,
        slice = slice_name,
        predicted = predicted,
        truth = truth,
        stringsAsFactors = FALSE
      )
    }
  }


  # Combine data


  df <- dplyr::bind_rows(df_list)

  df$method <- factor(
    df$method,
    levels = methods
  )

  df$slice <- factor(
    df$slice,
    levels = unique(df$slice)
  )


  # Resolve colour palettes


  make_palette <- function(values, supplied) {
    levels_use <- sort(
      unique(
        values[
          !is.na(values) & nzchar(values)
        ]
      )
    )

    n_levels <- length(levels_use)

    if (n_levels == 0) {
      return(character(0))
    }

    if (is.null(supplied)) {
      if (n_levels > 36) {
        stop(
          "More than 36 categories detected. ",
          "Please supply a custom palette."
        )
      }

      supplied <- Polychrome::palette36.colors(
        n_levels
      )

      names(supplied) <- levels_use
    } else if (is.null(names(supplied))) {
      if (length(supplied) < n_levels) {
        stop("Not enough colours supplied.")
      }

      supplied <- supplied[seq_len(n_levels)]

      names(supplied) <- levels_use
    } else {
      missing <- setdiff(
        levels_use,
        names(supplied)
      )

      if (length(missing) > 0) {
        stop(
          "Missing palette colours for: ",
          paste(missing, collapse = ", ")
        )
      }

      supplied <- supplied[levels_use]
    }

    supplied
  }


  # Predicted cluster colours


  if (show %in% c("both", "predicted")) {
    palette <- make_palette(
      df$predicted,
      palette
    )
  }


  # Ground-truth colours


  if (show %in% c("both", "truth")) {
    truth_palette <- make_palette(
      df$truth,
      truth_palette
    )
  }


  # Common spatial plotting helper


  make_spatial_plot <- function(data,
                                annotation_col,
                                colours,
                                legend_title) {
    data$annotation <- factor(
      data[[annotation_col]],
      levels = names(colours)
    )

    ggplot2::ggplot(
      data,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        colour = .data$annotation
      )
    ) +
      ggplot2::geom_point(
        size = point_size
      ) +
      ggplot2::coord_fixed() +
      ggplot2::scale_colour_manual(
        values = colours,
        na.value = "grey85",
        drop = FALSE
      ) +
      ggplot2::labs(
        x = "x",
        y = "y",
        colour = legend_title
      ) +
      .get_plot_theme(theme_name) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(
          size = 9,
          face = "bold"
        ),
        axis.text = ggplot2::element_text(
          size = 8
        ),
        legend.position = "right"
      ) +
      ggplot2::guides(
        colour = ggplot2::guide_legend(
          override.aes = list(
            size = legend_point_size
          )
        )
      )
  }


  # Ground truth: one plot per slice


  p_truth <- NULL

  if (show %in% c("both", "truth")) {
    # Use the first selected method as the reference
    # for each slice. Ground truth is not repeated
    # across normalisation methods.

    truth_df <- df[
      as.character(df$method) == methods[1], ,
      drop = FALSE
    ]

    if (nrow(truth_df) == 0) {
      stop("Ground-truth data could not be extracted.")
    }

    p_truth <- make_spatial_plot(
      data = truth_df,
      annotation_col = "truth",
      colours = truth_palette,
      legend_title = "Ground truth"
    )

    if (length(unique(truth_df$slice)) > 1) {
      p_truth <- p_truth +
        ggplot2::facet_wrap(
          ~slice,
          ncol = ncol
        )
    }

    p_truth <- p_truth +
      ggplot2::labs(
        title = "Ground truth"
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold"
        )
      )
  }


  # Predicted clusters: all selected methods


  p_predicted <- NULL

  if (show %in% c("both", "predicted")) {
    p_predicted <- make_spatial_plot(
      data = df,
      annotation_col = "predicted",
      colours = palette,
      legend_title = "BANKSY cluster"
    )

    if (length(unique(df$slice)) > 1) {
      p_predicted <- p_predicted +
        ggplot2::facet_grid(
          method ~ slice
        )
    } else {
      p_predicted <- p_predicted +
        ggplot2::facet_wrap(
          ~method,
          ncol = ncol
        )
    }

    p_predicted <- p_predicted +
      ggplot2::labs(
        title = "BANKSY clustering"
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold"
        )
      )
  }


  # Combine ground truth and predicted clusters


  if (show == "both") {
    combined_plot <- patchwork::wrap_plots(
      p_truth,
      p_predicted,
      ncol = 1,
      heights = c(1, 2)
    )
  } else if (show == "truth") {
    combined_plot <- p_truth
  } else {
    combined_plot <- p_predicted
  }


  # Return


  return(
    list(
      data = df,
      palette = palette,
      truth_palette = truth_palette,
      plots = list(
        truth = p_truth,
        predicted = p_predicted
      ),
      plot = combined_plot
    )
  )
}


# BANKSY clustering


#' Run BANKSY clustering independently for each slice
#'
#' Run BANKSY clustering independently across normalisation
#' methods and spatial slices.
#'
#' BANKSY cluster assignments are added to the colData of each
#' SpatialExperiment object so that the same objects can be used
#' in downstream clustering evaluation, plotting, differential
#' expression, and GO enrichment analyses.
#'
#' The function supports either:
#' \itemize{
#'   \item a named list of SpatialExperiment objects, or
#'   \item a named list where each method contains a named list
#'   of SpatialExperiment objects corresponding to multiple slices.
#' }
#'
#' @param norm_list Named list containing normalised
#'   SpatialExperiment objects.
#' @param methods Optional character vector specifying the
#'   normalisation methods to analyse.
#' @param assay_name Name of the assay used by BANKSY.
#' @param cluster_col Name of the colData column used to store
#'   BANKSY cluster assignments.
#' @param k_geom Number of spatial neighbours used by BANKSY.
#' @param lambda BANKSY spatial weighting parameter.
#' @param resolution Resolution parameter used for clustering.
#' @param seed Random seed for reproducibility.
#'
#' @return Named list with the same structure as `norm_list`.
#'   BANKSY cluster assignments are added to the colData of
#'   each SpatialExperiment object.
#'
#' @examples
#' data(norm_example)
#'
#' banksy_norm <- run_banksy_individual(
#'   norm_list = norm_example,
#'   methods = names(norm_example),
#'   assay_name = "logcounts",
#'   cluster_col = "banksy_cluster",
#'   k_geom = 15,
#'   lambda = 0.8,
#'   resolution = 0.2,
#'   seed = 123
#' )
#' @export
run_banksy_individual <- function(
  norm_list,
  methods = names(norm_list),
  assay_name = "logcounts",
  cluster_col = "banksy_cluster",
  k_geom = 15,
  lambda = 0.8,
  resolution = 0.2,
  seed = 123
) {
  # Check input


  .check_named_list(
    norm_list
  )


  missing_methods <- setdiff(
    methods,
    names(norm_list)
  )


  if (length(missing_methods) > 0) {
    stop(
      "The following methods are not present in 'norm_list': ",
      paste(
        missing_methods,
        collapse = ", "
      )
    )
  }


  # Copy input


  out <- norm_list


  # Run BANKSY


  for (method in methods) {
    method_data <- norm_list[[method]]


    # Single SpatialExperiment


    if (
      methods::is(
        method_data,
        "SpatialExperiment"
      )
    ) {
      spe_list <- list(
        slice1 = method_data
      )

      single_spe <- TRUE


      # Multiple SpatialExperiment slices
    } else if (is.list(method_data)) {
      spe_list <- method_data

      single_spe <- FALSE


      if (
        is.null(names(spe_list)) ||
          any(names(spe_list) == "")
      ) {
        stop(
          "Slice lists must be named for method: ",
          method
        )
      }
    } else {
      stop(
        "Each normalisation method must contain either a ",
        "SpatialExperiment object or a named list of ",
        "SpatialExperiment objects."
      )
    }


    # Process slices


    for (slice_name in names(spe_list)) {
      spe <- spe_list[[slice_name]]


      # Check SPE


      if (
        !methods::is(
          spe,
          "SpatialExperiment"
        )
      ) {
        stop(
          "Object for method '",
          method,
          "', slice '",
          slice_name,
          "' is not a SpatialExperiment."
        )
      }


      # Check assay


      if (
        !assay_name %in%
          SummarizedExperiment::assayNames(spe)
      ) {
        stop(
          "Assay '",
          assay_name,
          "' not found for method '",
          method,
          "', slice '",
          slice_name,
          "'."
        )
      }


      # Compute BANKSY


      spe <- Banksy::computeBanksy(
        spe,
        assay_name = assay_name,
        compute_agf = FALSE,
        k_geom = k_geom
      )


      # Run BANKSY PCA


      spe <- Banksy::runBanksyPCA(
        spe,
        lambda = lambda,
        use_agf = FALSE,
        seed = seed
      )


      # Run BANKSY clustering


      spe <- Banksy::clusterBanksy(
        spe,
        lambda = lambda,
        use_agf = FALSE,
        resolution = resolution,
        seed = seed
      )


      # Identify BANKSY cluster column


      cluster_name <- grep(
        "^clust_",
        colnames(
          SummarizedExperiment::colData(spe)
        ),
        value = TRUE
      )


      if (length(cluster_name) == 0) {
        stop(
          "BANKSY cluster column could not be identified ",
          "for method '",
          method,
          "', slice '",
          slice_name,
          "'."
        )
      }


      cluster_name <- utils::tail(
        cluster_name,
        1
      )


      # Store standard cluster column


      SummarizedExperiment::colData(spe)[[cluster_col]] <-
        as.factor(
          SummarizedExperiment::colData(spe)[[
            cluster_name
          ]]
        )


      # Store updated SPE


      spe_list[[slice_name]] <- spe
    }


    # Restore original structure


    if (single_spe) {
      out[[method]] <- spe_list[[1]]
    } else {
      out[[method]] <- spe_list
    }
  }


  # Return


  out
}


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
#'
#' @examples
#' data(norm_example)
#'
#' banksy_norm <- run_banksy_individual(
#'   norm_list = norm_example,
#'   seed = 123
#' )
#'
#' metrics_df <- calculate_clustering_metrics(
#'   norm_list = banksy_norm,
#'   ground_truth_col = "layer_guess",
#'   cluster_col = "banksy_cluster"
#' )
#'
#' metrics_df
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
        stop(
          "Missing colData columns for ", method, " / ", slice_name,
          ": ", paste(missing_cols, collapse = ", ")
        )
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
          call. = FALSE
        )
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
        withr::local_seed(seed)
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
          method, " / ", slice_name,
          call. = FALSE
        )
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
#' @examples
#' data(norm_example)
#'
#' banksy_norm <- run_banksy_individual(
#'   norm_list = norm_example,
#'   seed = 123
#' )
#'
#' metrics_df <- calculate_clustering_metrics(
#'   norm_list = banksy_norm,
#'   ground_truth_col = "layer_guess",
#'   cluster_col = "banksy_cluster"
#' )
#'
#' p <- plot_clustering_metrics(
#'   metrics_df = metrics_df,
#'   metrics = c("ARI", "NMI", "HOM", "COM", "ASW", "Overall"),
#'   point_size = 4,
#'   ncol = 2,
#'   scales = "free_y",
#'   theme_name = "bw"
#' )
#'
#' p$plot
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
      x = NULL, # "Normalisation method",
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
#' @examples
#' data(norm_example)
#'
#' banksy_norm <- run_banksy_individual(
#'   norm_list = norm_example,
#'   seed = 123
#' )
#'
#' metrics_df <- calculate_clustering_metrics(
#'   norm_list = banksy_norm,
#'   ground_truth_col = "layer_guess"
#' )
#'
#' ranking_df <- rank_normalization_methods(
#'   metrics_df = metrics_df
#' )
#'
#' ranking_df
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
#' @examples
#' data(norm_example)
#'
#' banksy_norm <- run_banksy_individual(
#'   norm_list = norm_example,
#'   seed = 123
#' )
#'
#' metrics_df <- calculate_clustering_metrics(
#'   norm_list = banksy_norm,
#'   ground_truth_col = "layer_guess"
#' )
#'
#' ranking_df <- rank_normalization_methods(
#'   metrics_df = metrics_df
#' )
#'
#' p <- plot_normalization_ranking(
#'   ranking_df = ranking_df,
#'   bar_width = 0.7,
#'   label_size = 3.5,
#'   theme_name = "classic"
#' )
#'
#' p$plot
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
    order(-ranking_df$Overall), ,
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
      y = NULL, # "Normalisation method",
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



#' Run BANKSY with different gene-selection strategies
#'
#' Runs BANKSY clustering independently for each normalisation
#' method using all genes, highly variable genes (HVGs), or
#' spatially variable genes (SVGs).
#'
#' HVGs are selected by expression variance within each
#' normalisation method. SVGs are selected using the gene
#' rankings supplied in svg_list.
#'
#' @param norm_list Named list of SpatialExperiment objects,
#'   or named lists of SpatialExperiment slices.
#' @param svg_list Optional named list of SVG results,
#'   corresponding to the normalisation methods.
#' @param methods Normalisation methods to analyse.
#' @param gene_selection Character vector containing any of
#'   "all", "hvg500", "hvg1000", "svg100", "svg200", "svg500".
#' @param rank_col Column containing SVG rankings.
#' @param assay_name Assay used for clustering and HVG selection.
#' @param cluster_col Column for BANKSY cluster assignments.
#' @param k_geom Number of spatial neighbours.
#' @param lambda BANKSY spatial weighting parameter.
#' @param resolution Clustering resolution.
#' @param seed Random seed.
#'
#' @return Nested list indexed by normalisation method,
#'   then gene-selection strategy. Each entry contains a
#'   clustered SpatialExperiment object or named slice list.
#'
#' @examples
#' data(norm_example)
#'
#' banksy_genes <- run_banksy_gene_selection(
#'   norm_list = norm_example,
#'   methods = names(norm_example),
#'   gene_selection = "all",
#'   assay_name = "logcounts",
#'   k_geom = 15,
#'   lambda = 0.8,
#'   resolution = 0.2,
#'   seed = 123
#' )
#' @export
run_banksy_gene_selection <- function(
  norm_list,
  svg_list = NULL,
  methods = names(norm_list),
  gene_selection = c(
    "all", "hvg500", "hvg1000",
    "svg100", "svg200", "svg500"
  ),
  rank_col = "rank",
  assay_name = "logcounts",
  cluster_col = "banksy_cluster",
  k_geom = 15,
  lambda = 0.8,
  resolution = 0.2,
  seed = 123
) {
  # Validate normalisation methods

  .check_named_list(norm_list)

  missing_methods <- setdiff(
    methods,
    names(norm_list)
  )

  if (length(missing_methods) > 0L) {
    stop(
      "Unknown normalisation methods: ",
      paste(missing_methods, collapse = ", ")
    )
  }

  methods <- unique(methods)

  if (length(methods) == 0L) {
    stop("No normalisation methods selected.")
  }

  # Validate gene-selection strategies

  allowed <- c(
    "all", "hvg500", "hvg1000",
    "svg100", "svg200", "svg500"
  )

  gene_selection <- unique(gene_selection)

  unknown <- setdiff(gene_selection, allowed)

  if (length(unknown) > 0L) {
    stop(
      "Unknown gene-selection strategies: ",
      paste(unknown, collapse = ", ")
    )
  }

  if (length(gene_selection) == 0L) {
    stop("No gene-selection strategies selected.")
  }

  needs_svg <- any(
    grepl("^svg", gene_selection)
  )

  if (needs_svg && is.null(svg_list)) {
    stop(
      "svg_list is required for SVG gene selection."
    )
  }

  if (needs_svg) {
    missing_svg <- setdiff(methods, names(svg_list))

    if (length(missing_svg) > 0L) {
      stop(
        "Missing SVG results for: ",
        paste(missing_svg, collapse = ", ")
      )
    }
  }

  # Helper: select genes from a single slice

  select_genes <- function(
    spe,
    strategy,
    svg_spe = NULL
  ) {
    genes <- rownames(spe)

    if (is.null(genes) ||
      anyDuplicated(genes) ||
      anyNA(genes)) {
      stop("Valid unique gene names are required.")
    }

    if (strategy == "all") {
      return(genes)
    }

    n_select <- as.integer(
      sub("^(hvg|svg)", "", strategy)
    )

    # Highly variable genes

    if (startsWith(strategy, "hvg")) {
      expr <- SummarizedExperiment::assay(
        spe,
        assay_name
      )

      variances <- MatrixGenerics::rowVars(
        expr,
        na.rm = TRUE
      )

      names(variances) <- genes

      variances[!is.finite(variances)] <- -Inf

      ordered <- order(
        variances,
        decreasing = TRUE
      )

      selected <- genes[
        utils::head(ordered, min(n_select, length(genes)))
      ]

      return(selected)
    }

    # Spatially variable genes

    if (is.null(svg_spe)) {
      stop("SVG results are missing.")
    }

    svg_data <- as.data.frame(
      SummarizedExperiment::rowData(svg_spe)
    )

    if (!rank_col %in% colnames(svg_data)) {
      stop(
        "SVG ranking column '",
        rank_col,
        "' not found."
      )
    }

    svg_genes <- rownames(svg_spe)

    ranks <- svg_data[[rank_col]]

    if (!is.numeric(ranks)) {
      stop("SVG rankings must be numeric.")
    }

    keep <- svg_genes %in% genes &
      is.finite(ranks)

    svg_genes <- svg_genes[keep]
    ranks <- ranks[keep]

    if (length(svg_genes) == 0L) {
      stop("No ranked SVGs overlap the expression data.")
    }

    ordered <- order(
      ranks,
      decreasing = FALSE
    )

    utils::head(
      svg_genes[ordered],
      min(n_select, length(svg_genes))
    )
  }

  # Run BANKSY for each method

  out <- list()

  for (method in methods) {
    method_data <- norm_list[[method]]

    single_spe <- methods::is(
      method_data,
      "SpatialExperiment"
    )

    if (single_spe) {
      slices <- list(
        slice1 = method_data
      )
    } else if (is.list(method_data)) {
      slices <- method_data

      if (
        length(slices) == 0L ||
          is.null(names(slices)) ||
          anyNA(names(slices)) ||
          any(!nzchar(names(slices))) ||
          anyDuplicated(names(slices))
      ) {
        stop(
          "Slice lists must have unique names: ",
          method
        )
      }
    } else {
      stop(
        "Invalid normalisation input: ",
        method
      )
    }

    out[[method]] <- list()

    for (strategy in gene_selection) {
      clustered_slices <- list()

      for (slice_name in names(slices)) {
        spe <- slices[[slice_name]]

        if (!methods::is(
          spe,
          "SpatialExperiment"
        )) {
          stop(
            "Invalid SpatialExperiment: ",
            method, " / ", slice_name
          )
        }

        if (!assay_name %in%
          SummarizedExperiment::assayNames(spe)) {
          stop(
            "Assay '", assay_name,
            "' missing for ", method
          )
        }

        # Resolve corresponding SVG results

        svg_spe <- NULL

        if (startsWith(strategy, "svg")) {
          svg_method <- svg_list[[method]]

          if (methods::is(
            svg_method,
            "SpatialExperiment"
          )) {
            svg_spe <- svg_method
          } else if (is.list(svg_method)) {
            svg_spe <- svg_method[[slice_name]]
          }

          if (is.null(svg_spe)) {
            stop(
              "Missing SVG results for ",
              method, " / ", slice_name
            )
          }
        }

        # Select genes

        selected_genes <- select_genes(
          spe = spe,
          strategy = strategy,
          svg_spe = svg_spe
        )

        if (length(selected_genes) < 2L) {
          stop(
            "Too few genes selected for ",
            method, " / ", strategy
          )
        }

        spe_sub <- spe[
          selected_genes, ,
          drop = FALSE
        ]

        # Run BANKSY using existing package function

        result <- run_banksy_individual(
          norm_list = list(
            selected = spe_sub
          ),
          methods = "selected",
          assay_name = assay_name,
          cluster_col = cluster_col,
          k_geom = k_geom,
          lambda = lambda,
          resolution = resolution,
          seed = seed
        )

        clustered_slices[[slice_name]] <-
          result[["selected"]]
      }

      # Preserve single- or multi-slice structure

      if (single_spe) {
        out[[method]][[strategy]] <-
          clustered_slices[[1L]]
      } else {
        out[[method]][[strategy]] <-
          clustered_slices
      }
    }
  }

  out
}


#' Calculate clustering metrics across gene-selection strategies
#'
#' Calculate ARI, NMI, HOM, COM, spatial ASW, and Overall
#' clustering agreement for each gene-selection strategy.
#'
#' @param banksy_list Output from run_banksy_gene_selection().
#' @param ground_truth_col Ground-truth annotation column.
#' @param cluster_col Predicted cluster column.
#' @param asw_max_cells Maximum cells used for spatial ASW.
#' @param seed Random seed.
#'
#' @return Data frame containing clustering metrics and gene sets.
#' #' @examples
#' data(norm_example)
#'
#' banksy_genes <- run_banksy_gene_selection(
#'   norm_list = norm_example,
#'   methods = names(norm_example),
#'   gene_selection = "all",
#'   assay_name = "logcounts",
#'   k_geom = 15,
#'   lambda = 0.8,
#'   resolution = 0.2,
#'   seed = 123
#' )
#'
#' metrics_genes <- calculate_gene_selection_metrics(
#'   banksy_list = banksy_genes,
#'   ground_truth_col = "layer_guess"
#' )
#'
#' @examples
#' data(norm_example)
#'
#' banksy_results <- run_banksy_individual(
#'   norm_list = norm_example,
#'   methods = names(norm_example),
#'   assay_name = "logcounts",
#'   cluster_col = "banksy_cluster",
#'   k_geom = 15,
#'   lambda = 0.8,
#'   resolution = 0.2,
#'   seed = 123
#' )
#'
#' banksy_list <- list(all = banksy_results)
#'
#' metrics_genes <- calculate_gene_selection_metrics(
#'   banksy_list = banksy_list,
#'   ground_truth_col = "layer_guess",
#'   cluster_col = "banksy_cluster",
#'   asw_max_cells = 3000,
#'   seed = 123
#' )
#'
#' head(metrics_genes)
#' @export
calculate_gene_selection_metrics <- function(
  banksy_list,
  ground_truth_col = "layer_guess",
  cluster_col = "banksy_cluster",
  asw_max_cells = 3000,
  seed = 123
) {
  results <- lapply(
    names(banksy_list),
    function(strategy) {
      df <- calculate_clustering_metrics(
        norm_list = banksy_list[[strategy]],
        ground_truth_col = ground_truth_col,
        cluster_col = cluster_col,
        asw_max_cells = asw_max_cells,
        seed = seed
      )

      df$gene_selection <- strategy

      df
    }
  )

  out <- do.call(rbind, results)

  rownames(out) <- NULL

  return(out)
}


#' Plot BANKSY gene-selection clustering metrics
#'
#' Visualise BANKSY clustering performance across
#' normalisation methods and gene-selection strategies.
#'
#' Points are coloured by normalisation method and
#' distinguished by shape according to gene selection.
#'
#' Metrics are displayed in separate facets.
#'
#' @param metrics_df Data frame returned by
#'   calculate_gene_selection_metrics().
#' @param metrics Character vector of clustering metrics
#'   to display.
#' @param point_size Point size.
#' @param point_alpha Point opacity.
#' @param ncol Number of facet columns.
#' @param scales Facet scales. One of "fixed", "free",
#'   "free_x", or "free_y".
#' @param palette Optional named colour vector for
#'   normalisation methods.
#' @param shapes Optional named numeric vector specifying
#'   plotting symbols for gene-selection strategies.
#' @param method_order Optional ordering of normalisation
#'   methods.
#' @param gene_selection_order Optional ordering of
#'   gene-selection strategies.
#' @param theme_name Plot theme. One of "classic", "bw",
#'   "minimal", or "gray".
#' @param rotate_x Angle of x-axis labels.
#'
#' @return A list containing the plotting data and
#'   the combined ggplot object.
#'
#' @examples
#' \dontshow{
#' # Examples require clustering metrics.
#' }
#'
#' @examples
#' data(norm_example)
#'
#' banksy_genes <- run_banksy_gene_selection(
#'   norm_list = norm_example,
#'   methods = names(norm_example),
#'   gene_selection = "all",
#'   assay_name = "logcounts",
#'   k_geom = 15,
#'   lambda = 0.8,
#'   resolution = 0.2,
#'   seed = 123
#' )
#'
#' metrics_genes <- calculate_gene_selection_metrics(
#'   banksy_list = banksy_genes,
#'   ground_truth_col = "layer_guess"
#' )
#'
#' p_genes <- plot_gene_selection_metrics(
#'   metrics_df = metrics_genes,
#'   point_size = 4,
#'   ncol = 3,
#'   scales = "free_y",
#'   theme_name = "classic"
#' )
#'
#' p_genes$plot
#' @export
plot_gene_selection_metrics <- function(
  metrics_df,
  metrics = c(
    "ARI", "NMI", "HOM",
    "COM", "ASW", "Overall"
  ),
  point_size = 4,
  point_alpha = 0.9,
  ncol = 3,
  scales = "free_y",
  palette = NULL,
  shapes = NULL,
  method_order = NULL,
  gene_selection_order = NULL,
  theme_name = "classic",
  rotate_x = 45
) {
  # Validate input

  if (
    !is.data.frame(metrics_df) ||
      nrow(metrics_df) == 0L
  ) {
    stop(
      "metrics_df must be a non-empty data frame."
    )
  }

  required_cols <- c(
    "method",
    "gene_selection",
    metrics
  )

  missing_cols <- setdiff(
    required_cols,
    colnames(metrics_df)
  )

  if (length(missing_cols) > 0L) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (length(metrics) == 0L ||
    anyDuplicated(metrics)) {
    stop(
      "metrics must contain unique metric names."
    )
  }

  if (!all(vapply(
    metrics_df[metrics],
    is.numeric,
    logical(1)
  ))) {
    stop("All selected metrics must be numeric.")
  }

  scales <- match.arg(
    scales,
    c("fixed", "free", "free_x", "free_y")
  )

  # Resolve method order

  available_methods <- unique(
    as.character(metrics_df$method)
  )

  if (is.null(method_order)) {
    method_order <- available_methods
  } else {
    method_order <- c(
      intersect(method_order, available_methods),
      setdiff(available_methods, method_order)
    )
  }

  # Resolve gene-selection order

  available_strategies <- unique(
    as.character(metrics_df$gene_selection)
  )

  if (is.null(gene_selection_order)) {
    gene_selection_order <- available_strategies
  } else {
    gene_selection_order <- c(
      intersect(
        gene_selection_order,
        available_strategies
      ),
      setdiff(
        available_strategies,
        gene_selection_order
      )
    )
  }

  # Prepare plotting data

  df_list <- lapply(
    metrics,
    function(metric) {
      data.frame(
        method = as.character(
          metrics_df$method
        ),
        gene_selection = as.character(
          metrics_df$gene_selection
        ),
        metric = metric,
        score = as.numeric(
          metrics_df[[metric]]
        ),
        stringsAsFactors = FALSE
      )
    }
  )

  df <- dplyr::bind_rows(df_list)

  df$method <- factor(
    df$method,
    levels = method_order
  )

  df$gene_selection <- factor(
    df$gene_selection,
    levels = gene_selection_order
  )

  df$metric <- factor(
    df$metric,
    levels = metrics
  )

  # Resolve shapes

  if (is.null(shapes)) {
    default_shapes <- c(
      16, 17, 15, 18, 3, 7, 8, 4
    )

    if (
      length(gene_selection_order) >
        length(default_shapes)
    ) {
      stop(
        "Too many gene-selection strategies ",
        "for the default shape palette. ",
        "Please supply 'shapes'."
      )
    }

    shapes <- stats::setNames(
      default_shapes[
        seq_along(gene_selection_order)
      ],
      gene_selection_order
    )
  } else {
    if (
      is.null(names(shapes)) ||
        !all(
          gene_selection_order %in% names(shapes)
        )
    ) {
      stop(
        "shapes must be a named vector ",
        "containing all gene-selection strategies."
      )
    }

    shapes <- shapes[
      gene_selection_order
    ]
  }

  # Construct plot

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$method,
      y = .data$score,
      colour = .data$method,
      shape = .data$gene_selection
    )
  ) +
    ggplot2::geom_point(
      size = point_size,
      alpha = point_alpha,
      na.rm = TRUE
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$metric),
      ncol = ncol,
      scales = scales
    ) +
    ggplot2::scale_shape_manual(
      values = shapes,
      name = "Gene selection",
      drop = FALSE
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Clustering metric score",
      colour = "Normalisation method"
    ) +
    .get_plot_theme(theme_name) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = rotate_x,
        hjust = 1
      ),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      legend.position = "right"
    )

  # Optional method colours

  if (!is.null(palette)) {
    if (
      is.null(names(palette)) ||
        !all(
          method_order %in% names(palette)
        )
    ) {
      stop(
        "palette must be a named colour vector ",
        "containing all normalisation methods."
      )
    }

    p <- p +
      ggplot2::scale_colour_manual(
        values = palette,
        breaks = method_order
      )
  }

  # Return

  list(
    data = df,
    plot = p
  )
}
