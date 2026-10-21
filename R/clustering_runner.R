# ============================================================
# BANKSY clustering
# ============================================================


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

  # ----------------------------------------------------------
  # Check input
  # ----------------------------------------------------------

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


  # ----------------------------------------------------------
  # Copy input
  # ----------------------------------------------------------

  out <- norm_list


  # ----------------------------------------------------------
  # Run BANKSY
  # ----------------------------------------------------------

  for (method in methods) {

    method_data <- norm_list[[method]]


    # --------------------------------------------------------
    # Single SpatialExperiment
    # --------------------------------------------------------

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


      # --------------------------------------------------------
      # Multiple SpatialExperiment slices
      # --------------------------------------------------------

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


    # --------------------------------------------------------
    # Process slices
    # --------------------------------------------------------

    for (slice_name in names(spe_list)) {

      spe <- spe_list[[slice_name]]


      # ------------------------------------------------------
      # Check SPE
      # ------------------------------------------------------

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


      # ------------------------------------------------------
      # Check assay
      # ------------------------------------------------------

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


      # ------------------------------------------------------
      # Compute BANKSY
      # ------------------------------------------------------

      spe <- Banksy::computeBanksy(
        spe,
        assay_name = assay_name,
        compute_agf = FALSE,
        k_geom = k_geom
      )


      # ------------------------------------------------------
      # Run BANKSY PCA
      # ------------------------------------------------------

      spe <- Banksy::runBanksyPCA(
        spe,
        lambda = lambda,
        use_agf = FALSE,
        seed = seed
      )


      # ------------------------------------------------------
      # Run BANKSY clustering
      # ------------------------------------------------------

      spe <- Banksy::clusterBanksy(
        spe,
        lambda = lambda,
        use_agf = FALSE,
        resolution = resolution,
        seed = seed
      )


      # ------------------------------------------------------
      # Identify BANKSY cluster column
      # ------------------------------------------------------

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


      # ------------------------------------------------------
      # Store standard cluster column
      # ------------------------------------------------------

      SummarizedExperiment::colData(spe)[[cluster_col]] <-
        as.factor(
          SummarizedExperiment::colData(spe)[[
            cluster_name
          ]]
        )


      # ------------------------------------------------------
      # Store updated SPE
      # ------------------------------------------------------

      spe_list[[slice_name]] <- spe
    }


    # --------------------------------------------------------
    # Restore original structure
    # --------------------------------------------------------

    if (single_spe) {

      out[[method]] <- spe_list[[1]]

    } else {

      out[[method]] <- spe_list
    }
  }


  # ----------------------------------------------------------
  # Return
  # ----------------------------------------------------------

  out
}



# ============================================================
# Adjusted Rand index
# ============================================================


#' Calculate adjusted Rand index
#'
#' Calculate the adjusted Rand index between clustering
#' assignments and user-specified ground-truth labels across
#' normalisation methods and spatial slices.
#'
#' @param norm_list Named list of SpatialExperiment objects
#'   containing clustering results.
#' @param ground_truth_col Name of the colData column containing
#'   ground-truth labels.
#' @param cluster_col Name of the colData column containing
#'   cluster assignments.
#' @param methods Optional character vector specifying the
#'   normalisation methods to analyse.
#'
#' @return A data frame containing ARI values for each
#'   normalisation method and spatial slice.
#'
#' @export
calculate_ari <- function(
    norm_list,
    ground_truth_col,
    cluster_col = "banksy_cluster",
    methods = names(norm_list)
) {

  # ----------------------------------------------------------
  # Check input
  # ----------------------------------------------------------

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


  # ----------------------------------------------------------
  # Initialise output
  # ----------------------------------------------------------

  out_list <- list()


  # ----------------------------------------------------------
  # Calculate ARI
  # ----------------------------------------------------------

  for (method in methods) {

    method_data <- norm_list[[method]]


    # --------------------------------------------------------
    # Single SPE
    # --------------------------------------------------------

    if (
      methods::is(
        method_data,
        "SpatialExperiment"
      )
    ) {

      spe_list <- list(
        slice1 = method_data
      )


      # --------------------------------------------------------
      # Multiple slices
      # --------------------------------------------------------

    } else if (is.list(method_data)) {

      spe_list <- method_data


    } else {

      stop(
        "Each method must contain either a SpatialExperiment ",
        "or a list of SpatialExperiment objects."
      )
    }


    # --------------------------------------------------------
    # Process slices
    # --------------------------------------------------------

    for (slice_name in names(spe_list)) {

      spe <- spe_list[[slice_name]]


      cd <- SummarizedExperiment::colData(
        spe
      )


      # ------------------------------------------------------
      # Check columns
      # ------------------------------------------------------

      required_cols <- c(
        ground_truth_col,
        cluster_col
      )


      missing_cols <- setdiff(
        required_cols,
        colnames(cd)
      )


      if (length(missing_cols) > 0) {

        stop(
          "Missing colData column(s) for method '",
          method,
          "', slice '",
          slice_name,
          "': ",
          paste(
            missing_cols,
            collapse = ", "
          )
        )
      }


      # ------------------------------------------------------
      # Extract labels
      # ------------------------------------------------------

      cluster <- cd[[cluster_col]]

      truth <- cd[[ground_truth_col]]


      # ------------------------------------------------------
      # Remove missing values
      # ------------------------------------------------------

      keep <- !is.na(cluster) &
        !is.na(truth)


      cluster <- cluster[
        keep
      ]

      truth <- truth[
        keep
      ]


      # ------------------------------------------------------
      # Check data
      # ------------------------------------------------------

      if (length(cluster) == 0) {

        warning(
          "No complete observations for method '",
          method,
          "', slice '",
          slice_name,
          "'.",
          call. = FALSE
        )

        next
      }


      # ------------------------------------------------------
      # Calculate ARI
      # ------------------------------------------------------

      ari_value <- mclust::adjustedRandIndex(
        cluster,
        truth
      )


      # ------------------------------------------------------
      # Store result
      # ------------------------------------------------------

      out_list[[
        paste(
          method,
          slice_name,
          sep = "_"
        )
      ]] <- data.frame(

        method = method,

        slice = slice_name,

        ARI = ari_value,

        n = length(cluster),

        stringsAsFactors = FALSE
      )
    }
  }


  # ----------------------------------------------------------
  # Check output
  # ----------------------------------------------------------

  if (length(out_list) == 0) {

    stop(
      "No ARI values could be calculated."
    )
  }


  # ----------------------------------------------------------
  # Combine
  # ----------------------------------------------------------

  out <- do.call(
    rbind,
    out_list
  )


  rownames(out) <- NULL


  # ----------------------------------------------------------
  # Return
  # ----------------------------------------------------------

  out
}



# ============================================================
# ARI dot plot
# ============================================================


#' Plot adjusted Rand index
#'
#' Create a dot plot of adjusted Rand index values across
#' normalisation methods.
#'
#' @param ari_df Data frame returned by `calculate_ari()`.
#' @param method_col Name of the column containing normalisation
#'   methods.
#' @param ari_col Name of the column containing ARI values.
#' @param method_order Optional character vector specifying the
#'   order of normalisation methods.
#' @param point_size Point size.
#' @param palette Optional named character vector specifying
#'   colours for normalisation methods.
#' @param theme_name Plot theme. One of `"classic"`, `"bw"`,
#'   `"minimal"`, or `"gray"`.
#'
#' @return A list containing plotting data and plot.
#'
#' @export
plot_ari <- function(
    ari_df,
    method_col = "method",
    ari_col = "ARI",
    method_order = NULL,
    point_size = 4,
    palette = NULL,
    theme_name = "classic"
) {

  # ----------------------------------------------------------
  # Check input
  # ----------------------------------------------------------

  if (!is.data.frame(ari_df)) {

    stop(
      "'ari_df' must be a data frame."
    )
  }


  # ----------------------------------------------------------
  # Check required columns
  # ----------------------------------------------------------

  required_cols <- c(
    method_col,
    ari_col
  )


  missing_cols <- setdiff(
    required_cols,
    colnames(ari_df)
  )


  if (length(missing_cols) > 0) {

    stop(
      "The following required column(s) are missing: ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  }


  # ----------------------------------------------------------
  # Resolve methods
  # ----------------------------------------------------------

  available_methods <- unique(
    as.character(
      ari_df[[method_col]]
    )
  )


  if (is.null(method_order)) {

    methods <- available_methods

  } else {

    methods <- method_order[
      method_order %in%
        available_methods
    ]
  }


  if (length(methods) == 0) {

    stop(
      "No valid normalisation methods were found."
    )
  }


  # ----------------------------------------------------------
  # Prepare plotting data
  # ----------------------------------------------------------

  df <- ari_df[
    ari_df[[method_col]] %in%
      methods,
    ,
    drop = FALSE
  ]


  df[[method_col]] <- factor(
    df[[method_col]],
    levels = methods
  )


  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[method_col]],
      y = .data[[ari_col]],
      colour = .data[[method_col]]
    )
  ) +

    ggplot2::geom_point(
      size = point_size
    ) +

    ggplot2::labs(
      x = "Normalisation method",
      y = "Adjusted Rand index",
      colour = "Method"
    ) +

    .get_plot_theme(
      theme_name
    ) +

    ggplot2::theme(

      panel.grid =
        ggplot2::element_blank(),

      legend.position = "right",

      axis.text.x =
        ggplot2::element_text(
          angle = 45,
          hjust = 1
        )
    )


  # ----------------------------------------------------------
  # Optional colour palette
  # ----------------------------------------------------------

  if (!is.null(palette)) {

    p <- p +

      ggplot2::scale_colour_manual(
        values = palette
      )
  }


  # ----------------------------------------------------------
  # Return
  # ----------------------------------------------------------

  list(
    data = df,
    plot = p
  )
}





# ============================================================
# Spatial clustering plot
# ============================================================


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

  # ----------------------------------------------------------
  # Check input
  # ----------------------------------------------------------

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


  # ----------------------------------------------------------
  # Resolve methods
  # ----------------------------------------------------------

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


  # ----------------------------------------------------------
  # Extract spatial annotations
  # ----------------------------------------------------------

  df_list <- list()

  for (method in methods) {

    method_data <- norm_list[[method]]


    # --------------------------------------------------------
    # Single SPE or multiple slices
    # --------------------------------------------------------

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


    # --------------------------------------------------------
    # Process each slice
    # --------------------------------------------------------

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


      # ------------------------------------------------------
      # Check required columns
      # ------------------------------------------------------

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


      # ------------------------------------------------------
      # Spatial coordinates
      # ------------------------------------------------------

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


      # ------------------------------------------------------
      # Construct plotting data
      # ------------------------------------------------------

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


  # ----------------------------------------------------------
  # Combine data
  # ----------------------------------------------------------

  df <- dplyr::bind_rows(df_list)

  df$method <- factor(
    df$method,
    levels = methods
  )

  df$slice <- factor(
    df$slice,
    levels = unique(df$slice)
  )


  # ----------------------------------------------------------
  # Resolve colour palettes
  # ----------------------------------------------------------

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


  # ----------------------------------------------------------
  # Predicted cluster colours
  # ----------------------------------------------------------

  if (show %in% c("both", "predicted")) {

    palette <- make_palette(
      df$predicted,
      palette
    )

  }


  # ----------------------------------------------------------
  # Ground-truth colours
  # ----------------------------------------------------------

  if (show %in% c("both", "truth")) {

    truth_palette <- make_palette(
      df$truth,
      truth_palette
    )

  }


  # ----------------------------------------------------------
  # Common spatial plotting helper
  # ----------------------------------------------------------

  make_spatial_plot <- function(
    data,
    annotation_col,
    colours,
    legend_title
  ) {

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


  # ----------------------------------------------------------
  # Ground truth: one plot per slice
  # ----------------------------------------------------------

  p_truth <- NULL

  if (show %in% c("both", "truth")) {

    # Use the first selected method as the reference
    # for each slice. Ground truth is not repeated
    # across normalisation methods.

    truth_df <- df[
      as.character(df$method) == methods[1],
      ,
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
          ~ slice,
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


  # ----------------------------------------------------------
  # Predicted clusters: all selected methods
  # ----------------------------------------------------------

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
          ~ method,
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


  # ----------------------------------------------------------
  # Combine ground truth and predicted clusters
  # ----------------------------------------------------------

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


  # ----------------------------------------------------------
  # Return
  # ----------------------------------------------------------

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
