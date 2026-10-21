# ============================================================
# Cell-type annotation
# ============================================================


#' Run cell-type annotation
#'
#' Run cell-type annotation for a named list of normalised
#' SpatialExperiment objects using SingleR.
#'
#' Predicted labels are added to the colData of each
#' SpatialExperiment object so that the same objects can be used
#' in downstream plotting and benchmarking functions.
#'
#' @param norm_list Named list of normalised SpatialExperiment
#'   objects.
#' @param ref Reference SingleCellExperiment or
#'   SummarizedExperiment object.
#' @param labels Character vector containing reference cell-type
#'   labels. Length must equal the number of columns in `ref`.
#' @param method Cell-type annotation method.
#'   Currently only `"SingleR"` is supported.
#' @param label_name Name of the colData column used to store
#'   predicted labels.
#' @param pruned_label_name Name of the colData column used to
#'   store pruned SingleR labels.
#' @param ... Additional parameters passed to
#'   `SingleR::SingleR()`.
#'
#' @return Named list of annotated SpatialExperiment objects.
#'
#' @export
run_annotation <- function(
    norm_list,
    ref,
    labels,
    method = "SingleR",
    label_name = "SingleR_label",
    pruned_label_name = "SingleR_pruned_label",
    ...
) {

  # ----------------------------------------------------------
  # Check input
  # ----------------------------------------------------------

  .check_named_list(
    norm_list
  )


  if (method != "SingleR") {

    stop(
      "Currently supported annotation method: SingleR"
    )
  }


  if (length(labels) != ncol(ref)) {

    stop(
      "Length of 'labels' must equal the number of columns in 'ref'."
    )
  }


  # ----------------------------------------------------------
  # Run annotation
  # ----------------------------------------------------------

  out <- lapply(
    names(norm_list),
    function(method_name) {

      message("Running SingleR: ", method_name)

      spe <- norm_list[[method_name]]
      # ------------------------------------------------------
      # Run SingleR
      # ------------------------------------------------------

      results <- SingleR::SingleR(
        test = spe,
        ref = ref,
        labels = labels,
        ...
      )


      # ------------------------------------------------------
      # Check number of predictions
      # ------------------------------------------------------

      if (nrow(results) != ncol(spe)) {

        stop(
          "Number of SingleR predictions does not match the number of spatial locations."
        )
      }


      # ------------------------------------------------------
      # Add labels to colData
      # ------------------------------------------------------

      SummarizedExperiment::colData(spe)[[label_name]] <-
        as.character(
          results$labels
        )


      # ------------------------------------------------------
      # Add pruned labels
      # ------------------------------------------------------

      if ("pruned.labels" %in% colnames(results)) {

        SummarizedExperiment::colData(spe)[[pruned_label_name]] <-
          as.character(
            results$pruned.labels
          )
      }


      # ------------------------------------------------------
      # Return SpatialExperiment
      # ------------------------------------------------------

      spe
    }
  )


  names(out) <- names(
    norm_list
  )


  out
}



# ============================================================
# Spatial cell-type annotation plot
# ============================================================


#' Plot spatial cell-type annotations
#'
#' Plot cell-type annotations stored in the colData of
#' SpatialExperiment objects across normalisation methods.
#'
#' @param norm_list Named list of annotated SpatialExperiment
#'   objects.
#' @param annotation_col Name of the colData column containing
#'   cell-type annotations.
#' @param method_order Optional character vector specifying the
#'   order of normalisation methods.
#' @param point_size Point size for spatial locations.
#' @param palette Optional named character vector specifying
#'   colours for cell types.
#' @param ncol Number of facet columns.
#' @param theme_name Plot theme. One of `"classic"`, `"bw"`,
#'   `"minimal"`, or `"gray"`.
#'
#' @return A list containing plotting data and plot.
#'
#' @export
plot_spatial_annotation <- function(
    norm_list,
    annotation_col = "SingleR_label",
    method_order = NULL,
    point_size = 0.5,
    palette = NULL,
    ncol = 1,
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

    stop(
      "No valid normalisation methods were found."
    )
  }


  # ----------------------------------------------------------
  # Extract coordinates and annotation
  # ----------------------------------------------------------

  df_list <- list()


  for (method in methods) {

    spe <- norm_list[[method]]


    # --------------------------------------------------------
    # Check annotation column
    # --------------------------------------------------------

    cd <- SummarizedExperiment::colData(
      spe
    )


    if (
      !annotation_col %in%
      colnames(cd)
    ) {

      warning(
        "Annotation column '",
        annotation_col,
        "' not found for method: ",
        method,
        call. = FALSE
      )

      next
    }


    # --------------------------------------------------------
    # Spatial coordinates
    # --------------------------------------------------------

    coords <- as.data.frame(
      SpatialExperiment::spatialCoords(
        spe
      )
    )


    if (ncol(coords) < 2) {

      stop(
        "spatialCoords(spe) must contain at least two columns."
      )
    }


    if (nrow(coords) != ncol(spe)) {

      stop(
        "Number of spatial coordinates does not match the number of columns in the SpatialExperiment."
      )
    }


    # --------------------------------------------------------
    # Construct data frame
    # --------------------------------------------------------

    temp <- data.frame(

      x = coords[, 1],

      y = coords[, 2],

      cell_type = as.character(
        cd[[annotation_col]]
      ),

      method = method,

      stringsAsFactors = FALSE
    )


    df_list[[method]] <- temp
  }


  # ----------------------------------------------------------
  # Check extracted data
  # ----------------------------------------------------------

  if (length(df_list) == 0) {

    stop(
      "No annotation data could be extracted."
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


  df$cell_type <- factor(
    df$cell_type
  )


  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$x,
      y = .data$y,
      colour = .data$cell_type
    )
  ) +

    ggplot2::geom_point(
      size = point_size
    ) +

    ggplot2::facet_wrap(
      ~ method,
      ncol = ncol
    ) +

    ggplot2::coord_fixed() +

    ggplot2::labs(
      x = "x",
      y = "y",
      colour = "Cell type"
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


  # ----------------------------------------------------------
  # Optional colour palette
  # ----------------------------------------------------------

  if (!is.null(palette)) {

    p <- p +

      ggplot2::scale_colour_manual(
        values = palette
      )
  }


  list(
    data = df,
    plot = p
  )
}



# ============================================================
# Cell-type proportion plot
# ============================================================


#' Plot cell-type proportions
#'
#' Plot cell-type proportions stored in the colData of annotated
#' SpatialExperiment objects across normalisation methods.
#'
#' @param norm_list Named list of annotated SpatialExperiment
#'   objects.
#' @param annotation_col Name of the colData column containing
#'   cell-type annotations.
#' @param method_order Optional character vector specifying the
#'   order of normalisation methods.
#' @param palette Optional named character vector specifying
#'   colours for cell types.
#' @param show Character specifying whether to show
#'   `"proportion"` or `"count"`.
#' @param theme_name Plot theme. One of `"classic"`, `"bw"`,
#'   `"minimal"`, or `"gray"`.
#'
#' @return A list containing summary data and plot.
#'
#' @export
plot_celltype_proportion <- function(
    norm_list,
    annotation_col = "SingleR_label",
    method_order = NULL,
    palette = NULL,
    show = "proportion",
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

    stop(
      "No valid normalisation methods were found."
    )
  }


  # ----------------------------------------------------------
  # Extract annotation labels
  # ----------------------------------------------------------

  df_list <- list()


  for (method in methods) {

    spe <- norm_list[[method]]


    cd <- SummarizedExperiment::colData(
      spe
    )


    # --------------------------------------------------------
    # Check annotation column
    # --------------------------------------------------------

    if (
      !annotation_col %in%
      colnames(cd)
    ) {

      warning(
        "Annotation column '",
        annotation_col,
        "' not found for method: ",
        method,
        call. = FALSE
      )

      next
    }


    # --------------------------------------------------------
    # Construct data frame
    # --------------------------------------------------------

    temp <- data.frame(

      method = method,

      cell_type = as.character(
        cd[[annotation_col]]
      ),

      stringsAsFactors = FALSE
    )


    df_list[[method]] <- temp
  }


  # ----------------------------------------------------------
  # Check extracted data
  # ----------------------------------------------------------

  if (length(df_list) == 0) {

    stop(
      "No annotation data could be extracted."
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
  # Remove missing labels
  # ----------------------------------------------------------

  df <- df[
    !is.na(df$cell_type) &
      df$cell_type != "",
    ,
    drop = FALSE
  ]


  if (nrow(df) == 0) {

    stop(
      "All annotation labels are missing."
    )
  }


  # ----------------------------------------------------------
  # Calculate counts and proportions
  # ----------------------------------------------------------

  summary_df <- df |>

    dplyr::count(
      .data$method,
      .data$cell_type,
      name = "n"
    ) |>

    dplyr::group_by(
      .data$method
    ) |>

    dplyr::mutate(

      proportion =
        .data$n /
        sum(.data$n)

    ) |>

    dplyr::ungroup()


  # ----------------------------------------------------------
  # Factor ordering
  # ----------------------------------------------------------

  summary_df$method <- factor(
    summary_df$method,
    levels = methods
  )


  summary_df$cell_type <- factor(
    summary_df$cell_type
  )


  # ==========================================================
  # Proportion
  # ==========================================================

  if (show == "proportion") {

    p <- ggplot2::ggplot(
      summary_df,
      ggplot2::aes(
        x = .data$method,
        y = .data$proportion,
        fill = .data$cell_type
      )
    ) +

      ggplot2::geom_col() +

      ggplot2::scale_y_continuous(
        labels =
          scales::percent_format()
      ) +

      ggplot2::labs(
        x = NULL,
        y = "Cell-type proportion",
        fill = "Cell type"
      )
  }


  # ==========================================================
  # Count
  # ==========================================================

  else if (show == "count") {

    p <- ggplot2::ggplot(
      summary_df,
      ggplot2::aes(
        x = .data$method,
        y = .data$n,
        fill = .data$cell_type
      )
    ) +

      ggplot2::geom_col() +

      ggplot2::labs(
        x = NULL,
        y = "Number of spots",
        fill = "Cell type"
      )
  }


  # ----------------------------------------------------------
  # Invalid option
  # ----------------------------------------------------------

  else {

    stop(
      "'show' must be either 'proportion' or 'count'."
    )
  }


  # ----------------------------------------------------------
  # Theme
  # ----------------------------------------------------------

  p <- p +

    .get_plot_theme(
      theme_name
    ) +

    ggplot2::theme(

      panel.grid =
        ggplot2::element_blank(),

      axis.text.x =
        ggplot2::element_text(
          angle = 45,
          hjust = 1
        ),

      legend.position = "right"
    )


  # ----------------------------------------------------------
  # Optional colour palette
  # ----------------------------------------------------------

  if (!is.null(palette)) {

    p <- p +

      ggplot2::scale_fill_manual(
        values = palette
      )
  }


  list(
    data = summary_df,
    plot = p
  )
}






# ============================================================
# Pairwise annotation agreement
# ============================================================


#' Calculate pairwise cell-type annotation agreement
#'
#' Calculate the proportion of spatial locations receiving the
#' same cell-type annotation between each pair of normalisation
#' methods.
#'
#' @param norm_list Named list of annotated SpatialExperiment
#'   objects.
#' @param annotation_col Name of the colData column containing
#'   cell-type annotations.
#' @param method_order Optional character vector specifying the
#'   order of normalisation methods.
#' @param remove_na Logical. If TRUE, spatial locations with a
#'   missing annotation in either method are removed before
#'   calculating agreement.
#'
#' @return A data frame containing pairwise annotation agreement.
#'
#' @export
calculate_annotation_agreement <- function(
    norm_list,
    annotation_col = "SingleR_label",
    method_order = NULL,
    remove_na = TRUE
) {

  # ----------------------------------------------------------
  # Resolve methods
  # ----------------------------------------------------------

  methods <- .resolve_methods(
    norm_list,
    method_order
  )

  if (length(methods) < 2) {

    stop(
      "At least two normalisation methods are required."
    )
  }


  # ----------------------------------------------------------
  # Extract annotations
  # ----------------------------------------------------------

  annotation_list <- lapply(
    methods,
    function(method) {

      spe <- norm_list[[method]]

      cd <- SummarizedExperiment::colData(
        spe
      )

      if (
        !annotation_col %in%
        colnames(cd)
      ) {

        stop(
          "Annotation column '",
          annotation_col,
          "' not found for method: ",
          method
        )
      }

      as.character(
        cd[[annotation_col]]
      )
    }
  )

  names(annotation_list) <- methods


  # ----------------------------------------------------------
  # Check same number of spatial locations
  # ----------------------------------------------------------

  n_spots <- vapply(
    annotation_list,
    length,
    integer(1)
  )

  if (length(unique(n_spots)) != 1) {

    stop(
      "All methods must contain the same number of spatial locations."
    )
  }


  # ----------------------------------------------------------
  # Pairwise combinations
  # ----------------------------------------------------------

  method_pairs <- utils::combn(
    methods,
    2,
    simplify = FALSE
  )


  agreement_list <- lapply(
    method_pairs,
    function(pair) {

      method1 <- pair[1]
      method2 <- pair[2]

      label1 <- annotation_list[[method1]]
      label2 <- annotation_list[[method2]]


      # ------------------------------------------------------
      # Remove missing labels if requested
      # ------------------------------------------------------

      if (remove_na) {

        keep <-
          !is.na(label1) &
          !is.na(label2) &
          label1 != "" &
          label2 != ""

        label1 <- label1[keep]
        label2 <- label2[keep]
      }


      if (length(label1) == 0) {

        agreement <- NA_real_
        n_compared <- 0L

      } else {

        agreement <- mean(
          label1 == label2
        )

        n_compared <- length(
          label1
        )
      }


      data.frame(
        method1 = method1,
        method2 = method2,
        agreement = agreement,
        n_compared = n_compared,
        stringsAsFactors = FALSE
      )
    }
  )


  do.call(
    rbind,
    agreement_list
  )
}




# ============================================================
# Pairwise annotation agreement heatmap
# ============================================================


#' Plot pairwise cell-type annotation agreement
#'
#' Plot pairwise annotation agreement between normalisation
#' methods as a heatmap.
#'
#' @param norm_list Named list of annotated SpatialExperiment
#'   objects.
#' @param annotation_col Name of the colData column containing
#'   cell-type annotations.
#' @param method_order Optional character vector specifying the
#'   order of normalisation methods.
#' @param remove_na Logical. If TRUE, missing annotations are
#'   removed before agreement is calculated.
#' @param show_values Logical. If TRUE, agreement values are
#'   displayed inside heatmap cells.
#' @param digits Number of decimal places shown in heatmap cells.
#' @param colour_option Viridis colour option.
#' @param theme_name Plot theme. One of `"classic"`, `"bw"`,
#'   `"minimal"`, or `"gray"`.
#'
#' @return A list containing agreement data and plot.
#'
#' @export
plot_annotation_agreement <- function(
    norm_list,
    annotation_col = "SingleR_label",
    method_order = NULL,
    remove_na = TRUE,
    show_values = TRUE,
    digits = 2,
    colour_option = "viridis",
    theme_name = "classic"
) {

  # ----------------------------------------------------------
  # Resolve methods
  # ----------------------------------------------------------

  methods <- .resolve_methods(
    norm_list,
    method_order
  )


  if (length(methods) < 2) {

    stop(
      "At least two normalisation methods are required."
    )
  }


  # ----------------------------------------------------------
  # Calculate pairwise agreement
  # ----------------------------------------------------------

  agreement_df <- calculate_annotation_agreement(
    norm_list = norm_list,
    annotation_col = annotation_col,
    method_order = methods,
    remove_na = remove_na
  )


  # ----------------------------------------------------------
  # Add diagonal
  # ----------------------------------------------------------

  diagonal_df <- data.frame(

    method1 = methods,

    method2 = methods,

    agreement = 1,

    n_compared = vapply(
      norm_list[methods],
      ncol,
      integer(1)
    ),

    stringsAsFactors = FALSE
  )


  # ----------------------------------------------------------
  # Add reversed comparisons
  # ----------------------------------------------------------

  reverse_df <- data.frame(

    method1 = agreement_df$method2,

    method2 = agreement_df$method1,

    agreement = agreement_df$agreement,

    n_compared = agreement_df$n_compared,

    stringsAsFactors = FALSE
  )


  # ----------------------------------------------------------
  # Combine
  # ----------------------------------------------------------

  heatmap_df <- rbind(
    agreement_df,
    reverse_df,
    diagonal_df
  )


  # ----------------------------------------------------------
  # Factor ordering
  # ----------------------------------------------------------

  heatmap_df$method1 <- factor(
    heatmap_df$method1,
    levels = methods
  )

  heatmap_df$method2 <- factor(
    heatmap_df$method2,
    levels = methods
  )


  # ----------------------------------------------------------
  # Create labels BEFORE constructing ggplot
  # ----------------------------------------------------------

  heatmap_df$label <- ifelse(
    is.na(heatmap_df$agreement),
    "",
    formatC(
      heatmap_df$agreement,
      format = "f",
      digits = digits
    )
  )


  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------

  p <- ggplot2::ggplot(
    heatmap_df,
    ggplot2::aes(
      x = .data$method1,
      y = .data$method2,
      fill = .data$agreement
    )
  ) +

    ggplot2::geom_tile() +

    viridis::scale_fill_viridis(
      option = colour_option,
      limits = c(0, 1),
      name = "Agreement"
    ) +

    ggplot2::coord_fixed() +

    ggplot2::labs(
      x = NULL,
      y = NULL
    ) +

    .get_plot_theme(
      theme_name
    ) +

    ggplot2::theme(

      panel.grid =
        ggplot2::element_blank(),

      axis.text.x =
        ggplot2::element_text(
          angle = 45,
          hjust = 1
        ),

      legend.position = "right"
    )


  # ----------------------------------------------------------
  # Add agreement values
  # ----------------------------------------------------------

  if (show_values) {

    p <- p +

      ggplot2::geom_text(
        ggplot2::aes(
          label = .data$label
        ),
        size = 3
      )
  }


  # ----------------------------------------------------------
  # Return
  # ----------------------------------------------------------

  list(
    data = heatmap_df,
    pairwise = agreement_df,
    plot = p
  )
}
