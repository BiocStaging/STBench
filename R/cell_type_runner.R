# Cell-type annotation


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
#' @examples
#' data(norm_example)
#'
#' # Use a subset of cells as an illustrative reference
#' ref <- norm_example$None[, seq_len(100)]
#' labels <- as.character(ref$cell_type)
#'
#' # Annotate the remaining cells
#' query <- lapply(
#'   norm_example,
#'   function(spe) spe[, seq.int(101, ncol(spe))]
#' )
#'
#' annotated_norm <- run_annotation(
#'   norm_list = query,
#'   ref = ref,
#'   labels = labels,
#'   method = "SingleR",
#'   label_name = "SingleR_label",
#'   pruned_label_name = "SingleR_pruned_label",
#'   de.method = "classic"
#' )
#'
#' table(annotated_norm$None$SingleR_label)
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
  # Check input


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


  # Run annotation


  out <- lapply(
    names(norm_list),
    function(method_name) {
      message("Running SingleR: ", method_name)

      spe <- norm_list[[method_name]]

      # Run SingleR


      results <- SingleR::SingleR(
        test = spe,
        ref = ref,
        labels = labels,
        ...
      )


      # Check number of predictions


      if (nrow(results) != ncol(spe)) {
        stop(
          "Number of SingleR predictions does not match the number of spatial locations."
        )
      }


      # Add labels to colData


      SummarizedExperiment::colData(spe)[[label_name]] <-
        as.character(
          results$labels
        )


      # Add pruned labels


      if ("pruned.labels" %in% colnames(results)) {
        SummarizedExperiment::colData(spe)[[pruned_label_name]] <-
          as.character(
            results$pruned.labels
          )
      }


      # Return SpatialExperiment


      spe
    }
  )


  names(out) <- names(
    norm_list
  )


  out
}


# Spatial cell-type annotation plot


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
#' @param palette Colour palette. One of "viridis", "magma",
#'   "inferno", "plasma", "cividis", "rocket", "mako",
#'   "turbo", or a custom character vector of colours.
#'   If NULL, the default ggplot2 palette is used.
#' @param ncol Number of facet columns.
#' @param theme_name Plot theme. One of "classic", "bw",
#'   "minimal", or "gray".
#'
#' @return A list containing plotting data and plot.
#'
#' @examples
#' data(norm_example)
#'
#' p <- plot_spatial_annotation(
#'   norm_list = norm_example,
#'   annotation_col = "cell_type",
#'   point_size = 0.5,
#'   theme_name = "classic"
#' )
#'
#' p$plot
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
  # Resolve methods

  methods <- .resolve_methods(
    norm_list,
    method_order
  )

  if (length(methods) == 0) {
    stop(
      "No valid normalisation methods were found."
    )
  }

  # Extract coordinates and annotation

  df_list <- list()

  for (method in methods) {
    spe <- norm_list[[method]]

    # Check annotation column

    cd <- SummarizedExperiment::colData(spe)

    if (!annotation_col %in% colnames(cd)) {
      warning(
        "Annotation column '",
        annotation_col,
        "' not found for method: ",
        method,
        call. = FALSE
      )

      next
    }

    # Spatial coordinates

    coords <- as.data.frame(
      SpatialExperiment::spatialCoords(spe)
    )

    if (ncol(coords) < 2) {
      stop(
        "spatialCoords(spe) must contain at least two columns."
      )
    }

    if (nrow(coords) != ncol(spe)) {
      stop(
          "Number of spatial coordinates does not match ",
          "the number of columns in the SpatialExperiment."
          )
    }

    # Construct data frame

    temp <- data.frame(
      x = coords[, 1],
      y = coords[, 2],
      cell_type = as.character(cd[[annotation_col]]),
      method = method,
      stringsAsFactors = FALSE
    )

    df_list[[method]] <- temp
  }

  # Check extracted data

  if (length(df_list) == 0) {
    stop(
      "No annotation data could be extracted."
    )
  }

  # Combine data

  df <- do.call(
    rbind,
    df_list
  )

  rownames(df) <- NULL

  # Factor ordering

  df$method <- factor(
    df$method,
    levels = methods
  )

  df$cell_type <- factor(
    df$cell_type
  )

  # Plot

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
      ~method,
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

  # Optional colour palette

  if (!is.null(palette)) {
    viridis_options <- c(
      "viridis" = "D",
      "magma" = "A",
      "inferno" = "B",
      "plasma" = "C",
      "cividis" = "E",
      "rocket" = "F",
      "mako" = "G",
      "turbo" = "H"
    )

    if (
      is.character(palette) &&
        length(palette) == 1 &&
        !is.na(palette) &&
        palette %in% names(viridis_options)
    ) {
      p <- p +
        ggplot2::scale_colour_viridis_d(
          option = viridis_options[[palette]],
          na.value = "grey80"
        )
    } else {
      p <- p +
        ggplot2::scale_colour_manual(
          values = palette,
          na.value = "grey80"
        )
    }
  }

  # Return
  p <- p +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(size = 3)
      )
    )
  list(
    data = df,
    plot = p
  )
}


# Cell-type proportion plot


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
#'
#' @return A list containing summary data and plot.
#' @examples
#' data(norm_example)
#'
#' p <- plot_celltype_proportion(
#'   norm_list = norm_example,
#'   annotation_col = "cell_type",
#'   method_order = NULL,
#'   palette = "turbo",
#'   show = "proportion",
#'   theme_name = "classic"
#' )
#'
#' p$plot
#' @export
plot_celltype_proportion <- function(
  norm_list,
  annotation_col = "SingleR_label",
  method_order = NULL,
  palette = NULL,
  show = "proportion",
  theme_name = "classic"
) {
  # Resolve methods

  methods <- .resolve_methods(
    norm_list,
    method_order
  )

  if (length(methods) == 0) {
    stop(
      "No valid normalisation methods were found."
    )
  }

  # Extract annotation labels

  df_list <- list()

  for (method in methods) {
    spe <- norm_list[[method]]

    cd <- SummarizedExperiment::colData(spe)

    # Check annotation column

    if (!annotation_col %in% colnames(cd)) {
      warning(
        "Annotation column '",
        annotation_col,
        "' not found for method: ",
        method,
        call. = FALSE
      )

      next
    }

    # Construct data frame

    temp <- data.frame(
      method = method,
      cell_type = as.character(
        cd[[annotation_col]]
      ),
      stringsAsFactors = FALSE
    )

    df_list[[method]] <- temp
  }

  # Check extracted data

  if (length(df_list) == 0) {
    stop(
      "No annotation data could be extracted."
    )
  }

  # Combine

  df <- do.call(
    rbind,
    df_list
  )

  rownames(df) <- NULL

  # Remove missing labels

  df <- df[
    !is.na(df$cell_type) &
      df$cell_type != "", ,
    drop = FALSE
  ]

  if (nrow(df) == 0) {
    stop(
      "All annotation labels are missing."
    )
  }

  # Calculate counts and proportions

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
      proportion = .data$n / sum(.data$n)
    ) |>
    dplyr::ungroup()

  # Factor ordering

  summary_df$method <- factor(
    summary_df$method,
    levels = methods
  )

  summary_df$cell_type <- factor(
    summary_df$cell_type
  )

  # Proportion plot

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
        labels = scales::percent_format()
      ) +
      ggplot2::labs(
        x = NULL,
        y = "Cell-type proportion",
        fill = "Cell type"
      )
  } else if (show == "count") {
    # Count plot

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
  } else {
    stop(
      "'show' must be either 'proportion' or 'count'."
    )
  }

  # Theme

  p <- p +
    .get_plot_theme(
      theme_name
    ) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      legend.position = "right"
    )

  # Optional colour palette

  if (!is.null(palette)) {
    viridis_options <- c(
      "viridis" = "D",
      "magma" = "A",
      "inferno" = "B",
      "plasma" = "C",
      "cividis" = "E",
      "rocket" = "F",
      "mako" = "G",
      "turbo" = "H"
    )

    if (
      is.character(palette) &&
        length(palette) == 1 &&
        !is.na(palette) &&
        palette %in% names(viridis_options)
    ) {
      p <- p +
        ggplot2::scale_fill_viridis_d(
          option = viridis_options[[palette]],
          na.value = "grey80"
        )
    } else {
      p <- p +
        ggplot2::scale_fill_manual(
          values = palette,
          na.value = "grey80"
        )
    }
  }

  # Return

  list(
    data = summary_df,
    plot = p
  )
}


# Pairwise annotation agreement


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
#' @examples
#' data(norm_example)
#'
#' agreement <- calculate_annotation_agreement(
#'   norm_list = norm_example,
#'   annotation_col = "cell_type",
#'   method_order = NULL,
#'   remove_na = TRUE
#' )
#'
#' agreement
#' @export
calculate_annotation_agreement <- function(
  norm_list,
  annotation_col = "SingleR_label",
  method_order = NULL,
  remove_na = TRUE
) {
  # Resolve methods


  methods <- .resolve_methods(
    norm_list,
    method_order
  )

  if (length(methods) < 2) {
    stop(
      "At least two normalisation methods are required."
    )
  }


  # Extract annotations


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


  # Check same number of spatial locations


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


  # Pairwise combinations


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


      # Remove missing labels if requested


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


# Pairwise annotation agreement heatmap


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
  # Resolve methods


  methods <- .resolve_methods(
    norm_list,
    method_order
  )


  if (length(methods) < 2) {
    stop(
      "At least two normalisation methods are required."
    )
  }


  # Calculate pairwise agreement


  agreement_df <- calculate_annotation_agreement(
    norm_list = norm_list,
    annotation_col = annotation_col,
    method_order = methods,
    remove_na = remove_na
  )


  # Add diagonal


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


  # Add reversed comparisons


  reverse_df <- data.frame(
    method1 = agreement_df$method2,
    method2 = agreement_df$method1,
    agreement = agreement_df$agreement,
    n_compared = agreement_df$n_compared,
    stringsAsFactors = FALSE
  )


  # Combine


  heatmap_df <- rbind(
    agreement_df,
    reverse_df,
    diagonal_df
  )


  # Factor ordering


  heatmap_df$method1 <- factor(
    heatmap_df$method1,
    levels = methods
  )

  heatmap_df$method2 <- factor(
    heatmap_df$method2,
    levels = methods
  )


  # Create labels BEFORE constructing ggplot


  heatmap_df$label <- ifelse(
    is.na(heatmap_df$agreement),
    "",
    formatC(
      heatmap_df$agreement,
      format = "f",
      digits = digits
    )
  )


  # Plot


  p <- ggplot2::ggplot(
    heatmap_df,
    ggplot2::aes(
      x = .data$method1,
      y = .data$method2,
      fill = .data$agreement
    )
  ) +
    ggplot2::geom_tile(alpha = 0.8) +
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


  # Add agreement values


  if (show_values) {
    p <- p +

      ggplot2::geom_text(
        ggplot2::aes(
          label = .data$label
        ),
        size = 3
      )
  }


  # Return


  list(
    data = heatmap_df,
    pairwise = agreement_df,
    plot = p
  )
}


.extract_marker_expression <- function(
  norm_list,
  markers,
  annotation_col,
  method_order = NULL,
  reference_method = NULL
) {
  if (is.null(names(norm_list)) ||
    anyNA(names(norm_list)) ||
    any(!nzchar(names(norm_list))) ||
    anyDuplicated(names(norm_list))) {
    stop("norm_list must have unique, non-empty method names.")
  }

  if (is.null(method_order)) {
    methods <- names(norm_list)
  } else {
    if (anyDuplicated(method_order) ||
      !all(method_order %in% names(norm_list))) {
      stop("method_order contains duplicate or unknown methods.")
    }
    methods <- method_order
  }

  if (length(methods) == 0) {
    stop("No normalisation methods selected.")
  }

  if (length(markers) == 0 || anyDuplicated(markers)) {
    stop("markers must contain unique gene names.")
  }

  # Optional fixed reference annotations

  ref_labels <- NULL

  if (!is.null(reference_method)) {
    if (!reference_method %in% names(norm_list)) {
      stop("reference_method not found in norm_list.")
    }

    ref_spe <- norm_list[[reference_method]]

    ref_cd <- as.data.frame(
      SummarizedExperiment::colData(ref_spe)
    )

    if (!annotation_col %in% names(ref_cd)) {
      stop("Annotation column missing from reference method.")
    }

    ref_ids <- colnames(ref_spe)

    if (is.null(ref_ids) ||
      anyNA(ref_ids) ||
      anyDuplicated(ref_ids)) {
      stop("Reference must have unique cell identifiers.")
    }

    ref_labels <- stats::setNames(
      as.character(ref_cd[[annotation_col]]),
      ref_ids
    )
  }

  # Extract expression

  df_list <- lapply(methods, function(method) {
    spe <- norm_list[[method]]

    missing_genes <- setdiff(markers, rownames(spe))

    if (length(missing_genes) > 0) {
      stop(
        "Missing genes in ", method, ": ",
        toString(missing_genes)
      )
    }

    if (!"logcounts" %in%
      SummarizedExperiment::assayNames(spe)) {
      stop("logcounts assay missing in method: ", method)
    }

    cell_ids <- colnames(spe)

    if (is.null(cell_ids) ||
      anyNA(cell_ids) ||
      anyDuplicated(cell_ids)) {
      stop("Unique cell identifiers required: ", method)
    }

    # Get annotation labels

    if (!is.null(reference_method)) {
      if (!setequal(cell_ids, names(ref_labels))) {
        stop(
          "Cell identifiers do not match the reference in: ",
          method
        )
      }

      cell_types <- unname(ref_labels[cell_ids])
    } else {
      cd <- as.data.frame(
        SummarizedExperiment::colData(spe)
      )

      if (!annotation_col %in% names(cd)) {
        stop("Annotation column missing in: ", method)
      }

      cell_types <- as.character(cd[[annotation_col]])
    }

    # Extract marker expression

    exprs <- SummarizedExperiment::assay(
      spe,
      "logcounts"
    )[markers, , drop = FALSE]

    exprs <- as.data.frame(
      t(as.matrix(exprs))
    )

    exprs$cell_id <- cell_ids
    exprs$cell_type <- cell_types
    exprs$method <- method

    tidyr::pivot_longer(
      exprs,
      cols = tidyselect::all_of(markers),
      names_to = "gene",
      values_to = "expression"
    )
  })

  df <- dplyr::bind_rows(df_list)

  df <- df[
    !is.na(df$cell_type) &
      nzchar(df$cell_type) &
      is.finite(df$expression), ,
    drop = FALSE
  ]

  if (nrow(df) == 0) {
    stop("No valid marker-expression data.")
  }

  df$method <- factor(
    df$method,
    levels = methods
  )

  df$gene <- factor(
    df$gene,
    levels = markers
  )

  df$cell_type <- factor(df$cell_type)

  df
}


#' Plot marker expression across normalisation methods
#'
#' Compare marker-gene expression distributions across
#' selected cell types and normalisation methods.
#'
#' @param norm_list Named list of SpatialExperiment objects.
#' @param markers Character vector of marker genes.
#' @param annotation_col Cell-type annotation column.
#' @param methods Optional character vector selecting
#'   normalisation methods. NULL includes all methods.
#' @param method_order Optional character vector specifying
#'   the display order of selected methods.
#' @param cell_types Optional character vector selecting
#'   cell types. NULL includes all cell types.
#' @param reference_method Optional method used to define
#'   fixed cell-type annotations across methods.
#' @param scale_expression If TRUE, apply pooled gene-wise
#'   Z-score standardisation across methods and cells.
#' @param winsorize If TRUE, cap expression at the specified
#'   upper quantile before plotting or scaling.
#' @param winsor_quantile Upper winsorisation quantile.
#' @param layout Facet layout. Either "method_gene" or
#'   "gene_celltype".
#' @param palette Optional named colour vector or viridis
#'   palette name.
#' @param box_width Width of boxplots.
#' @param theme_name Plot theme.
#'
#' @return A list containing plotting data and ggplot.
#'
#' @examples
#' data(norm_example)
#'
#' markers <- head(
#'   rownames(norm_example[[1]]),
#'   2
#' )
#'
#' p <- plot_marker_expression(
#'   norm_list = norm_example,
#'   annotation_col = "cell_type",
#'   markers = markers,
#'   palette = "turbo"
#' )
#'
#' p$plot
#' @export

plot_marker_expression <- function(
  norm_list,
  markers,
  annotation_col = "SingleR_label",
  methods = NULL,
  method_order = NULL,
  cell_types = NULL,
  reference_method = NULL,
  scale_expression = TRUE,
  winsorize = TRUE,
  winsor_quantile = 0.99,
  layout = "method_gene",
  palette = "viridis",
  box_width = 0.7,
  theme_name = "classic"
) {
  layout <- match.arg(
    layout,
    c("method_gene", "gene_celltype")
  )

  # Select normalisation methods

  available_methods <- names(norm_list)

  if (is.null(methods)) {
    methods <- available_methods
  }

  if (length(methods) == 0 ||
    anyDuplicated(methods) ||
    !all(methods %in% available_methods)) {
    stop(
      "methods must contain valid, unique ",
      "normalisation method names."
    )
  }

  # Set method order

  if (is.null(method_order)) {
    method_order <- methods
  } else {
    if (anyDuplicated(method_order) ||
      !setequal(method_order, methods)) {
      stop(
        "method_order must contain exactly ",
        "the selected methods."
      )
    }
  }

  # Validate winsorisation quantile

  if (length(winsor_quantile) != 1 ||
    !is.finite(winsor_quantile) ||
    winsor_quantile <= 0 ||
    winsor_quantile > 1) {
    stop(
      "winsor_quantile must be in (0, 1]."
    )
  }

  # Extract marker expression

  df <- .extract_marker_expression(
    norm_list = norm_list,
    markers = markers,
    annotation_col = annotation_col,
    method_order = method_order,
    reference_method = reference_method
  )

  # Winsorisation and scaling
  # Calculated before cell-type filtering

  df <- df |>
    dplyr::group_by(.data$gene) |>
    dplyr::mutate(
      expression_winsor = if (winsorize) {
        pmin(
          .data$expression,
          stats::quantile(
            .data$expression,
            winsor_quantile,
            na.rm = TRUE
          )
        )
      } else {
        .data$expression
      }
    ) |>
    dplyr::mutate(
      expression_scaled = if (scale_expression) {
        s <- stats::sd(.data$expression_winsor)

        if (is.na(s) || s == 0) {
          rep(0, dplyr::n())
        } else {
          (
            .data$expression_winsor -
              mean(.data$expression_winsor)
          ) / s
        }
      } else {
        .data$expression_winsor
      }
    ) |>
    dplyr::ungroup()

  # Select cell types

  if (!is.null(cell_types)) {
    available_cell_types <- unique(
      as.character(df$cell_type)
    )

    missing_types <- setdiff(
      cell_types,
      available_cell_types
    )

    if (length(missing_types) > 0) {
      stop(
        "Unknown cell types: ",
        toString(missing_types)
      )
    }

    df <- df |>
      dplyr::filter(
        as.character(.data$cell_type) %in% cell_types
      )

    df$cell_type <- factor(
      df$cell_type,
      levels = cell_types
    )
  } else {
    df$cell_type <- droplevels(
      df$cell_type
    )
  }

  # Method ordering

  df$method <- factor(
    df$method,
    levels = method_order
  )

  # Plot

  if (layout == "method_gene") {
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data$cell_type,
        y = .data$expression_scaled,
        fill = .data$cell_type
      )
    ) +
      ggplot2::facet_grid(
        method ~ gene
      )
  } else {
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data$method,
        y = .data$expression_scaled,
        fill = .data$method
      )
    ) +
      ggplot2::facet_grid(
        gene ~ cell_type
      )
  }

  p <- p +
    ggplot2::geom_boxplot(
      outlier.shape = NA, alpha = 0.7,
      width = box_width
    ) +
    ggplot2::labs(
      x = NULL,
      y = if (scale_expression) {
        "Gene-wise standardised log expression"
      } else {
        "Log expression"
      },
      fill = "Normalisation"
    ) +
    .get_plot_theme(theme_name) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      strip.background = ggplot2::element_blank(),
      legend.position = "none"
    )

  # Optional colour palette

  if (!is.null(palette)) {
    options <- c(
      viridis = "D",
      magma = "A",
      inferno = "B",
      plasma = "C",
      cividis = "E",
      rocket = "F",
      mako = "G",
      turbo = "H"
    )

    if (is.character(palette) &&
      length(palette) == 1 &&
      !is.na(palette) &&
      palette %in% names(options)) {
      p <- p +
        ggplot2::scale_fill_viridis_d(
          option = options[[palette]]
        )
    } else {
      p <- p +
        ggplot2::scale_fill_manual(
          values = palette
        )
    }
  }

  list(
    data = df,
    plot = p
  )
}


#' Plot marker-gene discrimination across methods
#'
#' Calculate ROC AUC for marker genes across
#' normalisation methods using target cell-type labels.
#'
#' @param norm_list Named list of SpatialExperiment objects.
#' @param markers Character vector of marker genes.
#' @param target_types Named character vector mapping marker
#'   genes to their expected target cell types.
#' @param annotation_col Cell-type annotation column.
#' @param method_order Optional normalisation method order.
#' @param reference_method Optional method used to define
#'   fixed cell-type labels across methods.
#' @param palette Optional named colour vector or
#'   viridis palette name.
#' @param point_size Size of AUC points.
#' @param ncol number of column of the panel
#' @param theme_name Plot theme.
#'
#' @return List containing AUC summary data and ggplot.
#'
#' @examples
#' data(norm_example)
#'
#' spe <- norm_example[[1]]
#'
#' labels <- as.character(
#'   SummarizedExperiment::colData(spe)$cell_type
#' )
#'
#' label_counts <- sort(
#'   table(labels),
#'   decreasing = TRUE
#' )
#'
#' target_label <- names(label_counts)[1]
#'
#' markers <- head(rownames(spe), 2)
#'
#' target_types <- stats::setNames(
#'   rep(target_label, length(markers)),
#'   markers
#' )
#'
#' p <- plot_marker_auc(
#'   norm_list = norm_example,
#'   markers = markers,
#'   target_types = target_types,
#'   annotation_col = "cell_type",
#'   reference_method = "None",
#'   palette = "magma",
#'   point_size = 3
#' )
#'
#' p$plot
#' @export

plot_marker_auc <- function(
  norm_list,
  markers,
  target_types,
  annotation_col = "SingleR_label",
  method_order = NULL,
  reference_method = NULL,
  palette = "viridis",
  point_size = 3,
  theme_name = "classic",
  ncol = 2
) {
  # Validate target cell types

  if (is.null(names(target_types)) ||
    anyDuplicated(names(target_types)) ||
    !all(markers %in% names(target_types))) {
    stop(
      "target_types must be a named vector ",
      "containing one target cell type per marker."
    )
  }

  target_types <- target_types[markers]

  if (anyNA(target_types) ||
    any(!nzchar(target_types))) {
    stop("target_types contains missing or empty labels.")
  }

  # Extract expression

  df <- .extract_marker_expression(
    norm_list = norm_list,
    markers = markers,
    annotation_col = annotation_col,
    method_order = method_order,
    reference_method = reference_method
  )

  # Define target versus other cell types

  df$target <- as.character(df$cell_type) ==
    unname(target_types[as.character(df$gene)])

  # Calculate AUC

  auc_df <- df |>
    dplyr::group_by(
      .data$method,
      .data$gene
    ) |>
    dplyr::summarise(
      n_target = sum(.data$target),
      n_other = sum(!.data$target),
      AUC = {
        y <- .data$target
        x <- .data$expression

        if (all(y) || !any(y)) {
          NA_real_
        } else {
          roc_obj <- pROC::roc(
            response = y,
            predictor = x,
            levels = c(FALSE, TRUE),
            direction = "<",
            quiet = TRUE
          )

          as.numeric(pROC::auc(roc_obj))
        }
      },
      .groups = "drop"
    )

  # Plot

  p <- ggplot2::ggplot(
    auc_df,
    ggplot2::aes(
      x = .data$method,
      y = .data$AUC,
      colour = .data$method
    )
  ) +
    ggplot2::geom_point(
      size = point_size,
      na.rm = TRUE
    ) +
    ggplot2::facet_wrap(
      ~gene,
      scales = "free", ncol = ncol
    ) +
    # ggplot2::scale_y_continuous(
    #   limits = c(0, 1),
    #   breaks = seq(0, 1, 0.2)
    # ) +
    ggplot2::labs(
      x = "Normalisation method",
      y = "AUC value",
      colour = "Normalisation"
    ) +
    .get_plot_theme(theme_name) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      strip.background = ggplot2::element_blank(),
      legend.position = "none"
    )


  list(
    data = auc_df,
    plot = p
  )
}
