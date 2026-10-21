
# ============================================================
# UMAP comparison across normalisation methods
# ============================================================


#' Calculate UMAP embeddings across normalisation methods
#'
#' Perform PCA and UMAP independently for each normalisation
#' method using the same sampled spots and common genes.
#'
#' All colData annotations from the reference object are
#' retained for subsequent plotting.
#'
#' @param norm Named list of normalised SpatialExperiment objects.
#' @param methods Character vector of methods to compare.
#'   If NULL, all available methods are used.
#' @param prop_spots Proportion of spots to sample.
#' @param seed Random seed.
#' @param assay_name Assay used for PCA and UMAP.
#' @param n_pcs Number of principal components.
#' @param n_neighbors Number of UMAP neighbours.
#'
#' @return Data frame containing UMAP coordinates, method,
#'   spot identifiers, and all colData annotations.
#'
#' @export
get_umap_data <- function(
    norm,
    methods = NULL,
    prop_spots = 1,
    seed = 123,
    assay_name = "logcounts",
    n_pcs = 30,
    n_neighbors = 15
) {

  # ----------------------------------------------------------
  # Select normalisation methods
  # ----------------------------------------------------------

  .check_named_list(norm)

  if (is.null(methods)) {
    methods <- names(norm)
  }

  unknown <- setdiff(
    methods,
    names(norm)
  )

  if (length(unknown) > 0) {
    stop(
      "Unknown normalisation methods: ",
      paste(unknown, collapse = ", ")
    )
  }

  methods <- unique(methods)

  if (length(methods) == 0) {
    stop("No normalisation methods selected.")
  }


  # ----------------------------------------------------------
  # Check sampling proportion
  # ----------------------------------------------------------

  if (
    length(prop_spots) != 1 ||
    !is.numeric(prop_spots) ||
    is.na(prop_spots) ||
    prop_spots <= 0 ||
    prop_spots > 1
  ) {
    stop("prop_spots must be between 0 and 1.")
  }


  # ----------------------------------------------------------
  # Reference object
  # ----------------------------------------------------------

  reference <- norm[[methods[1]]]


  # ----------------------------------------------------------
  # Identify common spots
  # ----------------------------------------------------------

  common_spots <- Reduce(
    intersect,
    lapply(
      norm[methods],
      colnames
    )
  )


  # ----------------------------------------------------------
  # Identify common genes
  # ----------------------------------------------------------

  common_genes <- Reduce(
    intersect,
    lapply(
      norm[methods],
      rownames
    )
  )


  if (length(common_spots) < 4) {
    stop("At least four common spots are required.")
  }

  if (length(common_genes) < 3) {
    stop("At least three common genes are required.")
  }


  # ----------------------------------------------------------
  # Sample identical spots across methods
  # ----------------------------------------------------------

  set.seed(seed)

  n_sample <- max(
    4,
    floor(
      prop_spots * length(common_spots)
    )
  )

  sampled_spots <- sample(
    common_spots,
    size = n_sample,
    replace = FALSE
  )


  # ----------------------------------------------------------
  # Extract all reference annotations
  # ----------------------------------------------------------

  annotation_df <- as.data.frame(
    SummarizedExperiment::colData(reference)
  )

  annotation_df <- annotation_df[
    match(
      sampled_spots,
      colnames(reference)
    ),
    ,
    drop = FALSE
  ]


  # ----------------------------------------------------------
  # Run PCA and UMAP for each method
  # ----------------------------------------------------------

  umap_list <- lapply(
    methods,
    function(method) {

      message("Running PCA and UMAP: ", method)

      spe <- norm[[method]]


      # ----------------------------------------------
      # Subset to common genes and sampled spots
      # ----------------------------------------------

      spe <- spe[
        common_genes,
        sampled_spots,
        drop = FALSE
      ]


      # ----------------------------------------------
      # Check expression assay
      # ----------------------------------------------

      if (
        !assay_name %in%
        SummarizedExperiment::assayNames(spe)
      ) {
        stop(
          "Assay '",
          assay_name,
          "' not found for method: ",
          method
        )
      }


      # ----------------------------------------------
      # Run PCA
      # ----------------------------------------------

      set.seed(seed)

      spe <- scater::runPCA(
        spe,
        exprs_values = assay_name,
        ncomponents = min(
          n_pcs,
          nrow(spe) - 1,
          ncol(spe) - 1
        )
      )


      # ----------------------------------------------
      # Run UMAP
      # ----------------------------------------------

      set.seed(seed)

      spe <- scater::runUMAP(
        spe,
        dimred = "PCA",
        n_neighbors = min(
          n_neighbors,
          ncol(spe) - 1
        )
      )


      # ----------------------------------------------
      # Extract UMAP coordinates
      # ----------------------------------------------

      umap_coords <- as.data.frame(
        SingleCellExperiment::reducedDim(
          spe,
          "UMAP"
        )
      )

      colnames(umap_coords) <- c(
        "UMAP1",
        "UMAP2"
      )


      # ----------------------------------------------
      # Add all annotations
      # ----------------------------------------------

      umap_coords$method <- method

      umap_coords$spot <- colnames(spe)

      umap_coords <- cbind(
        umap_coords,
        annotation_df
      )

      umap_coords
    }
  )


  # ----------------------------------------------------------
  # Combine results
  # ----------------------------------------------------------

  df <- dplyr::bind_rows(
    umap_list
  )

  df$method <- factor(
    df$method,
    levels = methods
  )


  # ----------------------------------------------------------
  # Return
  # ----------------------------------------------------------

  return(df)
}



# ============================================================
# Plot previously calculated UMAP embeddings
# ============================================================


#' Plot UMAP comparisons across normalisation methods
#'
#' Visualise previously calculated UMAP embeddings using one
#' or multiple colData annotations without recalculating PCA
#' or UMAP.
#'
#' A single annotation uses facet_wrap by method.
#' Multiple annotations use facet_grid with methods in rows
#' and annotations in columns.
#'
#' @param umap_data Data frame returned by get_umap_data().
#' @param colour_by Character vector containing one or more
#'   annotation columns.
#' @param point_size Point size for UMAP plots.
#' @param alpha Point transparency.
#' @param ncol Number of columns for single-annotation facets.
#' @param theme_name Plot theme.
#'
#' @return A list containing plotting data and ggplot.
#'
#' @export
plot_umap_comparison <- function(
    umap_data,
    colour_by = "region",
    point_size = 0.58,
    alpha = 0.6,
    ncol = 3,
    theme_name = "classic"
) {

  # ----------------------------------------------------------
  # Check annotations
  # ----------------------------------------------------------

  colour_by <- unique(colour_by)

  if (length(colour_by) == 0) {
    stop("At least one annotation column is required.")
  }

  missing_cols <- setdiff(
    colour_by,
    colnames(umap_data)
  )

  if (length(missing_cols) > 0) {
    stop(
      "Annotation column(s) not found: ",
      paste(missing_cols, collapse = ", ")
    )
  }


  # ----------------------------------------------------------
  # One annotation: facet_wrap
  # ----------------------------------------------------------

  if (length(colour_by) == 1) {

    colour_col <- colour_by[1]

    df <- umap_data

    df$annotation <- as.factor(
      df[[colour_col]]
    )

    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data$UMAP1,
        y = .data$UMAP2,
        colour = .data$annotation
      )
    ) +
      ggplot2::geom_point(
        size = point_size,
        alpha = alpha,
        stroke = 0
      ) +
      ggplot2::facet_wrap(
        ~ method,
        ncol = ncol,
        scales = "free"
      ) +
      ggplot2::scale_colour_viridis_d(
        option = "turbo",
        end = 0.95,
        na.value = "grey80"
      ) +
      ggplot2::labs(
        x = "UMAP 1",
        y = "UMAP 2",
        colour = colour_col
      ) +
      .get_plot_theme(theme_name) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(
          colour = "black",
          fill = NA,
          linewidth = 0.5
        ),
        strip.text = ggplot2::element_text(
          face = "bold"
        ),
        legend.position = "bottom"
      ) +
      ggplot2::guides(
        colour = ggplot2::guide_legend(
          override.aes = list(
            size = 4,
            alpha = 1
          ),
          nrow = 1
        )
      )

    return(
      list(
        data = df,
        plot = p
      )
    )
  }


  # ----------------------------------------------------------
  # Multiple annotations: prepare long data
  # ----------------------------------------------------------

  df_list <- lapply(
    colour_by,
    function(col) {

      temp <- umap_data[
        ,
        c(
          "UMAP1",
          "UMAP2",
          "method",
          "spot"
        ),
        drop = FALSE
      ]

      temp$annotation_type <- col

      temp$annotation <- as.character(
        umap_data[[col]]
      )

      temp
    }
  )

  df <- dplyr::bind_rows(
    df_list
  )

  df$annotation_type <- factor(
    df$annotation_type,
    levels = colour_by
  )

  df$annotation <- factor(
    df$annotation
  )


  # ----------------------------------------------------------
  # Multiple annotations: facet_grid
  # ----------------------------------------------------------

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$UMAP1,
      y = .data$UMAP2,
      colour = .data$annotation
    )
  ) +
    ggplot2::geom_point(
      size = point_size,
      alpha = alpha,
      stroke = 0
    ) +
    ggplot2::facet_grid(
      method ~ annotation_type,
      scales = "free"
    ) +
    ggplot2::scale_colour_viridis_d(
      option = "turbo",
      end = 0.95,
      na.value = "grey80"
    ) +
    ggplot2::labs(
      x = "UMAP 1",
      y = "UMAP 2",
      colour = "Annotation"
    ) +
    .get_plot_theme(theme_name) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.5
      ),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      legend.position = "bottom"
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(
          size = 4,
          alpha = 1
        ),
        nrow = 1
      )
    )


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
