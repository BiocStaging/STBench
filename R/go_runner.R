# GO analysis from cluster-specific DE genes


#' Calculate GO enrichment from cluster-specific DE genes
#'
#' Identify cluster-specific differentially expressed genes
#' across normalisation methods and perform Gene Ontology
#' enrichment for each cluster.
#'
#' Each normalisation method should contain one
#' SpatialExperiment object representing the same spatial sample.
#'
#' @param norm_list Named list of normalised SpatialExperiment
#'   objects.
#' @param cluster_col Name of the colData column containing
#'   cluster assignments.
#' @param methods Optional character vector specifying the
#'   normalisation methods to analyse.
#' @param assay_name Name of the assay used for differential
#'   expression analysis.
#' @param OrgDb Organism annotation database passed to
#'   `clusterProfiler::enrichGO()`.
#' @param key_type Gene identifier type.
#' @param ont Gene Ontology category. One of `"BP"`, `"MF"`,
#'   `"CC"`, or `"ALL"`.
#' @param fdr_cutoff FDR threshold used to define significant
#'   cluster-specific DE genes.
#' @param min_prop Minimum proportion passed to
#'   `scran::findMarkers()`.
#' @param min_spots Minimum number of spatial locations required
#'   for a cluster.
#' @param min_genes Minimum number of DE genes required to run
#'   GO enrichment.
#' @param pvalue_cutoff GO enrichment p-value cutoff.
#' @param qvalue_cutoff GO enrichment q-value cutoff.
#' @param p_adjust_method Multiple-testing correction method.
#' @param readable Logical indicating whether gene IDs should
#'   be converted to readable symbols.
#'
#' @param test_type Statistical test used for differential expression
#'   analysis in `scran::findMarkers()`. Default is `"wilcox"`.
#' @return A list containing:
#' \itemize{
#'   \item `markers`: full marker results for each method.
#'   \item `genes`: significant DE genes used for GO.
#'   \item `go`: GO enrichment results.
#'   \item `universe`: background genes for each method.
#' }
#'
#' @examples
#' data(norm_ccc_example)
#'
#' go_result <- calculate_go(
#'   norm_list = norm_ccc_example,
#'   cluster_col = "layer_guess",
#'   methods = names(norm_ccc_example),
#'   assay_name = "logcounts",
#'   OrgDb = org.Mm.eg.db::org.Mm.eg.db,
#'   key_type = "SYMBOL",
#'   ont = "BP",
#'   fdr_cutoff = 0.05,
#'   min_prop = 0.4,
#'   min_spots = 10,
#'   min_genes = 20,
#'   pvalue_cutoff = 0.05,
#'   qvalue_cutoff = 0.05,
#'   p_adjust_method = "BH",
#'   readable = TRUE,
#'   test_type = "t"
#' )
#' @export
calculate_go <- function(
  norm_list,
  cluster_col = "banksy_cluster",
  methods = names(norm_list),
  assay_name = "logcounts",
  OrgDb,
  key_type = "SYMBOL",
  ont = "BP",
  fdr_cutoff = 0.05,
  min_prop = 0.7,
  min_spots = 30,
  min_genes = 10,
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.05,
  p_adjust_method = "BH",
  test_type = "wilcox",
  readable = TRUE
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


  # Initialise output


  marker_results <- list()

  gene_list <- list()

  go_list <- list()

  universe_list <- list()


  # Process methods


  for (method in methods) {
    spe <- norm_list[[method]]


    # Check SPE


    if (
      !methods::is(
        spe,
        "SpatialExperiment"
      )
    ) {
      stop(
        "Each method must contain one SpatialExperiment object. ",
        "Problem detected for method: ",
        method
      )
    }


    # Check assay


    if (
      !assay_name %in%
        SummarizedExperiment::assayNames(
          spe
        )
    ) {
      stop(
        "Assay '",
        assay_name,
        "' not found for method: ",
        method
      )
    }


    # Check cluster column


    cd <- SummarizedExperiment::colData(
      spe
    )


    if (
      !cluster_col %in%
        colnames(cd)
    ) {
      stop(
        "Cluster column '",
        cluster_col,
        "' not found for method: ",
        method
      )
    }


    # Extract clusters


    clusters <- factor(
      cd[[cluster_col]]
    )


    # Remove missing cluster labels


    keep <- !is.na(
      clusters
    )


    spe_use <- spe[
      ,
      keep,
      drop = FALSE
    ]


    clusters <- droplevels(
      clusters[
        keep
      ]
    )


    # Filter small clusters


    cluster_sizes <- table(
      clusters
    )


    valid_clusters <- names(
      cluster_sizes[
        cluster_sizes >= min_spots
      ]
    )


    keep <- clusters %in%
      valid_clusters


    spe_use <- spe_use[
      ,
      keep,
      drop = FALSE
    ]


    clusters <- droplevels(
      clusters[
        keep
      ]
    )


    # Check number of clusters


    if (
      length(
        levels(
          clusters
        )
      ) < 2
    ) {
      warning(
        "Fewer than two valid clusters for method: ",
        method,
        call. = FALSE
      )

      next
    }


    # Store gene universe


    universe_list[[method]] <- rownames(
      spe_use
    )


    # Find cluster-specific DE genes


    markers <- scran::findMarkers(
      spe_use,
      groups = clusters,
      assay.type = assay_name,
      direction = "up",
      pval.type = "some",
      min.prop = min_prop,
      test.type = "wilcox"
    )


    marker_results[[method]] <- markers


    # Process clusters


    for (cluster_id in names(markers)) {
      marker_df <- as.data.frame(
        markers[[cluster_id]]
      )


      marker_df$gene <- rownames(
        marker_df
      )


      # Significant DE genes


      marker_df <- marker_df[
        !is.na(
          marker_df$FDR
        ) &
          marker_df$FDR < fdr_cutoff, ,
        drop = FALSE
      ]


      marker_df <- marker_df[
        order(
          marker_df$FDR
        ), ,
        drop = FALSE
      ]


      # Skip empty results


      if (nrow(marker_df) == 0) {
        next
      }


      # Store genes


      gene_df <- data.frame(
        method = method,
        cluster = as.character(
          cluster_id
        ),
        gene = marker_df$gene,
        FDR = marker_df$FDR,
        rank = seq_len(
          nrow(marker_df)
        ),
        stringsAsFactors = FALSE
      )


      if (
        "summary.logFC" %in%
          colnames(marker_df)
      ) {
        gene_df$logFC <-
          marker_df$summary.logFC
      }


      gene_list[[
        paste(
          method,
          cluster_id,
          sep = "_"
        )
      ]] <- gene_df


      # Minimum genes for GO


      genes <- marker_df$gene


      if (length(genes) < min_genes) {
        warning(
          "Skipping GO for method '",
          method,
          "', cluster '",
          cluster_id,
          "': fewer than ",
          min_genes,
          " DE genes.",
          call. = FALSE
        )

        next
      }


      # GO enrichment


      go_result <- clusterProfiler::enrichGO(
        gene = genes,
        universe = universe_list[[method]],
        OrgDb = OrgDb,
        keyType = key_type,
        ont = ont,
        pAdjustMethod = p_adjust_method,
        pvalueCutoff = pvalue_cutoff,
        qvalueCutoff = qvalue_cutoff,
        readable = readable
      )


      go_df <- as.data.frame(
        go_result
      )


      if (nrow(go_df) == 0) {
        next
      }


      # Add identifiers


      go_df$method <- method

      go_df$cluster <- as.character(
        cluster_id
      )


      # Numeric gene ratio


      ratio_parts <- strsplit(
        go_df$GeneRatio,
        "/",
        fixed = TRUE
      )


      go_df$GeneRatio_numeric <- vapply(
        ratio_parts,
        function(x) {
          as.numeric(
            x[1]
          ) /
            as.numeric(
              x[2]
            )
        },
        numeric(1)
      )


      # -log10 adjusted p-value


      go_df$log10_padj <- -log10(
        go_df$p.adjust
      )


      # Store GO results


      go_list[[
        paste(
          method,
          cluster_id,
          sep = "_"
        )
      ]] <- go_df
    }
  }


  # Combine DE genes


  if (length(gene_list) > 0) {
    genes <- do.call(
      rbind,
      gene_list
    )

    rownames(genes) <- NULL
  } else {
    genes <- data.frame()
  }


  # Combine GO results


  if (length(go_list) > 0) {
    go <- do.call(
      rbind,
      go_list
    )

    rownames(go) <- NULL
  } else {
    go <- data.frame()

    warning(
      "No significant GO enrichment results were identified.",
      call. = FALSE
    )
  }


  # Return


  list(
    markers = marker_results,
    genes = genes,
    go = go,
    universe = universe_list
  )
}


# Plot GO enrichment


#' Plot GO enrichment results
#'
#' Plot enriched Gene Ontology terms across normalisation
#' methods and spatial clusters.
#'
#' Point size represents gene count, colour represents
#' enrichment significance, and shape represents the
#' normalisation method.
#'
#' @param go_result Output from calculate_go() or GO data frame.
#' @param top_n Number of top GO terms per method and cluster.
#'   If NULL, all available GO terms are displayed.
#' @param method_col Column containing normalisation methods.
#' @param cluster_col Column containing cluster assignments.
#' @param method_order Optional normalisation method order.
#' @param methods Optional methods to display.
#' @param title Optional plot title.
#' @param subtitle Optional plot subtitle.
#' @param point_size_range Minimum and maximum point sizes.
#' @param ncol Number of facet columns.
#' @param theme_name Plot theme.
#'
#' @return A list containing plotting data, full GO data,
#'   and the plot.
#'
#' @examples
#' data(norm_ccc_example)
#'
#' go_result <- calculate_go(
#'   norm_list = norm_ccc_example,
#'   cluster_col = "layer_guess",
#'   methods = names(norm_ccc_example),
#'   assay_name = "logcounts",
#'   OrgDb = org.Mm.eg.db::org.Mm.eg.db,
#'   key_type = "SYMBOL",
#'   ont = "BP",
#'   fdr_cutoff = 0.05,
#'   min_prop = 0.4,
#'   min_spots = 10,
#'   min_genes = 5,
#'   pvalue_cutoff = 0.05,
#'   qvalue_cutoff = 0.05,
#'   p_adjust_method = "BH",
#'   readable = TRUE,
#'   test_type = "t"
#' )
#'
#' p_go <- plot_go(
#'   go_result = go_result,
#'   top_n = 5,
#'   method_col = "method",
#'   cluster_col = "cluster",
#'   point_size_range = c(2, 7),
#'   theme_name = "classic",
#'   ncol = 2
#' )
#'
#' p_go$plot
#' @export
plot_go <- function(
  go_result,
  top_n = 5,
  method_col = "method",
  cluster_col = "cluster",
  method_order = NULL,
  point_size_range = c(2, 7),
  theme_name = "classic",
  methods = NULL,
  title = NULL,
  subtitle = NULL,
  ncol = 3
) {
  # Extract GO data


  if (
    is.list(go_result) &&
      "go" %in% names(go_result)
  ) {
    df <- go_result$go
  } else if (is.data.frame(go_result)) {
    df <- go_result
  } else {
    stop(
      "'go_result' must be the output from ",
      "'calculate_go()' or a GO data frame."
    )
  }


  # Check data


  if (nrow(df) == 0) {
    stop("No GO enrichment results are available.")
  }

  required_cols <- c(
    method_col,
    cluster_col,
    "Description",
    "Count",
    "p.adjust",
    "GeneRatio_numeric",
    "log10_padj"
  )

  missing_cols <- setdiff(
    required_cols,
    colnames(df)
  )

  if (length(missing_cols) > 0) {
    stop(
      "Missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }


  # Select methods


  available_methods <- unique(
    as.character(df[[method_col]])
  )

  if (is.null(methods)) {
    methods <- available_methods
  }

  unknown <- setdiff(
    methods,
    available_methods
  )

  if (length(unknown) > 0) {
    stop(
      "Unknown methods: ",
      paste(unknown, collapse = ", ")
    )
  }

  if (!is.null(method_order)) {
    methods <- c(
      intersect(method_order, methods),
      setdiff(methods, method_order)
    )
  }

  df <- df[
    df[[method_col]] %in% methods, ,
    drop = FALSE
  ]


  # Retain full selected GO results


  full_df <- df


  # Select top GO terms


  if (!is.null(top_n)) {
    if (
      length(top_n) != 1 ||
        !is.numeric(top_n) ||
        is.na(top_n) ||
        top_n < 1 ||
        top_n != floor(top_n)
    ) {
      stop("top_n must be a positive integer or NULL.")
    }

    split_df <- split(
      df,
      interaction(
        df[[method_col]],
        df[[cluster_col]],
        drop = TRUE
      )
    )

    plot_list <- lapply(
      split_df,
      function(x) {
        x <- x[
          order(x$p.adjust), ,
          drop = FALSE
        ]

        utils::head(
          x,
          top_n
        )
      }
    )

    plot_df <- do.call(
      rbind,
      plot_list
    )
  } else {
    plot_df <- df
  }

  rownames(plot_df) <- NULL


  # Factor ordering


  plot_df[[method_col]] <- factor(
    plot_df[[method_col]],
    levels = methods
  )

  plot_df[[cluster_col]] <- factor(
    plot_df[[cluster_col]]
  )


  # Plot


  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$GeneRatio_numeric,
      y = .data$Description
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = .data$Count,
        colour = .data$log10_padj,
        shape = .data[[method_col]]
      )
    ) +
    ggplot2::facet_wrap(
      stats::as.formula(
        paste("~", cluster_col)
      ),
      scales = "free_y",
      ncol = ncol
    ) +
    ggplot2::scale_size_continuous(
      range = point_size_range
    ) +
    ggplot2::scale_colour_viridis_c() +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Gene ratio",
      y = "GO term",
      size = "Gene count",
      shape = "Normalisation method",
      colour = expression(
        -log[10]("Adjusted p-value")
      )
    ) +
    .get_plot_theme(theme_name) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        size = 9
      ),
      axis.text.y = ggplot2::element_text(
        size = 8
      ),
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      legend.position = "right"
    )


  # Return


  return(
    list(
      data = plot_df,
      full_data = full_df,
      plot = p
    )
  )
}


# Plot DE genes used for GO


#' Plot cluster-specific DE genes used for GO
#'
#' Plot the top cluster-specific DE genes detected under each
#' normalisation method.
#'
#' @param go_result Result returned by `calculate_go()`, or the
#'   gene data frame contained in `go_result$genes`.
#' @param top_n Number of top DE genes shown per normalisation
#'   method and cluster.
#' @param method_col Name of the column containing normalisation
#'   methods.
#' @param cluster_col Name of the column containing clusters.
#' @param gene_col Name of the column containing gene names.
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
#' @examples
#' data(norm_ccc_example)
#'
#' go_result <- calculate_go(
#'   norm_list = norm_ccc_example,
#'   cluster_col = "layer_guess",
#'   methods = names(norm_ccc_example),
#'   assay_name = "logcounts",
#'   OrgDb = org.Mm.eg.db::org.Mm.eg.db,
#'   key_type = "SYMBOL",
#'   ont = "BP",
#'   fdr_cutoff = 0.05,
#'   min_prop = 0.4,
#'   min_spots = 10,
#'   min_genes = 5,
#'   pvalue_cutoff = 0.05,
#'   qvalue_cutoff = 0.05,
#'   p_adjust_method = "BH",
#'   readable = TRUE,
#'   test_type = "t"
#' )
#'
#' if (nrow(go_result$genes) > 0L) {
#'   p_go_genes <- plot_go_genes(
#'     go_result = go_result,
#'     top_n = 5,
#'     method_col = "method",
#'     cluster_col = "cluster",
#'     gene_col = "gene",
#'     point_size = 3,
#'     theme_name = "classic"
#'   )
#'
#'   p_go_genes$plot
#' }
#' @export
plot_go_genes <- function(
  go_result,
  top_n = 5,
  method_col = "method",
  cluster_col = "cluster",
  gene_col = "gene",
  method_order = NULL,
  point_size = 3,
  palette = NULL,
  theme_name = "classic"
) {
  # Extract gene data


  if (
    is.list(go_result) &&
      "genes" %in%
        names(go_result)
  ) {
    df <- go_result$genes
  } else if (is.data.frame(go_result)) {
    df <- go_result
  } else {
    stop(
      "'go_result' must be the output from 'calculate_go()' ",
      "or a DE gene data frame."
    )
  }


  # Check data


  if (nrow(df) == 0) {
    stop(
      "No DE genes are available to plot."
    )
  }


  required_cols <- c(
    method_col,
    cluster_col,
    gene_col,
    "FDR"
  )


  missing_cols <- setdiff(
    required_cols,
    colnames(df)
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


  # Resolve methods


  available_methods <- unique(
    as.character(
      df[[method_col]]
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


  df <- df[
    df[[method_col]] %in%
      methods, ,
    drop = FALSE
  ]


  # Select top N DE genes


  split_df <- split(
    df,
    interaction(
      df[[method_col]],
      df[[cluster_col]],
      drop = TRUE
    )
  )


  plot_list <- lapply(
    split_df,
    function(x) {
      x <- x[
        order(
          x$FDR
        ), ,
        drop = FALSE
      ]


      utils::head(
        x,
        top_n
      )
    }
  )


  plot_df <- do.call(
    rbind,
    plot_list
  )


  rownames(plot_df) <- NULL


  # Factor ordering


  plot_df[[method_col]] <- factor(
    plot_df[[method_col]],
    levels = methods
  )


  plot_df[[cluster_col]] <- factor(
    plot_df[[cluster_col]]
  )


  # Plot


  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data[[method_col]],
      y = .data[[gene_col]],
      colour = .data[[method_col]]
    )
  ) +
    ggplot2::geom_point(
      size = point_size
    ) +
    ggplot2::facet_wrap(
      stats::as.formula(
        paste(
          "~",
          cluster_col
        )
      ),
      scales = "free_y"
    ) +
    ggplot2::labs(
      x = "Normalisation method",
      y = NULL,
      colour = "Method"
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
      axis.text.x =
        ggplot2::element_text(
          angle = 45,
          hjust = 1
        ),
      axis.text.y =
        ggplot2::element_text(
          size = 8
        ),
      legend.position = "right"
    )


  # Optional palette


  if (!is.null(palette)) {
    p <- p +

      ggplot2::scale_colour_manual(
        values = palette
      )
  }


  # Return


  list(
    data = plot_df,
    plot = p
  )
}
