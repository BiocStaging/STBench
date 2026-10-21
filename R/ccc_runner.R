# ============================================================
# CCC RUNNER
# ============================================================


#' Run CellChat across normalization methods
#'
#' Runs CellChat cell-cell communication analysis across multiple
#' normalized versions of the same SpatialExperiment object.
#'
#' @param spelist Named list of SpatialExperiment objects.
#'   Names should correspond to normalization methods.
#' @param celltype_col Character. Column in colData(spe)
#'   containing cell-type labels.
#' @param database CellChat database object, e.g. CellChatDB.mouse.
#' @param spatial_factors Data frame containing CellChat spatial
#'   factors, typically ratio and tol.
#' @param interaction_range Numeric. Maximum spatial interaction distance.
#' @param contact_range Numeric. Contact-dependent interaction distance.
#' @param assay_name Character. Assay containing normalized expression.
#'   Default is "logcounts".
#' @param nboot Integer. Number of bootstrap iterations.
#' @param min_cells Integer. Minimum number of spots/cells per group.
#' @param trim Numeric. Trimming proportion for truncated mean.
#' @param scale_distance Numeric. Distance scaling used by CellChat.
#' @param variable_both Logical. Passed to
#'   identifyOverExpressedInteractions().
#' @param verbose Logical. Print progress messages.
#'
#' @return Named list containing CellChat results for each method.
#'
#' @export
run_ccc <- function(
    spelist,
    celltype_col,
    database,
    spatial_factors,
    interaction_range,
    contact_range,
    assay_name = "logcounts",
    nboot = 100,
    min_cells = 10,
    trim = 0.1,
    scale_distance = 0.1,
    variable_both = FALSE,
    verbose = TRUE
) {
  
  # ----------------------------------------------------------
  # Dependency
  # ----------------------------------------------------------
  
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    stop(
      "Package 'CellChat' is required for run_ccc().",
      call. = FALSE
    )
  }
  
  
  # ----------------------------------------------------------
  # Input checks
  # ----------------------------------------------------------
  
  if (!is.list(spelist) || length(spelist) == 0) {
    stop(
      "'spelist' must be a non-empty named list.",
      call. = FALSE
    )
  }
  
  if (is.null(names(spelist)) || any(names(spelist) == "")) {
    stop(
      "'spelist' must be named by normalization method.",
      call. = FALSE
    )
  }
  
  if (!is.character(celltype_col) ||
      length(celltype_col) != 1) {
    
    stop(
      "'celltype_col' must be a single column name.",
      call. = FALSE
    )
  }
  
  if (missing(database) || is.null(database)) {
    stop(
      "'database' must be supplied.",
      call. = FALSE
    )
  }
  
  if (!is.data.frame(spatial_factors)) {
    stop(
      "'spatial_factors' must be a data.frame.",
      call. = FALSE
    )
  }
  
  if (!all(c("ratio", "tol") %in% names(spatial_factors))) {
    stop(
      "'spatial_factors' must contain columns 'ratio' and 'tol'.",
      call. = FALSE
    )
  }
  
  
  # ----------------------------------------------------------
  # Run all normalization methods
  # ----------------------------------------------------------
  
  results <- lapply(
    names(spelist),
    function(method_name) {
      
      if (verbose) {
        message(
          "Running CellChat | method: ",
          method_name
        )
      }
      
      .run_ccc_one(
        spe = spelist[[method_name]],
        method_name = method_name,
        celltype_col = celltype_col,
        database = database,
        spatial_factors = spatial_factors,
        interaction_range = interaction_range,
        contact_range = contact_range,
        assay_name = assay_name,
        nboot = nboot,
        min_cells = min_cells,
        trim = trim,
        scale_distance = scale_distance,
        variable_both = variable_both
      )
    }
  )
  
  names(results) <- names(spelist)
  
  results
}



# ============================================================
# INTERNAL: RUN CCC FOR ONE SPE
# ============================================================

.run_ccc_one <- function(
    spe,
    method_name,
    celltype_col,
    database,
    spatial_factors,
    interaction_range,
    contact_range,
    assay_name,
    nboot,
    min_cells,
    trim,
    scale_distance,
    variable_both
) {
  
  # ----------------------------------------------------------
  # Metadata
  # ----------------------------------------------------------
  
  meta <- as.data.frame(
    SummarizedExperiment::colData(spe)
  )
  
  if (!celltype_col %in% colnames(meta)) {
    stop(
      "Column '",
      celltype_col,
      "' was not found in colData for method '",
      method_name,
      "'.",
      call. = FALSE
    )
  }
  
  
  # ----------------------------------------------------------
  # Remove missing labels
  # ----------------------------------------------------------
  
  keep <- !is.na(meta[[celltype_col]]) &
    trimws(as.character(meta[[celltype_col]])) != ""
  
  spe <- spe[, keep, drop = FALSE]
  
  if (ncol(spe) == 0) {
    stop(
      "No labelled spots/cells remain for method '",
      method_name,
      "'.",
      call. = FALSE
    )
  }
  
  
  # ----------------------------------------------------------
  # Expression matrix
  # ----------------------------------------------------------
  
  if (!assay_name %in%
      SummarizedExperiment::assayNames(spe)) {
    
    stop(
      "Assay '",
      assay_name,
      "' was not found for method '",
      method_name,
      "'.",
      call. = FALSE
    )
  }
  
  data.input <- SummarizedExperiment::assay(
    spe,
    assay_name
  )
  
  if (!inherits(data.input, "dgCMatrix")) {
    
    data.input <- methods::as(
      as.matrix(data.input),
      "dgCMatrix"
    )
  }
  
  
  # ----------------------------------------------------------
  # Metadata for CellChat
  # ----------------------------------------------------------
  
  meta <- as.data.frame(
    SummarizedExperiment::colData(spe)
  )
  
  meta$labels <- droplevels(
    factor(
      meta[[celltype_col]]
    )
  )
  
  meta$samples <- factor(
    rep(
      "sample1",
      nrow(meta)
    )
  )
  
  rownames(meta) <- colnames(data.input)
  
  
  # ----------------------------------------------------------
  # Spatial coordinates
  # ----------------------------------------------------------
  
  coordinates <- as.data.frame(
    SpatialExperiment::spatialCoords(spe)
  )
  
  if (ncol(coordinates) < 2) {
    stop(
      "spatialCoords(spe) must contain at least two columns.",
      call. = FALSE
    )
  }
  
  coordinates <- coordinates[
    ,
    1:2,
    drop = FALSE
  ]
  
  colnames(coordinates) <- c(
    "x",
    "y"
  )
  
  rownames(coordinates) <- colnames(data.input)
  
  
  # ----------------------------------------------------------
  # Create CellChat object
  # ----------------------------------------------------------
  
  cellchat <- CellChat::createCellChat(
    object = data.input,
    meta = meta,
    group.by = "labels",
    datatype = "spatial",
    coordinates = coordinates,
    spatial.factors = spatial_factors
  )
  
  
  # ----------------------------------------------------------
  # Database
  # ----------------------------------------------------------
  
  cellchat@DB <- CellChat::subsetDB(
    database
  )
  
  
  # ----------------------------------------------------------
  # Preprocessing
  # ----------------------------------------------------------
  
  cellchat <- CellChat::subsetData(
    cellchat
  )
  
  cellchat <- CellChat::identifyOverExpressedGenes(
    cellchat
  )
  
  cellchat <- CellChat::identifyOverExpressedInteractions(
    cellchat,
    variable.both = variable_both
  )
  
  
  # ----------------------------------------------------------
  # Check LR pairs
  # ----------------------------------------------------------
  
  if (is.null(cellchat@LR$LRsig) ||
      nrow(cellchat@LR$LRsig) == 0) {
    
    stop(
      "No ligand-receptor pairs were identified for method '",
      method_name,
      "'.",
      call. = FALSE
    )
  }
  
  
  # ----------------------------------------------------------
  # Communication probability
  # ----------------------------------------------------------
  
  cellchat <- CellChat::computeCommunProb(
    cellchat,
    type = "truncatedMean",
    trim = trim,
    distance.use = TRUE,
    interaction.range = interaction_range,
    scale.distance = scale_distance,
    contact.dependent = TRUE,
    contact.range = contact_range,
    nboot = nboot
  )
  
  cellchat <- CellChat::filterCommunication(
    cellchat,
    min.cells = min_cells
  )
  
  
  # ----------------------------------------------------------
  # LR-level results
  # ----------------------------------------------------------
  
  df.net <- CellChat::subsetCommunication(
    cellchat
  )
  
  if (nrow(df.net) > 0) {
    df.net$method <- method_name
  }
  
  
  # ----------------------------------------------------------
  # Pathway-level results
  # ----------------------------------------------------------
  
  cellchat <- CellChat::computeCommunProbPathway(
    cellchat
  )
  
  df.netP <- CellChat::subsetCommunication(
    cellchat,
    slot.name = "netP"
  )
  
  if (nrow(df.netP) > 0) {
    df.netP$method <- method_name
  }
  
  
  # ----------------------------------------------------------
  # Aggregate network
  # ----------------------------------------------------------
  
  cellchat <- CellChat::aggregateNet(
    cellchat
  )
  
  
  # ----------------------------------------------------------
  # Return
  # ----------------------------------------------------------
  
  list(
    cellchat = cellchat,
    df.net = df.net,
    df.netP = df.netP
  )
}



# ============================================================
# INTERNAL: COMBINE CCC RESULTS
# ============================================================

.combine_ccc_results <- function(
    ccc_result,
    type = c(
      "pathway",
      "LR"
    )
) {
  
  type <- match.arg(type)
  
  if (!is.list(ccc_result) ||
      length(ccc_result) == 0) {
    
    stop(
      "'ccc_result' must be the output from run_ccc().",
      call. = FALSE
    )
  }
  
  
  element_name <- if (type == "pathway") {
    "df.netP"
  } else {
    "df.net"
  }
  
  
  out <- lapply(
    names(ccc_result),
    function(method_name) {
      
      x <- ccc_result[[method_name]][[element_name]]
      
      if (is.null(x) || nrow(x) == 0) {
        return(NULL)
      }
      
      if (!"method" %in% colnames(x)) {
        x$method <- method_name
      }
      
      x
    }
  )
  
  
  out <- Filter(
    Negate(is.null),
    out
  )
  
  
  if (length(out) == 0) {
    stop(
      "No ",
      type,
      " results were found.",
      call. = FALSE
    )
  }
  
  
  dplyr::bind_rows(out)
}



# ============================================================
# PLOT 1: PATHWAY STRENGTH
# ============================================================


#' Plot CCC signaling pathway strength
#'
#' @param ccc_result Output from run_ccc().
#' @param pathways Optional character vector of pathways.
#' @param top_n Number of strongest pathways when pathways = NULL.
#' @param method_order Optional method order.
#' @param method_cols Optional named vector of method colors.
#' @param log_scale Logical. Use log10 y-axis.
#'
#' @return A ggplot object.
#'
#' @export
plot_ccc_pathways <- function(
    ccc_result,
    pathways = NULL,
    top_n = 10,
    method_order = NULL,
    method_cols = NULL,
    log_scale = FALSE
) {
  
  df <- .combine_ccc_results(
    ccc_result,
    type = "pathway"
  )
  
  
  plot_df <- df %>%
    dplyr::group_by(
      method,
      pathway_name
    ) %>%
    dplyr::summarise(
      pathway_strength = sum(
        prob,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  
  # ----------------------------------------------------------
  # Select pathways
  # ----------------------------------------------------------
  
  if (is.null(pathways)) {
    
    pathways <- plot_df %>%
      dplyr::group_by(
        pathway_name
      ) %>%
      dplyr::summarise(
        mean_strength = mean(
          pathway_strength,
          na.rm = TRUE
        ),
        .groups = "drop"
      ) %>%
      dplyr::slice_max(
        mean_strength,
        n = top_n,
        with_ties = FALSE
      ) %>%
      dplyr::pull(
        pathway_name
      )
  }
  
  
  plot_df <- plot_df %>%
    dplyr::filter(
      pathway_name %in% pathways
    )
  
  
  # ----------------------------------------------------------
  # Method order
  # ----------------------------------------------------------
  
  if (is.null(method_order)) {
    method_order <- names(ccc_result)
  }
  
  plot_df$method <- factor(
    plot_df$method,
    levels = method_order
  )
  
  
  plot_df$pathway_name <- factor(
    plot_df$pathway_name,
    levels = rev(pathways)
  )
  
  
  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------
  
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = pathway_name,
      y = pathway_strength,
      fill = method
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(
        width = 0.8
      ),
      width = 0.7
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = "Signaling pathway",
      y = "Communication strength",
      fill = "Normalization method"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )
  
  
  if (!is.null(method_cols)) {
    
    p <- p +
      ggplot2::scale_fill_manual(
        values = method_cols
      )
  }
  
  
  if (log_scale) {
    
    p <- p +
      ggplot2::scale_y_log10()
  }
  
  
  p
}



# ============================================================
# PLOT 2: LIGAND-RECEPTOR STRENGTH
# ============================================================


#' Plot CCC ligand-receptor interactions
#'
#' @param ccc_result Output from run_ccc().
#' @param pathways Optional pathways to include.
#' @param interactions Optional LR interactions to include.
#' @param top_n Number of strongest LR interactions.
#' @param method_order Optional method order.
#' @param method_cols Optional named vector of method colors.
#' @param log_scale Logical. Use log10 y-axis.
#'
#' @return A ggplot object.
#'
#' @export
plot_ccc_lr <- function(
    ccc_result,
    pathways = NULL,
    interactions = NULL,
    top_n = 15,
    method_order = NULL,
    method_cols = NULL,
    log_scale = FALSE
) {
  
  df <- .combine_ccc_results(
    ccc_result,
    type = "LR"
  )
  
  
  # ----------------------------------------------------------
  # Pathway filter
  # ----------------------------------------------------------
  
  if (!is.null(pathways)) {
    
    df <- df %>%
      dplyr::filter(
        pathway_name %in% pathways
      )
  }
  
  
  # ----------------------------------------------------------
  # Summarise
  # ----------------------------------------------------------
  
  plot_df <- df %>%
    dplyr::group_by(
      method,
      pathway_name,
      interaction_name
    ) %>%
    dplyr::summarise(
      LR_strength = sum(
        prob,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  
  # ----------------------------------------------------------
  # Select strongest LR pairs
  # ----------------------------------------------------------
  
  if (is.null(interactions)) {
    
    interactions <- plot_df %>%
      dplyr::group_by(
        interaction_name
      ) %>%
      dplyr::summarise(
        mean_strength = mean(
          LR_strength,
          na.rm = TRUE
        ),
        .groups = "drop"
      ) %>%
      dplyr::slice_max(
        mean_strength,
        n = top_n,
        with_ties = FALSE
      ) %>%
      dplyr::pull(
        interaction_name
      )
  }
  
  
  plot_df <- plot_df %>%
    dplyr::filter(
      interaction_name %in% interactions
    )
  
  
  # ----------------------------------------------------------
  # Method order
  # ----------------------------------------------------------
  
  if (is.null(method_order)) {
    method_order <- names(ccc_result)
  }
  
  plot_df$method <- factor(
    plot_df$method,
    levels = method_order
  )
  
  
  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------
  
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = interaction_name,
      y = LR_strength,
      fill = method
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(
        width = 0.8
      ),
      width = 0.7
    ) +
    ggplot2::facet_wrap(
      ~ pathway_name,
      scales = "free_x"
    ) +
    ggplot2::labs(
      x = "Ligand-receptor pair",
      y = "Communication strength",
      fill = "Normalization method"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )
  
  
  if (!is.null(method_cols)) {
    
    p <- p +
      ggplot2::scale_fill_manual(
        values = method_cols
      )
  }
  
  
  if (log_scale) {
    
    p <- p +
      ggplot2::scale_y_log10()
  }
  
  
  p
}



# ============================================================
# PLOT 3: SENDER-RECEIVER COMMUNICATION
# ============================================================


#' Plot CCC sender-receiver communication
#'
#' @param ccc_result Output from run_ccc().
#' @param methods Optional methods to include.
#' @param method_order Optional method order.
#'
#' @return A ggplot object.
#'
#' @export
plot_ccc_sender_receiver <- function(
    ccc_result,
    methods = NULL,
    method_order = NULL
) {
  
  df <- .combine_ccc_results(
    ccc_result,
    type = "LR"
  )
  
  
  # ----------------------------------------------------------
  # Methods
  # ----------------------------------------------------------
  
  if (is.null(methods)) {
    methods <- names(ccc_result)
  }
  
  df <- df %>%
    dplyr::filter(
      method %in% methods
    )
  
  
  if (is.null(method_order)) {
    method_order <- methods
  }
  
  
  # ----------------------------------------------------------
  # Summarise
  # ----------------------------------------------------------
  
  plot_df <- df %>%
    dplyr::group_by(
      method,
      source,
      target
    ) %>%
    dplyr::summarise(
      communication_strength = sum(
        prob,
        na.rm = TRUE
      ),
      n_interactions = dplyr::n(),
      .groups = "drop"
    )
  
  
  plot_df$method <- factor(
    plot_df$method,
    levels = method_order
  )
  
  
  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------
  
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = source,
      y = target,
      size = communication_strength
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.8
    ) +
    ggplot2::facet_wrap(
      ~ method
    ) +
    ggplot2::labs(
      x = "Sender cell type",
      y = "Receiver cell type",
      size = "Communication\nstrength"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )
}



# ============================================================
# PLOT 4: CELLCHAT CIRCLE NETWORKS
# ============================================================


#' Plot CellChat circle networks
#'
#' Produces CellChat circle plots for one or more pathways
#' across normalization methods.
#'
#' Rows correspond to pathways and columns correspond to
#' normalization methods.
#'
#' @param ccc_result Output from run_ccc().
#' @param pathways Character vector of pathways.
#' @param methods Optional normalization methods.
#' @param method_order Optional method order.
#'
#' @return Invisibly returns NULL.
#'
#' @export
plot_ccc_circle <- function(
    ccc_result,
    pathways,
    methods = NULL,
    method_order = NULL
) {
  
  # ----------------------------------------------------------
  # Dependency
  # ----------------------------------------------------------
  
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    stop(
      "Package 'CellChat' is required for plot_ccc_circle().",
      call. = FALSE
    )
  }
  
  
  # ----------------------------------------------------------
  # Check result
  # ----------------------------------------------------------
  
  if (!is.list(ccc_result) ||
      length(ccc_result) == 0) {
    
    stop(
      "'ccc_result' must be the output from run_ccc().",
      call. = FALSE
    )
  }
  
  
  # ----------------------------------------------------------
  # Pathways
  # ----------------------------------------------------------
  
  if (missing(pathways) ||
      length(pathways) == 0) {
    
    stop(
      "'pathways' must contain at least one pathway.",
      call. = FALSE
    )
  }
  
  pathways <- as.character(pathways)
  
  
  # ----------------------------------------------------------
  # Methods
  # ----------------------------------------------------------
  
  if (is.null(methods)) {
    methods <- names(ccc_result)
  }
  
  methods <- intersect(
    methods,
    names(ccc_result)
  )
  
  
  if (length(methods) == 0) {
    
    stop(
      "No matching normalization methods were found.",
      call. = FALSE
    )
  }
  
  
  # ----------------------------------------------------------
  # Method order
  # ----------------------------------------------------------
  
  if (is.null(method_order)) {
    method_order <- names(ccc_result)
  }
  
  methods <- method_order[
    method_order %in% methods
  ]
  
  
  # ----------------------------------------------------------
  # Save graphical parameters
  # ----------------------------------------------------------
  
  old_par <- graphics::par(
    no.readonly = TRUE
  )
  
  on.exit(
    graphics::par(old_par),
    add = TRUE
  )
  
  
  # ----------------------------------------------------------
  # Layout
  #
  # rows    = pathways
  # columns = methods
  # ----------------------------------------------------------
  
  graphics::par(
    mfrow = c(
      length(pathways),
      length(methods)
    ),
    xpd = TRUE,
    mar = c(
      1,
      1,
      3,
      1
    )
  )
  
  
  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------
  
  for (pathway in pathways) {
    
    for (method_name in methods) {
      
      cellchat <- ccc_result[[method_name]]$cellchat
      
      
      # ------------------------------------------------------
      # Missing CellChat object
      # ------------------------------------------------------
      
      if (is.null(cellchat)) {
        
        graphics::plot.new()
        
        graphics::title(
          main = paste0(
            method_name,
            "\nNo CellChat result"
          )
        )
        
        next
      }
      
      
      # ------------------------------------------------------
      # Pathway not detected
      # ------------------------------------------------------
      
      if (!pathway %in% cellchat@netP$pathways) {
        
        graphics::plot.new()
        
        graphics::title(
          main = paste0(
            method_name,
            "\n",
            pathway,
            " not detected"
          )
        )
        
        next
      }
      
      
      # ------------------------------------------------------
      # Circle network
      # ------------------------------------------------------
      
      CellChat::netVisual_aggregate(
        cellchat,
        signaling = pathway,
        layout = "circle"
      )
      
      
      graphics::title(
        main = paste0(
          method_name,
          " | ",
          pathway
        )
      )
    }
  }
  
  
  invisible(NULL)
}