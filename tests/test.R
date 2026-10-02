# # ------------------------------------------------------------
# # Normalisation
# # ------------------------------------------------------------
# library(devtools)
# devtools::load_all()
# norm <- normalize_spe(
#   spe,
#   methods = c(
#     "none",
#     "scran"
#   )
# )
#
# names(norm)
#
#
# # ------------------------------------------------------------
# # SVG analysis
# # ------------------------------------------------------------
#
# svg <- run_svg(
#   norm,
#   method = "nnSVG",
#   verbose = TRUE
# )
#
# names(svg)
#
# genes_use <-rownames(spe)[1:5]
#
# svg_stat <- plot_svg_stat(
#   svg,
#   genes = genes_use,
#   stat = "LR_stat",
#   transform = "log"
# )
#
# svg_stat$plot
#
#
# svg_rank <- plot_svg_rank(
#   svg,
#   genes = genes_use
# )
#
# svg_rank$plot
#
#
# p_scaled <- plot_spatial_expression(
#   norm_list = norm,
#   genes = genes_use,
#   scale_expr = TRUE,
#   point_size = 0.5
# )
#
# p_scaled$plot
#
#
#
# p_original <- plot_spatial_expression(
#   norm_list = norm,
#   genes = genes_use,
#   scale_expr = FALSE,
#   point_size = 0.5,
#   ncol_genes = NULL
# )
#
# p_original$plot
#
# stable <- plot_svg_gene_set_distribution(
#   svg,
#   genes = Housekeeping_Genes$Gene.name,
#   stat = "LR_stat",
#   transform = "log",
#   plot_type = "boxplot"
# )
#
# stable$plot
# topn <- plot_svg_topn_gene_set(
#   svg,
#   genes = Housekeeping_Genes$Gene.name,
#   top_n = c(
#     10,
#     20,
#     50,
#     100
#   )
# )
#
# topn$plot
#
#
#
# annotation_results <- run_annotation(
#   norm_list = norm,
#   ref = ref,
#   labels = ref$cell_type_spec,
#   de.method = "wilcox"
# )
# #
# #
# #
# # p1 <- plot_spatial_annotation(
# #   annotation_list = annotation_results,
# #   label_col = "label",
# #   ncol = 2
# # )
# #
# # p1$plot
# #
# # p2 <- plot_celltype_proportion(
# #   annotation_list = annotation_results,
# #   label_col = "label"
# # )
# #
# # p2$plot
#
# p_agree <- plot_annotation_agreement(
#   norm_list = annotation_results,
#   annotation_col = "SingleR_label"
# )
#
# p_agree$plot
















# ============================================================
# Example workflow
# ============================================================


# Run BANKSY
banksy_norm <- run_banksy_individual(
  norm_list = norm,
  cluster_col = "banksy_cluster",
  k_geom = 15,
  lambda = 0.8,
  resolution = 0.2,
  seed = 123
)


# Check stored clusters
table(
  banksy_norm$none$banksy_cluster
)

table(
  banksy_norm$scran$banksy_cluster
)


# Calculate ARI
ari_results <- calculate_ari(
  norm_list = banksy_norm,
  ground_truth_col = "layer_guess",
  cluster_col = "banksy_cluster"
)


ari_results


# ARI plot
ari_plot <- plot_ari(
  ari_df = ari_results,
  method_order = c(
    "none",
    "scran"
  )
)


ari_plot$plot


# Spatial clustering plot
cluster_plot <- plot_spatial_clustering(
  norm_list = banksy_norm,
  cluster_col = "banksy_cluster",
  method_order = c(
    "none",
    "scran"
  ),
  point_size = 0.5
)


cluster_plot$plot








# ============================================================
go_results <- calculate_go(
  norm_list = banksy_norm,
  cluster_col = "banksy_cluster",
  methods = c(
    "none",
    "scran"
  ),
  OrgDb = org.Hs.eg.db::org.Hs.eg.db
)

plot_go_genes(
  go_results,
  top_n = 5,
  method_order = c(
    "none",
    "scran"
  )
)$plot

plot_go(
  go_results,
  top_n = 5,
  method_order = c(
    "none",
    "scran"
  )
)$plot
