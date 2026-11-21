
test_that("extract_svg_gene_set extracts correct statistics", {

  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(
      counts = matrix(
        1:12,
        nrow = 3,
        dimnames = list(
          c("GeneA", "GeneB", "GeneC"),
          paste0("Cell", 1:4)
        )
      )
    )
  )

  SummarizedExperiment::rowData(spe)$LR_stat <- c(2.5, 4.0, 1.2)

  svg_norm <- list(
    None = spe,
    scran = spe
  )

  result <- extract_svg_gene_set(
    svg_norm = svg_norm,
    genes = c("GeneA", "GeneB")
  )

  expect_s3_class(result, "data.frame")

  expect_named(
    result,
    c("gene", "svg_stat", "method")
  )

  expect_equal(nrow(result), 4L)

  expect_setequal(
    unique(result$gene),
    c("GeneA", "GeneB")
  )

  expect_setequal(
    unique(result$method),
    c("None", "scran")
  )

  expect_equal(
    result$svg_stat,
    c(2.5, 4.0, 2.5, 4.0)
  )
})


test_that("extract_svg_gene_set handles missing genes", {

  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(
      counts = matrix(
        1:6,
        nrow = 2,
        dimnames = list(
          c("GeneA", "GeneB"),
          paste0("Cell", 1:3)
        )
      )
    )
  )

  SummarizedExperiment::rowData(spe)$LR_stat <- c(2, 3)

  result <- extract_svg_gene_set(
    svg_norm = list(None = spe),
    genes = "UnknownGene"
  )

  expect_equal(nrow(result), 0L)
})


test_that("extract_svg_gene_set reports missing statistics", {

  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(
      counts = matrix(
        1:6,
        nrow = 2,
        dimnames = list(
          c("GeneA", "GeneB"),
          paste0("Cell", 1:3)
        )
      )
    )
  )

  expect_error(
    extract_svg_gene_set(
      svg_norm = list(None = spe),
      genes = "GeneA",
      stat = "MissingStat"
    ),
    "not found for method"
  )
})


test_that("extract_svg_gene_set respects method order", {

  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(
      counts = matrix(
        1:6,
        nrow = 2,
        dimnames = list(
          c("GeneA", "GeneB"),
          paste0("Cell", 1:3)
        )
      )
    )
  )

  SummarizedExperiment::rowData(spe)$LR_stat <- c(2, 3)

  result <- extract_svg_gene_set(
    svg_norm = list(
      None = spe,
      scran = spe
    ),
    genes = "GeneA",
    method_order = c("scran", "None")
  )

  expect_equal(
    result$method,
    c("scran", "None")
  )
})
