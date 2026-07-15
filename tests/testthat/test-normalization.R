test_that("PTM_normalization keeps Index aligned with its own site's values", {
  # PTM sites are deliberately interleaved across proteins in row order
  # (PROTA, PROTB, PROTA) so that grouping by protein during normalization
  # reorders rows relative to the original rowData order. This reproduces
  # the bug fixed in 0a0b322, where residuals ended up mislabeled with the
  # wrong Index (e.g. PROTB_S20 receiving PROTA_S30's residuals).
  samples <- paste0("S", 1:4)

  prot_assay <- matrix(
    c(1, 2, 3, 4,
      4, 3, 2, 1),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("PROTA", "PROTB"), samples)
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = prot_assay),
    rowData = S4Vectors::DataFrame(ProteinID = c("PROTA", "PROTB"), Index = c("PROTA", "PROTB")),
    colData = S4Vectors::DataFrame(sample_name = samples)
  )

  ptm_index <- c("PROTA_S10", "PROTB_S20", "PROTA_S30")
  ptm_prot <- c("PROTA", "PROTB", "PROTA")
  ptm_assay <- matrix(
    c(10, 20, 30, 40,
      1, 3, 5, 7,
      100, 90, 80, 70),
    nrow = 3, byrow = TRUE, dimnames = list(ptm_index, samples)
  )
  ptm_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = ptm_assay),
    rowData = S4Vectors::DataFrame(ProteinID = ptm_prot, Index = ptm_index),
    colData = S4Vectors::DataFrame(sample_name = samples)
  )

  result <- PTM_normalization(ptm_se, se)

  expect_setequal(rownames(result), ptm_index)

  expected <- list(
    PROTA_S10 = c(-29, -55 / 3, -23 / 3, 3),
    PROTA_S30 = c(61, 155 / 3, 127 / 3, 33),
    PROTB_S20 = c(-36, -104 / 3, -100 / 3, -32)
  )

  result_assay <- as.matrix(SummarizedExperiment::assay(result))
  for (site in names(expected)) {
    expect_equal(unname(result_assay[site, ]), expected[[site]], tolerance = 1e-6)
  }
})

test_that("PTM_normalization removes correlation with protein-level abundance", {
  samples <- paste0("S", 1:6)

  prot_assay <- matrix(
    c(1, 2, 3, 4, 5, 6,
      6, 5, 4, 3, 2, 1),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("PROTA", "PROTB"), samples)
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = prot_assay),
    rowData = S4Vectors::DataFrame(ProteinID = c("PROTA", "PROTB"), Index = c("PROTA", "PROTB")),
    colData = S4Vectors::DataFrame(sample_name = samples)
  )

  ptm_index <- c("PROTA_S1", "PROTB_S2")
  ptm_assay <- matrix(
    c(2, 5, 6, 9, 8, 13,
      13, 8, 9, 6, 5, 2),
    nrow = 2, byrow = TRUE, dimnames = list(ptm_index, samples)
  )
  ptm_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = ptm_assay),
    rowData = S4Vectors::DataFrame(ProteinID = c("PROTA", "PROTB"), Index = ptm_index),
    colData = S4Vectors::DataFrame(sample_name = samples)
  )

  result <- PTM_normalization(ptm_se, se)

  residual_mat <- as.matrix(SummarizedExperiment::assay(result))
  prot_ids <- SummarizedExperiment::rowData(result)$ProteinID
  prot_mat <- prot_assay[prot_ids, colnames(residual_mat), drop = FALSE]

  residuals <- as.numeric(residual_mat)
  prot_values <- as.numeric(prot_mat)

  expect_equal(sum(residuals), 0, tolerance = 1e-8)
  expect_equal(unname(cor(residuals, prot_values)), 0, tolerance = 1e-8)
})
