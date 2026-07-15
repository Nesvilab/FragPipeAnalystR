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
    c(1, 2, 3, 4, 5, 6),
    nrow = 1, dimnames = list("PROTA", samples)
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = prot_assay),
    rowData = S4Vectors::DataFrame(ProteinID = "PROTA", Index = "PROTA"),
    colData = S4Vectors::DataFrame(sample_name = samples)
  )

  ptm_assay <- matrix(
    c(2, 5, 6, 9, 8, 13),
    nrow = 1, dimnames = list("PROTA_S1", samples)
  )
  ptm_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = ptm_assay),
    rowData = S4Vectors::DataFrame(ProteinID = "PROTA", Index = "PROTA_S1"),
    colData = S4Vectors::DataFrame(sample_name = samples)
  )

  result <- PTM_normalization(ptm_se, se)

  residuals <- as.numeric(SummarizedExperiment::assay(result)["PROTA_S1", ])
  expect_equal(sum(residuals), 0, tolerance = 1e-8)
  expect_equal(unname(cor(residuals, as.numeric(prot_assay))), 0, tolerance = 1e-8)
})
