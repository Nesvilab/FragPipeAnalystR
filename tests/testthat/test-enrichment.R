# Helper: build a minimal site-level SE with SequenceWindow and limma columns
make_site_se <- function(seqs, fc_values, contrast = "A_vs_B") {
  n <- length(seqs)
  mat <- matrix(rnorm(n * 4), nrow = n,
                dimnames = list(seqs, paste0("s", 1:4)))
  rd <- S4Vectors::DataFrame(
    name           = seqs,
    ID             = seqs,
    SequenceWindow = seqs
  )
  rd[[paste0(contrast, "_diff")]]        <- fc_values
  rd[[paste0(contrast, "_significant")]] <- rep(TRUE, n)
  SummarizedExperiment::SummarizedExperiment(
    assays  = list(intensity = mat),
    rowData = rd
  )
}

test_that("run_kinase_zscore returns expected data frame structure", {
  seqs <- c("AAAGSSSGGGAAAA", "BBBGSSSGGBBBB", "CCCGSSSGGCCCCC",
            "DDDGSSSGGDDDDD", "EEEGSSSGGEEEE")
  dep <- make_site_se(seqs, fc_values = c(2, -1, 0.5, -0.5, 1.5))

  ks_lib <- data.frame(
    source   = c("KIN1", "KIN1", "KIN1", "KIN2", "KIN2", "KIN2"),
    sequence = c("AAAGSSSGGGAAAA", "BBBGSSSGGBBBB", "CCCGSSSGGCCCCC",
                 "CCCGSSSGGCCCCC", "DDDGSSSGGDDDDD", "EEEGSSSGGEEEE"),
    mor      = c(1, 1, 1, -1, 1, 1),
    target   = c("P1", "P2", "P3", "P3", "P4", "P5")
  )

  res <- run_kinase_zscore(dep, ks_lib, min_targets = 3L)

  expect_s3_class(res, "data.frame")
  expect_true(all(c("kinase", "n_substrates", "score", "p_value",
                    "contrast", "adj_p_value") %in% colnames(res)))
  expect_true(nrow(res) >= 1)
  expect_true(all(!is.na(res$score)))
  expect_true(all(res$p_value >= 0 & res$p_value <= 1))
  expect_true(all(res$adj_p_value >= 0 & res$adj_p_value <= 1))
})

test_that("run_kinase_zscore errors when SequenceWindow column is missing", {
  dep <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(intensity = matrix(1:4, 2, 2)),
    rowData = S4Vectors::DataFrame(name = c("A", "B"), ID = c("A", "B"))
  )
  ks_lib <- data.frame(source = "K", sequence = "A", mor = 1)
  expect_error(run_kinase_zscore(dep, ks_lib), "SequenceWindow")
})

test_that("run_kinase_zscore errors when no contrast columns exist", {
  seqs <- c("SEQ1", "SEQ2", "SEQ3")
  mat  <- matrix(1:9, 3, 3, dimnames = list(seqs, paste0("s", 1:3)))
  rd   <- S4Vectors::DataFrame(name = seqs, ID = seqs, SequenceWindow = seqs)
  dep  <- SummarizedExperiment::SummarizedExperiment(
    assays = list(intensity = mat), rowData = rd
  )
  ks_lib <- data.frame(source = "K", sequence = "SEQ1", mor = 1)
  expect_error(run_kinase_zscore(dep, ks_lib), "No contrast columns")
})

test_that("run_kinase_zscore errors when no library sequences match", {
  seqs <- c("AAAGSSSGGGAAAA", "BBBGSSSGGBBBB", "CCCGSSSGGCCCCC")
  dep  <- make_site_se(seqs, fc_values = c(1, -1, 0.5))
  ks_lib <- data.frame(
    source   = c("K1", "K1", "K1"),
    sequence = c("NOMATCHA", "NOMATCHB", "NOMATCHC"),
    mor      = c(1, 1, 1)
  )
  expect_error(run_kinase_zscore(dep, ks_lib), "No kinase-substrate library entries match")
})

test_that("run_kinase_zscore errors when no kinase meets min_targets threshold", {
  seqs <- c("AAAGSSSGGGAAAA", "BBBGSSSGGBBBB")
  dep  <- make_site_se(seqs, fc_values = c(1, -1))
  ks_lib <- data.frame(
    source   = c("K1", "K1"),
    sequence = c("AAAGSSSGGGAAAA", "BBBGSSSGGBBBB"),
    mor      = c(1, 1)
  )
  expect_error(run_kinase_zscore(dep, ks_lib, min_targets = 5L),
               "No kinases with")
})

test_that("visualize PTM-SEA result works", {
  expect_true(ggplot2::is.ggplot(visualize_PTMSEA(test_path("testdata", "ccRCC-combined.gct"),
                                                   "Tumor_vs_NAT_diff",
                                                   score_cutoff = 0, fdr_pvalue_cutoff = 1)))
  expect_true(ggplot2::is.ggplot(visualize_PTMSEA(test_path("testdata", "ccRCC-combined.gct"),
                   "Tumor_vs_NAT_diff",
                   selected_concepts = c("KINASE-PSP_PKCA/PRKCA", "KINASE-PSP_PKCB/PRKCB", "KINASE-PSP_PKCG/PRKCG", "KINASE-PSP_PKCB_iso2/PRKCB", "KINASE-PSP_PKCD/PRKCD", "KINASE-PSP_PKCI/PRKCI", "KINASE-PSP_PKCT/PRKCQ",
                                        "KINASE-PSP_PKCH/PRKCH", "KINASE-PSP_PKCE/PRKCE", "KINASE-PSP_PKCZ/PRKCZ"),
                   score_cutoff = 0, fdr_pvalue_cutoff = 1)))
})
