test_that("load data from RData works", {
  data <- readResultRData(test_path("testdata", "ccrcc.rda"))
  expect_equal(dim(data)[1], 10671)
  expect_equal(dim(data)[2], 24)
})

test_that("load experiment annotation file works", {
  data <- readExpDesign(test_path("testdata", "3plex_ccRCC_TMT", "experiment_annotation.tsv"), type = "TMT")
  expect_equal(dim(data)[1], 24)
})

test_that("load quant table works", {
  data <- readQuantTable(test_path("testdata", "3plex_ccRCC_TMT", "ccRCC_prot_abundance_MD_3plex.tsv"), type = "TMT")
  expect_equal(dim(data)[2], 29)
})

test_that("load data from files works", {
  data <- make_se_from_files(test_path("testdata", "3plex_ccRCC_TMT", "ccRCC_prot_abundance_MD_3plex.tsv"),
    test_path("testdata", "3plex_ccRCC_TMT", "experiment_annotation.tsv"),
    type = "TMT", level = "peptide"
  )
  expect_s4_class(data, "SummarizedExperiment")
  expect_equal(dim(data)[1], 10671)
  expect_equal(dim(data)[2], 24)
})

test_that("load data from files works (nccrcc plexes in ccRCC cohort)", {
  data <- make_se_from_files(test_path("testdata", "3plex_panRCC_TMT", "abundance_gene_MD.tsv"),
    test_path("testdata", "3plex_panRCC_TMT", "experiment_annotation.tsv"),
    type = "TMT", level = "peptide"
  )
  expect_s4_class(data, "SummarizedExperiment")
  expect_equal(dim(data)[1], 10608)
  expect_equal(dim(data)[2], 27)
})

test_that("load data from files works (nccrcc plexes in ccRCC cohort, peptide level)", {
  data <- make_se_from_files(test_path("testdata", "3plex_panRCC_TMT", "abundance_peptide_MD.tsv"),
    test_path("testdata", "3plex_panRCC_TMT", "experiment_annotation.tsv"),
    type = "TMT", level = "peptide"
  )
  expect_s4_class(data, "SummarizedExperiment")
  expect_equal(dim(data)[1], 153408)
  expect_equal(dim(data)[2], 27)

  cols <- colnames(rowData(data))
  expect_true("Gene" %in% cols)
  expect_true("Peptide" %in% cols)
  expect_true(all(!is.na(rowData(data)["Peptide"])))
})
