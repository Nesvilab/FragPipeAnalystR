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

test_that("readExpDesign stops when sample_name column is missing", {
  f <- tempfile(fileext = ".tsv")
  write.table(
    data.frame(plex = 1, channel = "126", sample = "S1", condition = "A", replicate = 1),
    f, sep = "\t", quote = FALSE, row.names = FALSE
  )
  expect_error(readExpDesign(f, type = "TMT"), "'sample_name' column is not present")
  expect_error(readExpDesign(f, type = "LFQ"), "'sample_name' column is not present")
})

test_that("readExpDesign stops when sample column is missing", {
  f <- tempfile(fileext = ".tsv")
  write.table(
    data.frame(plex = 1, channel = "126", sample_name = "S1", condition = "A", replicate = 1),
    f, sep = "\t", quote = FALSE, row.names = FALSE
  )
  expect_error(readExpDesign(f, type = "TMT"), "'sample' column is not present")
  expect_error(readExpDesign(f, type = "LFQ"), "'sample' column is not present")
})

test_that("readExpDesign stops on duplicated sample_name", {
  f <- tempfile(fileext = ".tsv")
  write.table(
    data.frame(plex = c(1, 1), channel = c("126", "127N"), sample = c("S1", "S2"),
               sample_name = c("dup", "dup"), condition = c("A", "A"), replicate = c(1, 2)),
    f, sep = "\t", quote = FALSE, row.names = FALSE
  )
  expect_error(readExpDesign(f, type = "TMT"), "Duplicated 'sample_name'")
  expect_error(readExpDesign(f, type = "LFQ"), "Duplicated 'sample_name'")
})

test_that("readExpDesign stops on duplicated sample", {
  f <- tempfile(fileext = ".tsv")
  write.table(
    data.frame(plex = c(1, 1), channel = c("126", "127N"), sample = c("dup", "dup"),
               sample_name = c("S1", "S2"), condition = c("A", "A"), replicate = c(1, 2)),
    f, sep = "\t", quote = FALSE, row.names = FALSE
  )
  expect_error(readExpDesign(f, type = "TMT"), "Duplicated 'sample'")
  expect_error(readExpDesign(f, type = "LFQ"), "Duplicated 'sample'")
})

test_that("readExpDesign stops when file column is missing for DIA", {
  f <- tempfile(fileext = ".tsv")
  write.table(
    data.frame(sample = "S1", sample_name = "S1", condition = "A", replicate = 1),
    f, sep = "\t", quote = FALSE, row.names = FALSE
  )
  expect_error(readExpDesign(f, type = "DIA"), "'file' column is not present")
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
