test_that("run analysis using test ccrcc data", {
  data <- readResultRData(test_path("testdata", "ccrcc.rda"))
  expect_true(ggplot2::is.ggplot(plot_pca(data)))
  plot_missval_heatmap(data)
  plot_correlation_heatmap(data)
  expect_true(ggplot2::is.ggplot(plot_feature_numbers(data)))
})
