#' Filter proteins by missing values globally
#'
#' \code{global_filter} removes proteins that have missing values in more than
#' a given percentage of all samples.
#'
#' @param se SummarizedExperiment, proteomics data.
#' @param percentage Numeric, maximum percentage of missing values allowed
#'   across all samples (0-100). Default is 50.
#' @return A filtered SummarizedExperiment object.
#' @examples
#' data("ccrcc", package = "FragPipeAnalystR")
#' filtered <- global_filter(ccrcc, percentage = 50)
#'
#' @export
global_filter <- function(se, percentage=50){
  percentage <- percentage / 100
  ridx <- rowSums(is.na(assay(se))) / ncol(assay(se)) <= percentage
  se <- se[ridx,]
  return(se)
}

#' Filter proteins by missing values per condition
#'
#' \code{filter_by_condition} retains proteins that are observed in at least
#' \code{min_percentage} of samples in at least one condition.
#'
#' @param se SummarizedExperiment, proteomics data. Must have a \code{condition}
#'   column in \code{colData}.
#' @param min_percentage Numeric, minimum percentage of non-missing values
#'   required within at least one condition (0-100). Default is 50.
#' @return A filtered SummarizedExperiment object.
#' @examples
#' data("ccrcc", package = "FragPipeAnalystR")
#' filtered <- filter_by_condition(ccrcc, min_percentage = 50)
#'
#' @export
filter_by_condition <- function(se, min_percentage=50) {
  min_percentage <- min_percentage / 100
  conditions <- unique(colData(se)$condition)
  row_ids <- rep(0, nrow(assay(se)))
  for (c in conditions){
    se_c <- se[,colData(se)$condition == c]
    ridx <- rowSums(!is.na(assay(se_c))) / ncol(assay(se_c)) >= min_percentage
    row_ids <- row_ids + ridx
  }
  se <- se[row_ids > 0,]
  return(se)
}
