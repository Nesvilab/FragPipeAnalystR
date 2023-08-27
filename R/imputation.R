#' Impute missing values
#' customized from https://github.com/arnesmits/DEP/blob/b425d8d0db67b15df4b8bcf87729ef0bf5800256/R/functions.R
#'
#' \code{impute} imputes missing values in a proteomics dataset.
#'
#' @param se SummarizedExperiment,
#' Proteomics data (output from \code{\link{make_se}()} or
#' \code{\link{make_se_parse}()}). It is adviced to first remove
#' proteins with too many missing values using \code{\link{filter_missval}()}
#' and normalize the data using \code{\link{normalize_vsn}()}.
#' @param fun "bpca", "knn", "QRILC", "MLE", "MinDet",
#' "MinProb", "man", "min", "zero", "mixed" or "nbavg",
#' Function used for data imputation based on \code{\link{manual_impute}}
#' and \code{\link[MSnbase:impute-methods]{impute}}.
#' @param ... Additional arguments for imputation functions as depicted in
#' \code{\link{manual_impute}} and \code{\link[MSnbase:impute-methods]{impute}}.
#' @return An imputed SummarizedExperiment object.
#' @examples
#' # Load example
#' data("ccrcc", package = "FragPipeR")
#' ccrcc_impute <- impute(ccrcc, "knn")
#'
#' @export
impute <- function(se, fun = c(
                     "bpca", "knn", "QRILC", "MLE",
                     "MinDet", "MinProb", "man", "min", "zero",
                     "mixed", "nbavg"
                   ), ...) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(
    inherits(se, "SummarizedExperiment"),
    is.character(fun)
  )

  # Show error if inputs do not contain required columns
  fun <- match.arg(fun)

  if (any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '",
      deparse(substitute(se)),
      "'\nRun make_unique() and make_se() to obtain the required columns",
      call. = FALSE
    )
  }

  # Annotate whether or not there are missing values and how many
  rowData(se)$imputed <- apply(is.na(assay(se)), 1, any)
  rowData(se)$num_NAs <- rowSums(is.na(assay(se)))
  se <- se[!rowData(se)$num_NAs == dim(se)[2],]

  # Show error if there are no missing values
  if (!any(is.na(assay(se)))) {
    warning("No missing values in '", deparse(substitute(se)), "'. ",
      "Returning the unchanged object.",
      call. = FALSE
    )
    return(se)
  }

  # if the "man" function is selected, use the manual imputation method
  if (fun == "man") {
    se <- manual_impute(se, ...)
  }
  # else use the MSnSet::impute function
  else {
    MSnSet_data <- as(se, "MSnSet")
    MSnSet_imputed <- MSnbase::impute(MSnSet_data, method = fun, ...)
    assay(se) <- MSnbase::exprs(MSnSet_imputed)
  }

  return(se)
}

#' Imputation by random draws from a manually defined distribution
#' customized from https://github.com/arnesmits/DEP/blob/b425d8d0db67b15df4b8bcf87729ef0bf5800256/R/functions.R
#'
#' \code{manual_impute} imputes missing values in a proteomics dataset
#' by random draws from a manually defined distribution.
#'
#' @param se SummarizedExperiment,
#' Proteomics data (output from \code{\link{make_se}()} or
#' \code{\link{make_se_parse}()}). It is adviced to first remove
#' proteins with too many missing values using \code{\link{filter_missval}()}
#' and normalize the data using \code{\link{normalize_vsn}()}.
#' @param shift Numeric(1),
#' Sets the left-shift of the distribution (in standard deviations) from
#' the median of the original distribution.
#' @param scale Numeric(1),
#' Sets the width of the distribution relative to the
#' standard deviation of the original distribution.
#' @return An imputed SummarizedExperiment object.
#' @export
manual_impute <- function(se, scale = 0.3, shift = 1.8) {
  if (is.integer(scale)) scale <- is.numeric(scale)
  if (is.integer(shift)) shift <- is.numeric(shift)
  # Show error if inputs are not the required classes
  assertthat::assert_that(
    inherits(se, "SummarizedExperiment"),
    is.numeric(scale),
    length(scale) == 1,
    is.numeric(shift),
    length(shift) == 1
  )

  se_assay <- assay(se)

  # Show error if there are no missing values
  if (!any(is.na(se_assay))) {
    stop("No missing values in '", deparse(substitute(se)), "'",
      call. = FALSE
    )
  }

  # Get descriptive parameters of the current sample distributions
  stat <- se_assay %>%
    data.frame(check.names = F) %>%
    rownames_to_column() %>%
    gather(samples, value, -rowname) %>%
    filter(!is.na(value)) %>%
    group_by(samples) %>%
    summarise(
      mean = mean(value),
      median = median(value),
      sd = sd(value),
      n = n(),
      infin = nrow(se_assay) - n
    )
  # Impute missing values by random draws from a distribution
  # which is left-shifted by parameter 'shift' * sd and scaled by parameter 'scale' * sd.
  for (a in seq_len(nrow(stat))) {
    assay(se)[is.na(assay(se)[, stat$samples[a]]), stat$samples[a]] <-
      rnorm(stat$infin[a],
        mean = stat$median[a] - shift * stat$sd[a],
        sd = stat$sd[a] * scale
      )
  }
  return(se)
}
