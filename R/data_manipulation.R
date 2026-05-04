#' Remove samples from SummarizedExperiment
#'
#' Removes specified samples (columns) from a SummarizedExperiment object.
#'
#' @param se A \code{SummarizedExperiment} object.
#' @param samples Character vector of sample names (column names) to remove.
#'
#' @return A \code{SummarizedExperiment} object with the specified samples removed.
#'
#' @examples
#' \dontrun{
#' # Remove outlier samples
#' se_filtered <- remove_samples(se, c("Sample1", "Sample2"))
#' }
#'
#' @seealso \code{\link{average_samples}}
#'
#' @export
remove_samples <- function(se, samples){
  return(se[, !colnames(se) %in% samples])
}

#' Average replicate samples
#'
#' Averages intensity values across specified samples and keeps the result
#' under a single sample label, removing the other samples.
#'
#' @param se A \code{SummarizedExperiment} object.
#' @param samples Character vector of sample names to average together.
#' @param keep_label Character string specifying the sample name under which
#'   to store the averaged values. Must be one of the samples in \code{samples}.
#'
#' @return A \code{SummarizedExperiment} object with averaged values stored
#'   in the \code{keep_label} column and other specified samples removed.
#'
#' @details
#' The row means are calculated using \code{na.rm = TRUE}, so missing values
#' in some replicates will not result in NA for the averaged value. This is
#' useful for combining technical replicates or problematic sample pairs.
#'
#' @examples
#' \dontrun{
#' # Average three technical replicates
#' se_averaged <- average_samples(se,
#'                                samples = c("Sample1_R1", "Sample1_R2", "Sample1_R3"),
#'                                keep_label = "Sample1_R1")
#' }
#'
#' @seealso \code{\link{remove_samples}}
#'
#' @importFrom SummarizedExperiment assay assay<-
#'
#' @export
average_samples <- function(se, samples, keep_label){
  assay(se)[, keep_label] <- rowMeans(assay(se)[,samples], na.rm = T)
  rm_samples <- samples[!samples %in% keep_label]
  se <- remove_samples(se, rm_samples)
  return(se)
}
