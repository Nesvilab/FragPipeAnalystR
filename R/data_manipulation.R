# function to remove sample by column name
#' @export
remove_samples <- function(se, samples){
  return(se[, !colnames(se) %in% samples])
}

#' @export
average_samples <- function(se, samples, keep_label){
  assay(se)[, keep_label] <- rowMeans(assay(se)[,samples], na.rm = T)
  rm_samples <- samples[!samples %in% keep_label]
  se <- remove_samples(se, rm_samples)
  return(se)
}
