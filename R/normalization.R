#' Global normalization (median centering with MAD scaling)
#'
#' Performs global normalization on proteomics data by median centering
#' followed by median absolute deviation (MAD) scaling.
#'
#' @param se A \code{SummarizedExperiment} object containing log2-transformed
#'   proteomics data.
#'
#' @return A \code{SummarizedExperiment} object with normalized assay values.
#'
#' @details
#' The normalization consists of two steps:
#' \enumerate{
#'   \item Median centering: Subtract the column median from each value
#'   \item MAD scaling: Divide by the column MAD, then multiply by the median MAD
#' }
#' This ensures all samples have similar center and spread.
#'
#' @examples
#' \dontrun{
#' se_normalized <- GN_normalization(se)
#' }
#'
#' @seealso \code{\link{MD_normalization}}, \code{\link{VSN_normalization}}
#'
#' @importFrom SummarizedExperiment assay assay<-
#' @importFrom matrixStats colMedians
#' @importFrom stats mad median
#'
#' @export
GN_normalization <- function(se) {
  data <- assay(se)
  MD <- sweep(data, 2, colMedians(data, na.rm = T))
  MAD <- apply(MD, 2, mad, na.rm = T)
  MAD_0 <- median(MAD)
  assay(se) <- sweep(MD, 2, MAD, FUN = "/") * MAD_0
  return(se)
}

#' Median centering normalization
#'
#' Performs median centering normalization on proteomics data by subtracting
#' the column median from each value.
#'
#' @param se A \code{SummarizedExperiment} object containing log2-transformed
#'   proteomics data.
#'
#' @return A \code{SummarizedExperiment} object with median-centered assay values.
#'
#' @details
#' This is a simple normalization that centers each sample around zero by
#' subtracting the sample median. It assumes that the majority of proteins
#' are not changing between samples.
#'
#' @examples
#' \dontrun{
#' se_normalized <- MD_normalization(se)
#' }
#'
#' @seealso \code{\link{GN_normalization}}, \code{\link{VSN_normalization}}
#'
#' @importFrom SummarizedExperiment assay assay<-
#' @importFrom matrixStats colMedians
#'
#' @export
MD_normalization <- function(se) {
  data <- assay(se)
  assay(se) <- sweep(data, 2, colMedians(data, na.rm = T))
  return(se)
}

#' Variance stabilizing normalization (VSN)
#'
#' Performs variance stabilizing normalization on proteomics data using the
#' vsn package. This method is particularly useful for LFQ and DIA data.
#'
#' @param se A \code{SummarizedExperiment} object containing log2-transformed
#'   proteomics data. Must have \code{exp} set to "LFQ" or "DIA" in metadata.
#'
#' @return A \code{SummarizedExperiment} object with VSN-normalized assay values.
#'
#' @details
#' VSN transforms the data to stabilize the variance across the intensity range.
#' The function first converts log2 intensities back to linear scale, applies
#' VSN, and stores the transformed values. This normalization is only applied
#' to LFQ and DIA data types.
#'
#' @examples
#' \dontrun
#' se_normalized <- VSN_normalization(se)
#' }
#'
#' @seealso \code{\link{GN_normalization}}, \code{\link{MD_normalization}}
#'
#' @importFrom SummarizedExperiment assay assay<-
#' @importFrom S4Vectors metadata
#' @importFrom vsn vsnMatrix predict
#' @importFrom assertthat assert_that
#'
#' @export
VSN_normalization <- function(se) {
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  data <- assay(se)
  if (metadata(se)$exp %in% c("LFQ", "DIA")) {
    vsn.fit <- vsn::vsnMatrix(2 ^ data)
    assay(se) <- vsn::predict(vsn.fit, 2 ^ data)
  }
  return(se)
}

#' PTM normalization using protein-level data
#'
#' Normalizes post-translational modification (PTM) site data by regressing
#' out the protein-level abundance changes, isolating the PTM-specific signal.
#'
#' @param ptm_se A \code{SummarizedExperiment} object containing PTM-level data
#'   (e.g., phosphosite quantification). Must have "Index" and "ProteinID"
#'   columns in rowData.
#' @param se A \code{SummarizedExperiment} object containing protein-level data.
#'   Must have "Index" and "ProteinID" columns in rowData.
#' @param print_progress Logical indicating whether to print progress messages
#'   during normalization. Default is \code{FALSE}.
#'
#' @return A \code{SummarizedExperiment} object containing the protein-normalized
#'   PTM data. The assay values represent the residuals from regressing PTM
#'   intensities on protein intensities.
#'
#' @details
#' This normalization is useful for phosphoproteomics and other PTM studies
#' where you want to identify PTM changes that are independent of protein
#' abundance changes. The method:
#' \enumerate{
#'   \item Matches PTM sites to their parent proteins
#'   \item Fits a linear model: PTM ~ Protein
#'   \item Returns the residuals as the normalized PTM values
#' }
#'
#' The correlation between normalized PTM values and protein values should be
#' near zero, indicating successful removal of protein-level effects.
#'
#' @examples
#' \dontrun{
#' # Normalize phosphosite data using protein data
#' phospho_normalized <- PTM_normalization(phospho_se, protein_se)
#' }
#'
#' @seealso \code{\link{GN_normalization}}, \code{\link{MD_normalization}}
#'
#' @importFrom SummarizedExperiment assay assay<- rowData colData
#' @importFrom dplyr filter left_join mutate select
#' @importFrom stats lm cor
#'
#' @export
PTM_normalization <- function(ptm_se, se, print_progress=F) {
  pprot <- gsub("_.*", "", rowData(ptm_se)$Index)
  inter_prot <- intersect(pprot, rowData(se)$Index)

  prot_df <- cbind(rowData(se)[c("ProteinID")], assay(se))
  psite_df <- cbind(rowData(ptm_se)[c("ProteinID", "Index")], assay(ptm_se))
  psite_df <- psite_df[psite_df$ProteinID %in% inter_prot,]

  inter_sample <- intersect(colData(ptm_se)$sample_name, colData(se)$sample_name)

  psite_df <- as.data.frame(psite_df[,c("ProteinID", "Index", inter_sample)])
  colnames(psite_df) <- c("ProteinID", "Index", inter_sample)
  prot_df <- as.data.frame(prot_df[,c("ProteinID", inter_sample)])
  colnames(prot_df) <- c("ProteinID", inter_sample)

  prot_sort_df <- data.frame(ProteinID = inter_prot, stringsAsFactors = F) %>%
    dplyr::left_join(prot_df, by = "ProteinID")

  tl <- nrow(psite_df)* length(inter_sample)
  all_psite <- rep(0, tl)
  all_prot <- rep(0, tl)
  all_index <- rep(NA_character_, tl)
  record_pos_start <- rep(0, length(inter_prot))
  record_pos_end <- rep(0, length(inter_prot))
  prot1 <- inter_prot[1]

  sel_phos1 <- psite_df %>%
    dplyr::filter(ProteinID == prot1)

  phos_d = c(t(as.matrix(sel_phos1[,inter_sample])))
  prot_d = rep(c(as.matrix(prot_sort_df[1, inter_sample])), nrow(sel_phos1))
  index_d = rep(sel_phos1$Index, each = length(inter_sample))

  record_pos_start[1] = 1
  record_pos_end[1] = length(phos_d)
  all_psite[record_pos_start[1]: record_pos_end[1]] = phos_d
  all_prot[record_pos_start[1]: record_pos_end[1]] = prot_d
  all_index[record_pos_start[1]: record_pos_end[1]] = index_d

  for(i in 2:length(inter_prot))
  {
    prot = inter_prot[i]
    sel_phos = psite_df%>%
      dplyr::filter(ProteinID == prot)
    phos_d = c(t(as.matrix(sel_phos[,inter_sample])))
    prot_d = rep(c(as.matrix(prot_sort_df[i, inter_sample])), nrow(sel_phos))
    index_d = rep(sel_phos$Index, each = length(inter_sample))
    record_pos_start[i] = record_pos_end[i-1]+1
    record_pos_end[i] = record_pos_start[i]-1+ length(phos_d)
    all_psite[record_pos_start[i]: record_pos_end[i]] = phos_d
    all_prot[record_pos_start[i]: record_pos_end[i]] = prot_d
    all_index[record_pos_start[i]: record_pos_end[i]] = index_d
    if(print_progress & i%%1000 ==0) {
      print(paste0(i, "/", length(inter_prot),"\n"))
    }
  }


  rm(prot_df)

  p_df = data.frame(phospho = all_psite, prot = all_prot,
                    caseID = rep(inter_sample, nrow(psite_df)),
                    Index = all_index,
                    stringsAsFactors = F)%>%
    na.omit()

  rm(psite_df)
  rm(all_psite)
  rm(all_prot)
  rm(all_index)

  p = lm(phospho~prot, data = p_df)

  #
  # Call:
  #   lm(formula = phospho ~ prot, data = p_df)
  #
  # Coefficients:
  #   (Intercept)         prot
  # -0.02729      0.70707
  #

  res = p$residuals
  p_df = p_df%>%
    dplyr::mutate(subPsite = res)

  rm(res)

  result_cor = cor(p_df$subPsite, p_df$prot)
  cat("correlation with prot after normalization: ",  result_cor, "\n")

  # -4.872899e-16

  sites_windows = p_df%>%
    dplyr::select(Index)%>%
    unique()



  cat("get sites_windows", "\n")

  rm(p)

  dm = matrix(NA, nrow = nrow(sites_windows), ncol = length(inter_sample))

  for(i in 1:length(inter_sample))
  {
    this_name = inter_sample[i]
    t = p_df%>%
      dplyr::filter(caseID == this_name)

    m_df = sites_windows%>%
      dplyr::left_join(t, by = c("Index"))

    dm[,i] = m_df$subPsite
  }
  cat("get dm", "\n")
  rm(p_df)

  subpsite = data.frame(sites_windows,  dm, stringsAsFactors = F)
  rm(dm)
  colnames(subpsite) = c("Index", inter_sample)
  rownames(subpsite) <- subpsite$Index
  cat("get subpsite", "\n")

  normalized_se <- ptm_se
  normalized_se <- normalized_se[subpsite$Index, inter_sample]
  assay(normalized_se) <- subpsite[,-c(1)]

  return(normalized_se)
}
