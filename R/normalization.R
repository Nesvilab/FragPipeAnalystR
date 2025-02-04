#' @export
GN_normalization <- function(se) {
  data <- assay(se)
  MD <- sweep(data, 2, colMedians(data, na.rm = T))
  MAD <- apply(MD, 2, mad, na.rm = T)
  MAD_0 <- median(MAD)
  assay(se) <- sweep(MD, 2, MAD, FUN = "/") * MAD_0
  return(se)
}

#' @export
MD_normalization <- function(se) {
  data <- assay(se)
  assay(se) <- sweep(data, 2, colMedians(data, na.rm = T))
  return(se)
}

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
  record_pos_start <- rep(0, length(inter_prot))
  record_pos_end <- rep(0, length(inter_prot))
  prot1 <- inter_prot[1]

  sel_phos1 <- psite_df %>%
    dplyr::filter(ProteinID == prot1)

  phos_d = c(t(as.matrix(sel_phos1[,inter_sample])))
  prot_d = rep(c(as.matrix(prot_sort_df[1, inter_sample])), nrow(sel_phos1))

  record_pos_start[1] = 1
  record_pos_end[1] = length(phos_d)
  all_psite[record_pos_start[1]: record_pos_end[1]] = phos_d
  all_prot[record_pos_start[1]: record_pos_end[1]] = prot_d

  for(i in 2:length(inter_prot))
  {
    prot = inter_prot[i]
    sel_phos = psite_df%>%
      dplyr::filter(ProteinID == prot)
    phos_d = c(t(as.matrix(sel_phos[,inter_sample])))
    prot_d = rep(c(as.matrix(prot_sort_df[i, inter_sample])), nrow(sel_phos))
    record_pos_start[i] = record_pos_end[i-1]+1
    record_pos_end[i] = record_pos_start[i]-1+ length(phos_d)
    all_psite[record_pos_start[i]: record_pos_end[i]] = phos_d
    all_prot[record_pos_start[i]: record_pos_end[i]] = prot_d
    if(print_progress & i%%1000 ==0) {
      print(paste0(i, "/", length(inter_prot),"\n"))
    }
  }


  rm(prot_df)

  p_df = data.frame(phospho = all_psite, prot = all_prot,
                    caseID = rep(inter_sample, nrow(psite_df)),
                    Index = rep(psite_df$Index, each = length(inter_sample)),
                    stringsAsFactors = F)%>%
    na.omit()

  rm(psite_df)
  rm(all_psite)
  rm(all_prot)

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
