sumna <- function(x) {
  s <- sum(is.na(x))
  return(s)
}

#' @export
bar_number <- function(prot_data, tag = "Number of quantified proteins") {
  # prot_data = proteome_data se object
  # tag = "Number of quantified proteins"

  plexes <- unique(colData(prot_data)$plex)
  nplex <- length(plexes)
  prot_nu <- rep(0, nplex)

  for (i in 1:nplex)
  {
    cols <- colData(prot_data)$sample_name[which(colData(prot_data)$plex == plexes[i])]

    d <- assay(prot_data)[, cols]
    sumnad <- apply(d, 1, sumna)
    notAllNA <- which(sumnad < length(cols))

    prot_nu[i] <- length(notAllNA)
  }
  prot_nu_med <- median(prot_nu)
  data <- data.frame(plex = plexes, num = prot_nu)
  p <- ggplot(data, aes(x = plex, y = num)) +
    geom_bar(stat = "identity") +
    ggtitle(paste0("Median: ", prot_nu_med)) +
    theme_bw()
  return(p)
}

# plot density plot for each sample, colored by condition
#' @export
get_density <- function(se,
                        tag = "") {
  # #
  # prot_d = proteome_ratio_data_5keep_imputed_combat, se object
  # prot_annot = proteomic_annot_noRef
  # tag = "Proteome gene"
  # pdf_output_filename = "G:/My Drive/CPTAC_AML/MS_QC_20230123/filter5_knn_combat/proteome_gene_ratio_density.pdf"
  #

  prot_d <- assay(se)
  prot_annot <- colData(se)
  is.na(prot_d) <- sapply(prot_d, is.infinite)

  ### always need to deal with  -Inf, it might be a result of log2 on 0

  ym <- c()

  for (i in 1:ncol(prot_d))
  {
    d <- density(prot_d[, i], na.rm = T)

    ym <- c(ym, max(d$y))
  }




  d <- density(prot_d[, 1], na.rm = T)
  plot(density(prot_d[, 1], na.rm = T),
    ylim = c(0, max(ym)), xlim = c(min(prot_d, na.rm = T), max(prot_d, na.rm = T)), col = "white",
    main = paste0(tag, " density")
  )

  ### stratify on condition
  conditions <- unique(prot_annot$condition)
  for (i in 1:length(conditions)) {
    c <- conditions[i]
    cols <- prot_annot$sample_name[prot_annot$condition == c]
    data <- prot_d[, cols]
    for (j in 1:length(cols)) {
      this_dens <- density(data[, j], na.rm = T)
      lines(this_dens, col = palette()[i])
    }
  }
  legend(
    x = "topright", # Position
    legend = conditions, # Legend texts
    col = palette()[1:length(conditions)], # Line colors
    lwd = 2
  ) # Line width
  p <- recordPlot()
  return(p)
}
