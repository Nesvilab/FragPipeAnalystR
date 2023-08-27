
# generate a barplot for number of PSM across TMT plex sets
#' @export
PSM_barplot <- function(result_dir) {
  total_count <- count()
  files <- Sys.glob(paste0(result_dir, "/*/", "psm.tsv"))
  for (i in 1:length(files)){
    temp <- fread(files[i], data.table = F)
    total_count <- c(total_count, nrow(temp))
  }

  temp <- data.frame(plex=1:length(files), total_psm=total_count)
  psm_bar <- ggplot(temp, aes(x=as.factor(plex), y=total_psm)) +
    geom_bar(stat = "identity") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.ticks.x = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank(),
          axis.text.x = element_text(vjust = 13), plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0("median: ", median(temp$total_psm))) +
    ylab("Number of PSM") +
    xlab("Plex")
  return(psm_bar)
}

# generate a barplot for number of glyco PSM across TMT plex sets
#' @export
glycoPSM_barplot <- function(result_dir, qval_filter=F, qval_threshould=0.01) {
  total_count <- c()
  total_glyco_count <- c()
  if (qval_filter) {
    for (i in 1:length(files)){
      temp <- fread(files[i], data.table = F)
      total_count <- c(total_count, nrow(temp))
      temp <- temp[!is.na(temp$`Glycan q-value`),]
      temp <- temp[temp$`Glycan q-value` <= qval_threshould,]
      total_glyco_count <- c(total_glyco_count, nrow(temp))
    }
  } else {
    for (i in 1:length(files)){
      temp <- fread(files[i], data.table = F)
      total_count <- c(total_count, nrow(temp))
      temp <- temp[!is.na(temp$`Glycan q-value`),]
      total_glyco_count <- c(total_glyco_count, nrow(temp))
    }
  }
  
  temp <- data.frame(plex=1:length(files), total_psm=total_count, glyco_psm=total_glyco_count)
  psm_bar <- ggplot(temp, aes(x=as.factor(plex), y=glyco_psm)) +
    geom_bar(stat = "identity") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.ticks.x = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank(),
          axis.text.x = element_text(vjust = 13), plot.title = element_text(hjust = 0.5)) +
    ylab("Number of Glyco PSM") +
    xlab("Plex")

  if (qval_filter) {
    psm_bar <- psm_bar + ggtitle(paste0("median: ", median(temp$glyco_psm)))
  } else {
    psm_bar <- psm_bar + ggtitle(paste0("median: ", median(temp$glyco_psm), ", q value <= ", qval_threshould))
  }
  return(psm_bar)
}
