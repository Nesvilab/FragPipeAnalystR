#' Plot PSM counts across TMT plex sets
#'
#' Creates a bar plot showing the total number of peptide-spectrum matches
#' (PSMs) identified in each TMT plex set by reading psm.tsv files from
#' the FragPipe output directory.
#'
#' @param result_dir Character string specifying the path to the FragPipe
#'   results directory containing subdirectories with psm.tsv files.
#'
#' @return A \code{ggplot} object showing a bar plot where:
#'   \itemize{
#'     \item X-axis: Plex number
#'     \item Y-axis: Number of PSMs
#'     \item Title: Displays the median PSM count across plexes
#'   }
#'
#' @details
#' The function reads all psm.tsv files from subdirectories within the
#' specified result directory. This is useful for quality control to
#' identify plexes with unusually low or high PSM counts.
#'
#' @examples
#' \dontrun{
#' # Plot PSM counts from FragPipe output
#' PSM_barplot("/path/to/fragpipe/results")
#' }
#'
#' @seealso \code{\link{glycoPSM_barplot}}
#'
#' @importFrom ggplot2 ggplot aes geom_bar theme ggtitle ylab xlab
#'   element_blank element_line element_text
#' @importFrom data.table fread
#'
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

#' Plot glyco-PSM counts across TMT plex sets
#'
#' Creates a bar plot showing the number of glycopeptide-spectrum matches
#' identified in each TMT plex set, optionally filtered by q-value.
#'
#' @param result_dir Character string specifying the path to the FragPipe
#'   results directory containing subdirectories with psm.tsv files.
#' @param qval_filter Logical indicating whether to apply q-value filtering
#'   to glyco-PSMs. Default is \code{FALSE}.
#' @param qval_threshould Numeric value for the glycan q-value threshold
#'   when \code{qval_filter = TRUE}. Default is 0.01.
#'
#' @return A \code{ggplot} object showing a bar plot where:
#'   \itemize{
#'     \item X-axis: Plex number
#'     \item Y-axis: Number of glyco-PSMs
#'     \item Title: Displays median count and q-value threshold if applied
#'   }
#'
#' @details
#' Glyco-PSMs are identified as PSMs with non-NA values in the "Glycan q-value"
#' column. This function is useful for quality control of glycoproteomics
#' experiments to ensure consistent glycopeptide identification across plexes.
#'
#' @examples
#' \dontrun{
#' # Plot all glyco-PSMs
#' glycoPSM_barplot("/path/to/fragpipe/results")
#'
#' # Plot filtered glyco-PSMs (q-value <= 0.01)
#' glycoPSM_barplot("/path/to/fragpipe/results", qval_filter = TRUE,
#'                  qval_threshould = 0.01)
#' }
#'
#' @seealso \code{\link{PSM_barplot}}, \code{\link{plot_glycan_distribution}}
#'
#' @importFrom ggplot2 ggplot aes geom_bar theme ggtitle ylab xlab
#'   element_blank element_line element_text
#' @importFrom data.table fread
#'
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
