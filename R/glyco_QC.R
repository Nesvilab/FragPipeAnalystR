N_glycan_property <- function(glycan_string){
  # from https://www.sciencedirect.com/science/article/pii/S1535947620351070, Table II

  glycan_composition <- list()
  monosaccharides <- c('HexNAc', 'Hex', 'Fuc',  'NeuAc', 'NeuGc')
  for (i in 1:length(monosaccharides)){

    temp <- as.numeric(gsub(paste0(monosaccharides[i],"(\\d+)(\\w+)?"), "\\1", glycan_string))
    if (is.na(temp)) {
      glycan_composition[monosaccharides[i]] <- 0
    } else {
      glycan_composition[monosaccharides[i]] <- temp
    }
    glycan_string <- gsub(paste0(monosaccharides[i],"\\d+"), "", glycan_string)
  }

  if (glycan_composition$Hex >= 5 & glycan_composition$HexNAc <= 2 & glycan_composition$Fuc <= 1){
    return('oligomannose')
  }
  if (glycan_composition$Fuc > 0){
    if (glycan_composition$NeuAc > 0 | glycan_composition$NeuGc > 0){
      return('fuco-sialylated')
    } else {
      return('fucosylated')
    }
  } else {
    if (glycan_composition$NeuAc > 0 | glycan_composition$NeuGc > 0)
      return('sialylated')
    else{
      return('neutral')
    }
  }
}

#' Plot glycan distribution by category
#'
#' Creates a bar plot showing the distribution of N-glycan types based on
#' their structural properties (sialylation, fucosylation, etc.).
#'
#' @param se A \code{SummarizedExperiment} object containing glycoproteomics
#'   data at the glycan level.
#' @param legacy Logical indicating whether to use legacy parsing for FragPipe
#'   versions before 23.0. Default is \code{FALSE}.
#'
#' @return A \code{ggplot} object showing a bar plot with glycan categories:
#'   \itemize{
#'     \item \code{sialylated}: Contains NeuAc or NeuGc, no fucose
#'     \item \code{fuco-sialylated}: Contains both fucose and sialic acid
#'     \item \code{fucosylated}: Contains fucose, no sialic acid
#'     \item \code{neutral}: No fucose or sialic acid, not oligomannose
#'     \item \code{oligomannose}: High mannose (Hex>=5, HexNAc<=2, Fuc<=1)
#'   }
#'
#' @details
#' The glycan classification is based on the monosaccharide composition
#' following standard N-glycan nomenclature. The color scheme follows
#' conventions used in glycomics visualization.
#'
#' @examples
#' \dontrun{
#' # Plot glycan distribution
#' plot_glycan_distribution(glyco_se)
#'
#' # For older FragPipe versions
#' plot_glycan_distribution(glyco_se, legacy = TRUE)
#' }
#'
#' @seealso \code{\link{plot_glycan_feature_numbers}}
#'
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual theme_bw theme
#'   element_blank element_line
#'
#' @export
plot_glycan_distribution <- function(se, legacy=F) {
  if (legacy) { # before FragPipe 23.0
    df <- as.data.frame(table(sapply(gsub(" _.*", "", gsub(".*_Hex", "Hex", rownames(se))), N_glycan_property)))
  } else {
    df <- as.data.frame(table(sapply( gsub(".*_", "", gsub("_\\d+\\.\\d+$", "", rownames(se))), N_glycan_property)))
  }
  colnames(df)[1] <- "Category"
  df$Category <- factor(df$Category, levels = c("sialylated",
                                                "fuco-sialylated",
                                                "fucosylated",
                                                "neutral",
                                                "oligomannose"))
  p <- ggplot(df, aes(x=Category, y=Freq, fill=Category)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("sialylated"="#BC8867",
                                 "fuco-sialylated"="#5C7099",
                                 "fucosylated"="#699870",
                                 "neutral"="#A45F61",
                                 "oligomannose"="#8077A1")) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  return(p)
}

#' Plot glycan feature numbers across TMT plexes
#'
#' Creates a stacked bar plot showing the number of glycoforms identified
#' in each TMT plex set, colored by glycan category.
#'
#' @param se A \code{SummarizedExperiment} object containing glycoproteomics
#'   data at the glycan level. Must have \code{level = "glycan"} in metadata.
#' @param exp Character string specifying the experiment type. If \code{NULL},
#'   uses the value from metadata. Currently only "TMT" is supported.
#' @param feature Character string for the feature label in the plot title.
#'   Default is \code{NULL}.
#' @param legacy Logical indicating whether to use legacy parsing for FragPipe
#'   versions before 23.0. Default is \code{TRUE}.
#'
#' @return A \code{ggplot} object showing a stacked bar plot where:
#'   \itemize{
#'     \item X-axis: TMT plex identifiers
#'     \item Y-axis: Number of glycoforms
#'     \item Fill colors: Glycan categories (sialylated, fucosylated, etc.)
#'   }
#'   Only glycoforms with complete quantification across all samples within
#'   each plex are counted.
#'
#' @examples
#' \dontrun{
#' # Plot glycan numbers across plexes
#' plot_glycan_feature_numbers(glyco_se)
#' }
#'
#' @seealso \code{\link{plot_glycan_distribution}}, \code{\link{plot_feature_numbers}}
#'
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual scale_y_continuous
#'   labs theme_bw theme element_blank element_line element_text
#' @importFrom SummarizedExperiment assay colData metadata
#' @importFrom dplyr filter if_all
#'
#' @export
plot_glycan_feature_numbers <- function(se, exp=NULL, feature=NULL, legacy=T) {
  # Show error if input is not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  assertthat::assert_that(metadata(se)$level == "glycan")
  
  if (is.null(exp)) {
    exp <- metadata(se)$exp
  }
  
  if (exp == "TMT") {
    unique_plexes <- unique(colData(se)$plex)
    n_plex <- length(unique_plexes)
    combined_df <- data.frame()
    for(i in 1:length(unique_plexes)){
      glycan_strings <- assay(se[, se$plex == unique_plexes[i]]) %>%
        data.frame() %>%
        filter(if_all(everything(), ~!is.na(.))) %>%
        rownames()
      if (legacy) {
        df <- as.data.frame(table(sapply(gsub(" _.*", "", gsub(".*_Hex", "Hex", glycan_strings)), N_glycan_property)))
      } else {
        df <- as.data.frame(table(sapply(gsub(".*_", "", gsub("_\\d+\\.\\d+$", "", glycan_strings)), N_glycan_property)))
      }
      colnames(df)[1] <- "Category"
      df$Category <- factor(df$Category, levels = c("sialylated",
                                                    "fuco-sialylated",
                                                    "fucosylated",
                                                    "neutral",
                                                    "oligomannose"))
      df$plex <- unique_plexes[i]
      combined_df <- rbind(combined_df, df)
    }
    p <- ggplot(combined_df, aes(x = plex, y = Freq, fill = Category)) +
      geom_bar(stat="identity", width = 0.85, position = "stack") +
      scale_fill_manual(values = c("sialylated"="#BC8867",
                                   "fuco-sialylated"="#5C7099",
                                   "fucosylated"="#699870",
                                   "neutral"="#A45F61",
                                   "oligomannose"="#8077A1")) +
      scale_y_continuous(expand = c(0,0)) +
      labs(title = paste0("Number of Glycan", feature , " across Plex Sets"),
           x = "Plex", y = paste0("Number of ", feature)) +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  } else {
    p <- ggplot() # currently not supported
  }
  return (p)
}
