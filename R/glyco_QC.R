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

# generate a barplot for number of glycoforms based on categories
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

# generate a barplot for number of glycoforms across TMT plex set, similar to plot_feature_numbers
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
