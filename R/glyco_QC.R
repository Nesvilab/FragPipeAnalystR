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

