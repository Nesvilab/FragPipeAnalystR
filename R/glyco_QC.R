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

#' Get glycoforms per glycosite
#'
#' @param se A \code{SummarizedExperiment} object at the glycan level.
#' @return A data frame with columns \code{site} and \code{n_glycoforms}.
#' @export
get_glycoforms_per_site <- function(se) {
  ids_no_mass <- gsub(" _.*$", "", rownames(se))
  site_ids <- sapply(strsplit(ids_no_mass, "_"), function(x) paste(x[1], x[length(x) - 1], sep = "_"))
  gpsite <- as.data.frame(table(site_ids))
  colnames(gpsite) <- c("site", "n_glycoforms")
  gpsite$n_glycoforms <- as.integer(gpsite$n_glycoforms)
  gpsite
}

#' Plot glycoforms per site distribution
#'
#' @param se A \code{SummarizedExperiment} object at the glycan level.
#' @param title Plot title. Default is \code{""}.
#' @return A \code{ggplot} object.
#' @importFrom ggplot2 ggplot aes geom_col geom_text scale_x_continuous
#'   scale_y_continuous theme_classic labs expansion
#' @export
plot_glycoforms_per_site <- function(se, title = "") {
  gpsite <- get_glycoforms_per_site(se)
  count_df <- as.data.frame(table(gpsite$n_glycoforms))
  colnames(count_df) <- c("n_glycoforms", "count")
  count_df$n_glycoforms <- as.integer(as.character(count_df$n_glycoforms))
  x_max <- max(count_df$n_glycoforms)
  ggplot(count_df, aes(x = n_glycoforms, y = count)) +
    geom_col(fill = "steelblue", width = 0.7) +
    geom_text(aes(label = count), vjust = -0.3, size = 3) +
    scale_x_continuous(breaks = 1:x_max, limits = c(0.5, x_max + 0.5)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_classic() +
    labs(title = title, x = "Glycoforms per site", y = "# sites")
}

#' Get glycoforms per protein
#'
#' @param se A \code{SummarizedExperiment} object at the glycan level.
#' @return A data frame with columns \code{protein} and \code{n_glycoforms}.
#' @export
get_glycoforms_per_protein <- function(se) {
  protein_ids <- sapply(strsplit(rownames(se), "_"), `[[`, 1)
  gprot <- as.data.frame(table(protein_ids))
  colnames(gprot) <- c("protein", "n_glycoforms")
  gprot$n_glycoforms <- as.integer(gprot$n_glycoforms)
  gprot
}

#' Plot glycoforms per protein distribution
#'
#' @param se A \code{SummarizedExperiment} object at the glycan level.
#' @param title Plot title. Default is \code{""}.
#' @return A \code{ggplot} object.
#' @importFrom ggplot2 ggplot aes geom_col geom_text scale_x_continuous
#'   scale_y_continuous theme_classic labs expansion
#' @export
plot_glycoforms_per_protein <- function(se, title = "") {
  gprot <- get_glycoforms_per_protein(se)
  count_df <- as.data.frame(table(gprot$n_glycoforms))
  colnames(count_df) <- c("n_glycoforms", "count")
  count_df$n_glycoforms <- as.integer(as.character(count_df$n_glycoforms))
  x_max <- max(count_df$n_glycoforms)
  ggplot(count_df, aes(x = n_glycoforms, y = count)) +
    geom_col(fill = "steelblue", width = 0.7) +
    geom_text(aes(label = count), vjust = -0.3, size = 3) +
    scale_x_continuous(breaks = 1:x_max, limits = c(0.5, x_max + 0.5)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_classic() +
    labs(title = title, x = "Glycoforms per protein", y = "# proteins")
}

#' Get glycosites per protein
#'
#' @param se A \code{SummarizedExperiment} object at the glycosite level.
#' @return A data frame with columns \code{protein} and \code{n_glycosites}.
#' @export
get_glycosite_per_protein <- function(se) {
  protein_ids <- sapply(strsplit(rownames(se), "_"), `[[`, 1)
  gprot <- as.data.frame(table(protein_ids))
  colnames(gprot) <- c("protein", "n_glycosites")
  gprot$n_glycosites <- as.integer(gprot$n_glycosites)
  gprot
}

#' Plot glycosites per protein distribution
#'
#' @param se A \code{SummarizedExperiment} object at the glycosite level.
#' @param title Plot title. Default is \code{""}.
#' @return A \code{ggplot} object.
#' @importFrom ggplot2 ggplot aes geom_col geom_text scale_x_continuous
#'   scale_y_continuous theme_classic labs expansion
#' @export
plot_glycosite_per_protein <- function(se, title = "") {
  gprot <- get_glycosite_per_protein(se)
  count_df <- as.data.frame(table(gprot$n_glycosites))
  colnames(count_df) <- c("n_glycosites", "count")
  count_df$n_glycosites <- as.integer(as.character(count_df$n_glycosites))
  x_max <- max(count_df$n_glycosites)
  ggplot(count_df, aes(x = n_glycosites, y = count)) +
    geom_col(fill = "steelblue", width = 0.7) +
    geom_text(aes(label = count), vjust = -0.3, size = 3) +
    scale_x_continuous(breaks = 1:x_max, limits = c(0.5, x_max + 0.5)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_classic() +
    labs(title = title, x = "Glycosites per protein", y = "# proteins")
}

#' Plot glycan distribution across multiple datasets
#'
#' @param se_list A named list of \code{SummarizedExperiment} objects, all at
#'   the same metadata level.
#' @param legacy Logical; use legacy row name parsing for FragPipe < 23.0.
#'   Default is \code{FALSE}.
#' @return A \code{ggplot} object with glycan categories faceted by dataset.
#' @importFrom ggplot2 ggplot aes geom_bar theme_bw theme element_blank
#'   element_line
#' @importFrom SummarizedExperiment metadata
#' @export
plot_glycan_distribution_combined <- function(se_list, legacy = FALSE) {
  lvls <- sapply(se_list, function(se) metadata(se)$level)
  if (length(unique(lvls)) > 1) {
    stop("All SE objects must have the same metadata level. Got: ",
         paste(names(lvls), lvls, sep = "=", collapse = ", "))
  }

  df <- do.call(rbind, lapply(names(se_list), function(nm) {
    se <- se_list[[nm]]
    compositions <- if (legacy) {
      gsub(" _.*", "", gsub(".*_Hex", "Hex", rownames(se)))
    } else {
      gsub(".*_", "", gsub(" _.*$", "", rownames(se)))
    }
    d <- as.data.frame(table(sapply(compositions, N_glycan_property)))
    colnames(d) <- c("Category", "Freq")
    d$Dataset <- nm
    d
  }))

  df$Category <- factor(df$Category,
                        levels = c("sialylated", "fuco-sialylated",
                                   "fucosylated", "neutral", "oligomannose"))
  df$Dataset <- factor(df$Dataset, levels = names(se_list))

  ggplot(df, aes(x = Category, y = Freq, fill = Dataset)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw() +
    theme(panel.border     = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line        = element_line(colour = "black"))
}
