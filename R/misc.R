#' Plot coefficient of variation histogram
#'
#' Creates a histogram visualization of the coefficient of variation (CV)
#' for each protein across replicates within each condition.
#'
#' @param se A \code{SummarizedExperiment} object containing log2-transformed
#'   proteomics data.
#' @param id Character string specifying the column in colData to use for
#'   sample identification. Default is "sample_name".
#' @param scale Logical indicating whether to fix the x-axis scale from 0 to 1
#'   (\code{TRUE}) or use automatic scaling (\code{FALSE}). Default is \code{TRUE}.
#' @param check.names Logical passed to \code{data.frame()}. Default is \code{FALSE}.
#'
#' @return A \code{ggplot} object showing faceted histograms of CV values for
#'   each condition. The median CV is indicated by a dashed vertical line and
#'   displayed as text.
#'
#' @details
#' The CV is calculated as the standard deviation divided by the mean for each
#' protein within each condition, using the back-transformed (linear scale)
#' intensity values. Lower CVs indicate better reproducibility across replicates.
#'
#' @examples
#' \dontrun{
#' # Basic CV plot
#' plot_cvs(se)
#'
#' # Without fixed scaling
#' plot_cvs(se, scale = FALSE)
#' }
#'
#' @seealso \code{\link{plot_feature_numbers}}
#'
#' @importFrom ggplot2 ggplot aes geom_histogram facet_wrap geom_vline
#'   scale_x_continuous labs theme_bw theme element_text geom_text
#' @importFrom SummarizedExperiment assay colData
#' @importFrom dplyr left_join group_by summarise mutate
#' @importFrom tidyr gather
#' @importFrom tibble rownames_to_column
#' @importFrom scales percent
#'
#' @export
plot_cvs <- function(se, id="sample_name", scale=T, check.names=F) {

  ## backtransform data
  untransformed_intensity<- 2^(assay(se))
  exp_design<-colData(se)

  ### merge untransformed to exp design and calculate cvs
  if (id == "ID") {
    cvs_group<- untransformed_intensity %>% data.frame() %>%
      tibble::rownames_to_column() %>%
      tidyr::gather("ID", "Intensity", -rowname) %>%
      dplyr::left_join(.,data.frame(exp_design), by="ID") %>%
      dplyr::group_by(rowname,condition) %>%
      dplyr::summarise(cvs=coef_variation(Intensity)) %>%
      dplyr::group_by(condition)%>%
      dplyr::mutate(condition_median=median(cvs))
  } else {
    cvs_group<- untransformed_intensity %>% data.frame(check.names=check.names) %>%
      tibble::rownames_to_column() %>%
      tidyr::gather("ID", "Intensity", -rowname) %>%
      dplyr::left_join(.,data.frame(exp_design), by=c("ID"=id)) %>%
      dplyr::group_by(rowname,condition) %>%
      dplyr::summarise(cvs=coef_variation(Intensity)) %>%
      dplyr::group_by(condition)%>%
      dplyr::mutate(condition_median=median(cvs, na.rm = T))
  }
  if (scale) {
    p1 <- ggplot(cvs_group, aes(cvs, color=condition, fill=condition)) +
      geom_histogram(alpha=.5, bins= 20, show.legend = FALSE) +
      facet_wrap(~condition) +
      geom_vline(aes(xintercept=condition_median, group=condition),
                 color='grey40',
                 linetype="dashed") +
      scale_x_continuous(labels = scales::percent, limits=c(0, 1)) +
      labs(title= 'Sample Coefficient of Variation', x="Coefficient of Variation", y="Count") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5,face = "bold"), axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
  } else {
    p1 <- ggplot(cvs_group, aes(cvs, color=condition, fill=condition)) +
      geom_histogram(alpha=.5, bins= 20, show.legend = FALSE) +
      facet_wrap(~condition) +
      geom_vline(aes(xintercept=condition_median, group=condition),
                 color='grey40',
                 linetype="dashed") +
      scale_x_continuous(labels = scales::percent) +
      labs(title= 'Sample Coefficient of Variation', x="Coefficient of Variation", y="Count") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5,face = "bold"), axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
  }


  p <- p1 + geom_text(aes(x=0.9,
                          y=max(ggplot_build(p1)$data[[1]]$ymax*1.1),
                          label=paste0("Median =",round(condition_median,2)*100,"%",by="")),
                      show.legend = FALSE, size=4)
  return(p)
}

coef_variation <- function(x) {
  coef <- sd(x) / mean(x)
}

#' Plot individual protein/feature abundance
#'
#' Creates boxplot or violin plot visualizations of abundance values for
#' specific proteins or features across conditions.
#'
#' @param se A \code{SummarizedExperiment} object containing proteomics data.
#' @param protein Character vector of protein/feature identifiers to plot.
#'   Must match rownames or values in the \code{index} column.
#' @param index Character string specifying the column in rowData to use for
#'   matching protein identifiers. Default is \code{NULL} (use rownames).
#' @param type Character string specifying the plot type. Options are
#'   \code{"boxplot"} (default) or \code{"violin"}.
#' @param id Character string specifying the column in colData to use for
#'   sample identification. Default is "sample_name".
#' @param interactive Logical indicating whether to create an interactive
#'   plotly visualization. Default is \code{FALSE}.
#'
#' @return A \code{ggplot} object (when \code{interactive = FALSE}) or
#'   a \code{plotly} object (when \code{interactive = TRUE}) showing the
#'   abundance distribution for the specified proteins across conditions.
#'
#' @details
#' When multiple proteins are specified, they are displayed in separate facets.
#' Points are colored by replicate to help identify sample-specific effects.
#' Interactive plots provide hover information showing sample names.
#'
#' @examples
#' \dontrun{
#' # Boxplot for a single protein
#' plot_feature(se, protein = "EGFR")
#'
#' # Violin plot for multiple proteins
#' plot_feature(se, protein = c("EGFR", "AKT1", "MTOR"), type = "violin")
#'
#' # Interactive boxplot
#' plot_feature(se, protein = "EGFR", interactive = TRUE)
#' }
#'
#' @seealso \code{\link{plot_pca}}, \code{\link{plot_feature_numbers}}
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_violin geom_jitter
#'   position_dodge facet_wrap labs theme_bw theme element_blank element_line
#'   scale_color_brewer
#' @importFrom plotly plot_ly layout
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom dplyr left_join group_by summarize mutate
#' @importFrom tidyr gather
#' @importFrom tibble rownames_to_column
#' @importFrom readr parse_factor
#'
#' @export
plot_feature <- function(se, protein, index=NULL,
                         type="boxplot", id="sample_name", interactive=F) {
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.character(protein),
                          is.character(type))
  if (is.null(index)) {
    subset <- se[protein,]
    df_reps <- data.frame(assay(subset), check.names = F) %>%
      rownames_to_column() %>%
      gather(ID, val, -rowname) %>%
      left_join(., data.frame(colData(subset)), by = c("ID"=id))

    df_reps$rowname <- parse_factor(as.character(df_reps$rowname), levels = protein)
  } else {
    subset <- se[rowData(se)[,index] %in% protein,]
    df_reps <- data.frame(assay(subset), check.names = F) %>%
      rownames_to_column() %>%
      gather(ID, val, -rowname) %>%
      left_join(., data.frame(colData(subset)), by = c("ID"=id)) %>%
      left_join(., data.frame(rowData(subset)[,c("Index", index)]), by = c("rowname"="Index"))
    df_reps$rowname <- parse_factor(as.character(df_reps[,index]), levels = protein)
  }



  df_CI <- df_reps %>%
    group_by(condition, rowname) %>%
    summarize(mean = mean(val, na.rm = TRUE),
              sd = sd(val, na.rm = TRUE),
              n = n()) %>%
    mutate(error = qnorm(0.975) * sd / sqrt(n),
           CI.L = mean - error,
           CI.R = mean + error) %>%
    as.data.frame()
  df_CI$rowname <- parse_factor(as.character(df_CI$rowname), levels = protein)

  df_reps$condition <- as.factor(df_reps$condition)
  df_reps <- df_reps[!is.na(df_reps$val),]
  df_reps$replicate <- as.character(df_reps$replicate)
  if (interactive) {
    if(type=="violin"){
      if (max(df_reps$replicate) == 1){
        p <- plot_ly(df_reps,
                     x = ~rowname,
                     y = ~val,
                     color = ~condition,
                     text = ~sample_name,
                     hoverinfo = "text",
                     type = "violin") %>%
          plotly::layout(violinmode = "group",
                         xaxis = list(title = ''),
                         yaxis = list(title = 'Abundance'), showlegend=T)
      } else {
        p <- plot_ly(df_reps,
                     x = ~rowname,
                     y = ~val,
                     color = ~condition,
                     text = ~sample_name,
                     hoverinfo = "text",
                     type = "violin") %>%
          plotly::layout(violinmode = "group",
                         xaxis = list(title = ''),
                         yaxis = list(title = 'Abundance'), showlegend=T)
      }
      return(p)
    } else if(type=="boxplot"){
      if (max(df_reps$replicate) == 1){
        p <- plot_ly(df_reps,
                     x = ~rowname,
                     y = ~val,
                     color = ~condition,
                     text = ~sample_name,
                     hoverinfo = "text",
                     type = "box",
                     boxpoints = "all", jitter = 0.3, pointpos = 0) %>%
          plotly::layout(boxmode = "group",
                         xaxis = list(title = ''),
                         yaxis = list(title = 'Abundance'), showlegend=T)
      } else {
        p <- plot_ly(df_reps,
                     x = ~rowname,
                     y = ~val,
                     color = ~condition,
                     text = ~sample_name,
                     hoverinfo = "text",
                     type = "box",
                     boxpoints = "all", jitter = 0.3, pointpos = 0) %>%
          plotly::layout(boxmode = "group",
                         xaxis = list(title = ''),
                         yaxis = list(title = 'Abundance'), showlegend=T)
      }
    }
  } else {
    if(type=="violin"){
      if (max(df_reps$replicate) == 1){
        p <- ggplot(df_reps, aes(condition, val)) +
          geom_violin(fill="grey90", scale = "width",
                      draw_quantiles = 0.5,
                      trim =TRUE) +
          geom_jitter(size = 3, position = position_dodge(width=0.3)) +
          labs(
            y = expression(log[2]~"Intensity")) +
          facet_wrap(~rowname) +
          theme_bw() +
          theme(axis.title.x = element_blank(),
                panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
      } else {
        p<-ggplot(df_reps, aes(condition, val))+
          geom_violin(fill="grey90", scale = "width",
                      draw_quantiles = 0.5,
                      trim =TRUE) +
          geom_jitter(aes(color = factor(replicate)),
                      size = 3, position = position_dodge(width=0.3)) +
          labs(
            y = expression(log[2]~"Intensity"),
            col = "Replicates") +
          facet_wrap(~rowname) +
          scale_color_brewer(palette = "Dark2") +
          theme_bw() +
          theme(axis.title.x = element_blank(),
                panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
      }
    } else if(type=="boxplot") {
      if (max(df_reps$replicate) == 1){
        p <- ggplot(df_reps, aes(condition, val))+
          geom_boxplot()+
          geom_jitter(size = 3, position = position_dodge(width=0.3)) +
          labs(y = expression(log[2]~"Intensity")) +
          facet_wrap(~rowname) +
          theme_bw() +
          theme(axis.title.x = element_blank(),
                panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
      } else {
        p <- ggplot(df_reps, aes(condition, val))+
          geom_boxplot()+
          geom_jitter(aes(color = factor(replicate)),
                      size = 3, position = position_dodge(width=0.3)) +
          labs(
            y = expression(log[2]~"Intensity"),
            col = "Replicates") +
          facet_wrap(~rowname) +
          scale_color_brewer(palette = "Dark2") +
          theme_bw() +
          theme(axis.title.x = element_blank(),
                panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
      }
    }
  }
  return(p)
}

#' Plot number of identified features per sample
#'
#' Creates a bar plot showing the number of identified proteins/peptides/features
#' in each sample or across TMT plexes.
#'
#' @param se A \code{SummarizedExperiment} object containing proteomics data.
#' @param exp Character string specifying the experiment type. If \code{NULL},
#'   uses the value from metadata. Options are "TMT", "LFQ", "DIA".
#' @param feature Character string specifying the feature label for the plot
#'   (e.g., "Proteins", "Peptides", "Sites"). If \code{NULL}, automatically
#'   determined from metadata level.
#' @param fill Character string specifying the colData column to use for
#'   bar coloring. Default is "condition".
#'
#' @return A \code{ggplot} object showing:
#'   \itemize{
#'     \item For TMT: Stacked bar plot showing plex-wise vs project-wise
#'       quantified features
#'     \item For LFQ/DIA: Bar plot showing feature counts per sample colored
#'       by the specified fill variable
#'   }
#'
#' @details
#' For TMT experiments, the plot distinguishes between features quantified
#' only within a single plex ("Plex-wise") vs features quantified across
#' all plexes ("Project-wise").
#'
#' For LFQ/DIA experiments, each bar represents a sample and shows the
#' number of non-missing features.
#'
#' @examples
#' \dontrun{
#' # Basic feature number plot
#' plot_feature_numbers(se)
#'
#' # Color by replicate instead of condition
#' plot_feature_numbers(se, fill = "replicate")
#'
#' # Specify custom feature label
#' plot_feature_numbers(se, feature = "Phosphosites")
#' }
#'
#' @seealso \code{\link{plotCumulativeMissingPercent}}, \code{\link{plot_cvs}}
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_col scale_fill_manual
#'   scale_y_continuous labs theme_bw theme element_blank element_line
#'   element_text
#' @importFrom SummarizedExperiment assay colData metadata
#' @importFrom dplyr group_by summarize left_join filter if_all mutate
#' @importFrom tidyr gather
#' @importFrom tibble rownames_to_column
#'
#' @export
plot_feature_numbers <- function(se, exp=NULL, feature=NULL, fill="condition") {
  # Show error if input is not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))

  if (is.null(feature)) {
    if (metadata(se)$level == "protein") {
      feature <- "Proteins"
    } else if (metadata(se)$level == "gene") {
      feature <- "Peptides"
    } else if (metadata(se)$level == "peptide") {
      feature <- "Peptides"
    } else if (metadata(se)$level == "site") {
      feature <- "Sites"
    } else if (metadata(se)$level == "glycan") {
      feature <- "Glycans"
    }
  }


  if (is.null(exp)) {
    exp <- metadata(se)$exp
  }

  if (exp == "TMT") {
    unique_plexes <- unique(colData(se)$plex)
    n_plex <- length(unique_plexes)
    prot_v <- c()
    for(i in 1:length(unique_plexes)){
      n_prot <- assay(se[, se$plex == unique_plexes[i]]) %>%
        data.frame() %>%
        filter(if_all(everything(), ~!is.na(.))) %>%
        nrow()
      prot_v <- c(prot_v, n_prot)
    }
    prot_c <- nrow(na.omit(assay(se)))
    df_prot <- data.frame(plex=factor(rep(unique_plexes,2)),
                          num_protein=c(prot_v-prot_c,rep(prot_c,n_plex)),
                          Status=c(rep("Plex-wise",n_plex),rep("Project-wise",n_plex)))
    p <- ggplot(df_prot, aes(x = plex, y = num_protein, fill = Status)) +
      geom_bar(stat="identity", width = 0.85, position = "stack") +
      scale_fill_manual(values = c("Plex-wise"="#75E6DA","Project-wise"="#189AB4")) +
      scale_y_continuous(expand = c(0,0)) +
      labs(title = paste0("Number of ", feature , " across Plex Sets"),
           x = "Plex", y = paste0("Number of ", feature)) +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  } else {
    # Make a binary long data.frame (1 = valid value, 0 = missing value)
    df <- assay(se) %>%
      data.frame(check.names = F) %>%
      rownames_to_column() %>%
      gather(ID, bin, -rowname) %>%
      mutate(bin = ifelse(is.na(bin), 0, 1))

    # Summarize the number of proteins identified
    # per sample and generate a barplot
    stat <- df %>%
      group_by(ID) %>%
      summarize(n = n(), sum = sum(bin)) %>%
      left_join(., data.frame(colData(se)), by = c("ID"="sample_name"))

    p <- ggplot(stat, aes(x = ID, y = sum, fill = .data[[fill]])) +
      geom_col() +
      labs(title = paste0("Number of ", feature, " per Sample (Total Number: ", dim(assay(se))[1], ")"),
           x = "", y = paste0("Number of ", feature)) +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  return (p)
}

#' Plot cumulative feature count vs missing percentage
#'
#' Creates a line plot showing the cumulative number of features as a function
#' of the maximum allowed missing percentage. This helps visualize data
#' completeness and inform filtering decisions.
#'
#' @param se A \code{SummarizedExperiment} object containing proteomics data.
#' @param smp_lst Character vector of sample names to include in the analysis.
#'   If \code{NULL} (default), all samples are included.
#' @param title Character string for the plot title. Default is empty.
#'
#' @return A \code{ggplot} object showing a cumulative line plot where:
#'   \itemize{
#'     \item X-axis: Missing percentage threshold
#'     \item Y-axis: Cumulative number of features with missing percentage
#'       at or below the threshold
#'   }
#'   The plot includes annotations showing the number of complete features,
#'   features with less than 50% missing, and total identified features.
#'
#' @details
#' This plot is useful for deciding on missing value thresholds for filtering.
#' Features are sorted by their missing percentage, and the cumulative count
#' shows how many features would be retained at each threshold.
#'
#' @examples
#' \dontrun{
#' # Basic cumulative missing plot
#' plotCumulativeMissingPercent(se)
#'
#' # With custom title
#' plotCumulativeMissingPercent(se, title = "Data Completeness")
#'
#' # For specific samples only
#' plotCumulativeMissingPercent(se, smp_lst = c("Sample1", "Sample2", "Sample3"))
#' }
#'
#' @seealso \code{\link{plot_feature_numbers}}, \code{\link{plot_missval_heatmap}}
#'
#' @importFrom ggplot2 ggplot aes geom_line annotate labs theme_classic theme
#'   element_text
#' @importFrom SummarizedExperiment assay
#' @importFrom dplyr arrange mutate
#' @importFrom data.table setDT
#'
#' @export
#' Plot Venn diagram of features across SE objects
#'
#' @param se_list A named list of \code{SummarizedExperiment} objects, all at
#'   the same metadata level.
#' @param title Plot title. Default is \code{""}.
#' @param fill_color Character vector of fill colors. Default is
#'   \code{c("#4E79A7", "#F28E2B")}.
#' @param stroke_size Numeric stroke size. Default is \code{0.5}.
#' @param text_size Numeric text size. Default is \code{4}.
#' @return A \code{ggplot} object.
#' @importFrom ggvenn ggvenn
#' @importFrom ggplot2 ggtitle
#' @importFrom SummarizedExperiment metadata
#' @export
plot_venn_se <- function(se_list, title = "",
                         fill_color = c("#4E79A7", "#F28E2B"),
                         stroke_size = 0.5, text_size = 4) {
  lvls <- sapply(se_list, function(se) metadata(se)$level)
  if (length(unique(lvls)) > 1) {
    stop("All SE objects must have the same metadata level. Got: ",
         paste(names(lvls), lvls, sep = "=", collapse = ", "))
  }
  ggvenn(lapply(se_list, rownames),
         fill_color = fill_color,
         stroke_size = stroke_size,
         text_size   = text_size) +
    ggtitle(title)
}

plotCumulativeMissingPercent = function(se,
                                        smp_lst = NULL,
                                        title = ""){

  df <- data.frame(assay(se))
  if (is.null(smp_lst)) {
    smp_lst <- colnames(df)
  } else {
    smp_lst <- intersect(smp_lst, colnames(df))
  }


  # NA counting for each feature
  df <- data.table::setDT(df)
  na_df <- df[, .(Index = rownames(df), NA_count = rowSums(is.na(.SD))), .SDcols = smp_lst]

  # prepare data to plot
  na_df_srt <- na_df %>%
    dplyr::arrange(NA_count) %>%
    dplyr::mutate(missPercent = NA_count / length(smp_lst),
                  FeatureCount = seq_len(nrow(.)))

  # find row index with missPercent == 0.5
  row_vline <- sum(na_df_srt$missPercent <= 0.5)
  xinter <- na_df_srt$missPercent[row_vline]
  yinter <- na_df_srt$FeatureCount[row_vline]
  test_annot <- yinter

  # plot
  p <- ggplot(na_df_srt, aes(x=missPercent,y=FeatureCount)) +
    geom_line(color="#F8766D", linewidth=1.5) +
    annotate(
      "text",
      x = 0.6,
      y = nrow(df) / 5,
      label = paste0(
        "Complete data: ",
        nrow(na.omit(df)),
        "\nLess than 50% missing: ",
        yinter,
        "\nAll identified: ",
        nrow(df)
      ),
      size = 4,
      color = "black"
    ) +
    labs(
      x = "Missing percentage",
      y = "Cumulative feature count",
      title = title
    ) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10, face = "bold"),
      plot.title = element_text(size = 12)
    )

  return(p)

}

