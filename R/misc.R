#' @export
plot_cvs <- function(se, id="ID", scale=T, check.names=T) {

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

#' @export
plot_feature <- function(se, protein, index=NULL,
                         type="boxplot", id="label", interactive=F) {
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.character(protein),
                          is.character(type))
  if (is.null(index)) {
    subset <- se[protein]
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

#' @export
plot_feature_numbers <- function(se, exp=NULL, feature="Proteins") {
  # Show error if input is not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))

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
      left_join(., data.frame(colData(se)), by = c("ID"="label"))

    p <- ggplot(stat, aes(x = sample_name, y = sum, fill = condition)) +
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
