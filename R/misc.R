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
plot_feature <- function(se, protein, type="boxplot", id="label", interactive=F) {
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.character(protein),
                          is.character(type))
  subset <- se[protein]

  df_reps <- data.frame(assay(subset), check.names = F) %>%
    rownames_to_column() %>%
    gather(ID, val, -rowname) %>%
    left_join(., data.frame(colData(subset)), by = c("ID"=id))

  df_reps$rowname <- parse_factor(as.character(df_reps$rowname), levels = protein)

  df_CI<- df_reps %>%
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
          theme(axis.title.x = element_blank()) +
          theme_bw()
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
          theme(axis.title.x = element_blank()) +
          theme_bw()
      }
    } else if(type=="boxplot") {
      if (max(df_reps$replicate) == 1){
        p <- ggplot(df_reps, aes(condition, val))+
          geom_boxplot()+
          geom_jitter(size = 3, position = position_dodge(width=0.3)) +
          labs(y = expression(log[2]~"Intensity")) +
          facet_wrap(~rowname) +
          theme(axis.title.x = element_blank()) +
          theme_bw()
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
          theme(axis.title.x = element_blank()) +
          theme_bw()
      }
    }
  }
  return(p)
}
