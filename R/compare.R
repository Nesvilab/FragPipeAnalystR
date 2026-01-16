compare_CumulativeMissingPercent <- function(se1, se2,
                                             smp_lst1 = NULL, smp_lst2 = NULL,
                                             se1_name = NULL, se2_name = NULL,
                                             title = "") {
  # Helper function to get cumulative missing percent dataframe
  get_missing_df <- function(df, smp_lst) {
    if (is.null(smp_lst)) {
      smp_lst <- colnames(df)
    } else {
      smp_lst <- intersect(smp_lst, colnames(df))
    }
    df <- data.table::setDT(df)
    na_df <- df[, .(Index = rownames(df), NA_count = rowSums(is.na(.SD))), .SDcols = smp_lst]
    na_df_srt <- na_df %>%
      dplyr::arrange(NA_count) %>%
      dplyr::mutate(missPercent = NA_count / length(smp_lst),
                    FeatureCount = seq_len(nrow(.)))
    return(na_df_srt)
  }
  df1 <- data.frame(assay(se1))
  df2 <- data.frame(assay(se2))
  if (is.null(se1_name)) {
    se1_name <- metadata(se1)$exp
  }
  if (is.null(se2_name)) {
    se2_name <- metadata(se2)$exp
  }
  na_df_srt1 <- get_missing_df(df1, smp_lst1) %>% dplyr::mutate(Source = se1_name)
  na_df_srt2 <- get_missing_df(df2, smp_lst2) %>% dplyr::mutate(Source = se2_name)
  # Combine for plotting
  plot_df <- dplyr::bind_rows(na_df_srt1, na_df_srt2)
  # Plot
  p <- ggplot(plot_df, aes(x=missPercent, y=FeatureCount, color=Source)) +
    geom_line(linewidth=1.5) +
    labs(
      x = "Missing percentage",
      y = "Cumulative feature count",
      title = title,
      color = "Dataset"
    ) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10, face = "bold"),
      plot.title = element_text(size = 12)
    )
  return(p)
}

