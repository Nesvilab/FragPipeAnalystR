# ── Feature-wise correlation ──────────────────────────────────────────────────

#' Calculate feature-wise correlation between two SummarizedExperiment objects
#'
#' @param se1 A \code{SummarizedExperiment} object.
#' @param se2 A \code{SummarizedExperiment} object.
#' @param assay1 Assay index or name for \code{se1}. Default is \code{1}.
#' @param assay2 Assay index or name for \code{se2}. Default is \code{1}.
#' @param feature_id_col Column name in \code{rowData} to use as feature IDs.
#'   If \code{NULL}, \code{rownames} are used.
#' @param use Character string passed to \code{cor}. Default is
#'   \code{"pairwise.complete.obs"}.
#' @param cor_method Correlation method passed to \code{cor}. Default is
#'   \code{"spearman"}.
#' @param missing_value_cutoff Fraction of missing values above which a feature
#'   is excluded. Default is \code{0.5}.
#' @return A data frame with columns \code{ID}, \code{sd_x}, \code{sd_y},
#'   \code{cor}, and \code{pvalue}.
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom tidyr gather
#' @importFrom dplyr group_by summarise filter
#' @export
calculate_feature_wise_correlation_SE <- function(
    se1, se2, assay1 = 1, assay2 = 1,
    feature_id_col = NULL,
    use = "pairwise.complete.obs",
    cor_method = "spearman",
    missing_value_cutoff = 0.5) {

  mat1 <- assay(se1, assay1)
  mat2 <- assay(se2, assay2)
  samples1 <- colnames(mat1)
  samples2 <- colnames(mat2)

  ids1 <- if (is.null(feature_id_col)) rownames(mat1) else rowData(se1)[[feature_id_col]]
  ids2 <- if (is.null(feature_id_col)) rownames(mat2) else rowData(se2)[[feature_id_col]]

  common_ids     <- intersect(ids1, ids2)
  print(length(common_ids))
  common_samples <- intersect(samples1, samples2)
  print(length(common_samples))
  mat1 <- mat1[match(common_ids, ids1), common_samples, drop = FALSE]
  mat2 <- mat2[match(common_ids, ids2), common_samples, drop = FALSE]
  ids  <- common_ids

  para1 <- data.frame(ID = ids, mat1, row.names = NULL, stringsAsFactors = FALSE)
  para2 <- data.frame(ID = ids, mat2, row.names = NULL, stringsAsFactors = FALSE)
  para1 <- para1[rowSums(is.na(para1[, -1])) < ((ncol(para1) - 1) * missing_value_cutoff), ]
  para2 <- para2[rowSums(is.na(para2[, -1])) < ((ncol(para2) - 1) * missing_value_cutoff), ]
  ids_filtered <- intersect(para1$ID, para2$ID)
  para1 <- para1[para1$ID %in% ids_filtered, ]
  para2 <- para2[para2$ID %in% ids_filtered, ]

  x1 <- tidyr::gather(para1, key = "sample", value = "value", -ID)
  x2 <- tidyr::gather(para2, key = "sample", value = "value", -ID)
  m  <- merge(x1, x2, by = c("ID", "sample"))
  names(m)[3:4] <- c("x", "y")
  m  <- m[!apply(m[, 3:4], 1, function(z) any(is.na(z))), ]
  n3 <- m %>% group_by(ID) %>%
    summarise(n = sum(!is.na(x) & !is.na(y))) %>% filter(n >= 3)
  m  <- m %>% filter(ID %in% n3$ID)

  get_cor_pvalue <- function(x, y, cor_method = "spearman", use = "pairwise.complete.obs") {
    if (sum(!is.na(x) & !is.na(y)) < 3) return(NA)
    cor.test(x, y, method = cor_method, use = use)$p.value
  }
  m %>% group_by(ID) %>%
    summarise(
      sd_x   = sd(x, na.rm = TRUE),
      sd_y   = sd(y, na.rm = TRUE),
      cor    = cor(x, y, use = use, method = cor_method),
      pvalue = get_cor_pvalue(x, y, cor_method = cor_method, use = use)
    )
}

# ── Sample-wise correlation ───────────────────────────────────────────────────

#' Calculate sample-wise correlation between two SummarizedExperiment objects
#'
#' @param se1 A \code{SummarizedExperiment} object.
#' @param se2 A \code{SummarizedExperiment} object.
#' @param feature_id_col Column name in \code{rowData} to use as feature IDs.
#'   If \code{NULL}, \code{rownames} are used.
#' @param use Character string passed to \code{cor}. Default is
#'   \code{"pairwise.complete.obs"}.
#' @param cor_method Correlation method passed to \code{cor}. Default is
#'   \code{"spearman"}.
#' @param missing_value_cutoff Fraction of missing values above which a feature
#'   is excluded. Default is \code{0.5}.
#' @return A data frame with columns \code{sample} and \code{correlation}.
#' @importFrom SummarizedExperiment assay rowData
#' @export
calculate_sample_wise_correlation_SE <- function(
    se1, se2,
    feature_id_col = NULL,
    use = "pairwise.complete.obs",
    cor_method = "spearman",
    missing_value_cutoff = 0.5) {

  mat1 <- assay(se1)
  mat2 <- assay(se2)

  ids1 <- if (is.null(feature_id_col)) rownames(mat1) else rowData(se1)[[feature_id_col]]
  ids2 <- if (is.null(feature_id_col)) rownames(mat2) else rowData(se2)[[feature_id_col]]
  common_ids     <- intersect(ids1, ids2)
  common_samples <- intersect(colnames(mat1), colnames(mat2))
  mat1 <- mat1[match(common_ids, ids1), common_samples, drop = FALSE]
  mat2 <- mat2[match(common_ids, ids2), common_samples, drop = FALSE]

  keep1 <- rowSums(is.na(mat1)) < (ncol(mat1) * missing_value_cutoff)
  keep2 <- rowSums(is.na(mat2)) < (ncol(mat2) * missing_value_cutoff)
  mat1  <- mat1[keep1 & keep2, , drop = FALSE]
  mat2  <- mat2[keep1 & keep2, , drop = FALSE]

  data.frame(
    sample      = common_samples,
    correlation = sapply(common_samples, function(s)
      cor(mat1[, s], mat2[, s], use = use, method = cor_method))
  )
}

# ── Plot: feature-wise correlation histogram ──────────────────────────────────

#' Plot feature-wise correlation as a histogram
#'
#' @param feature_cor_df Data frame returned by
#'   \code{\link{calculate_feature_wise_correlation_SE}}.
#' @param alpha Significance threshold for BH-adjusted p-values. Default is
#'   \code{0.01}.
#' @return A \code{ggplot} object.
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline annotate xlab ylab
#'   theme scale_fill_manual
#' @importFrom ggpubr theme_pubr
#' @importFrom dplyr mutate
#' @export
plot_feature_wise_correlation <- function(feature_cor_df, alpha = 0.01) {
  x <- feature_cor_df
  x$pvalue <- p.adjust(x$pvalue, method = "BH")
  br <- c(seq(min(x$cor, na.rm = TRUE) - 0.1, -1e-9, by = 0.01),
          seq(0, max(x$cor, na.rm = TRUE), by = 0.01))
  max_count             <- max(table(cut(x$cor, breaks = br)))
  positive_cor_ratio    <- sum(x$cor > 0, na.rm = TRUE) / nrow(x)
  sig_positive_cor_ratio <- sum(x$pvalue <= alpha, na.rm = TRUE) / nrow(x)
  mean_cor              <- mean(x$cor, na.rm = TRUE)
  n_pairs               <- nrow(x)

  x <- x %>% mutate(col = ifelse(cor > 0, "red", "blue"))
  ggplot(x, aes(x = cor, fill = col)) +
    geom_histogram(breaks = br, color = "black", size = 0.1) +
    xlab("Spearman's correlation") + ylab("Frequency") +
    ggpubr::theme_pubr() +
    theme(legend.position = "none") +
    scale_fill_manual(breaks = c("red", "blue"), values = c("red", "blue")) +
    geom_vline(xintercept = mean_cor, linetype = 2) +
    annotate("text", x = mean_cor, y = 1.05 * max_count,
             label = paste("Mean =", sprintf("%.2f", mean_cor)),
             hjust = -0.05, size = 3.5) +
    annotate("text", x = min(br) + 0.01, y = 1.05 * max_count,
             label = paste0(n_pairs, " pairs\n",
                            sprintf("%.2f%%", 100 * positive_cor_ratio), " positive correlation\n",
                            sprintf("%.2f%%", 100 * sig_positive_cor_ratio),
                            " significant\npositive correlation\n(adjusted P <= ", alpha, ")"),
             vjust = 1, hjust = 0, size = 3.5)
}

# ── Plot: feature-wise correlation density ────────────────────────────────────

#' Plot feature-wise correlation as a density plot (two-color)
#'
#' @param feature_cor_df Data frame returned by
#'   \code{\link{calculate_feature_wise_correlation_SE}}.
#' @param alpha Significance threshold for BH-adjusted p-values. Default is
#'   \code{0.01}.
#' @return A \code{ggplot} object.
#' @importFrom ggplot2 ggplot aes geom_density geom_vline annotate xlab ylab
#'   theme scale_fill_manual
#' @importFrom ggpubr theme_pubr
#' @importFrom dplyr mutate
#' @export
plot_feature_wise_correlation_density <- function(feature_cor_df, alpha = 0.01) {
  x <- feature_cor_df
  x$pvalue <- p.adjust(x$pvalue, method = "BH")
  positive_cor_ratio     <- sum(x$cor > 0, na.rm = TRUE) / nrow(x)
  sig_positive_cor_ratio <- sum(x$pvalue <= alpha, na.rm = TRUE) / nrow(x)
  mean_cor               <- mean(x$cor, na.rm = TRUE)
  n_pairs                <- nrow(x)
  max_density            <- max(density(x$cor, na.rm = TRUE)$y)

  x <- x %>% mutate(col = ifelse(cor > 0, "red", "blue"))
  ggplot(x, aes(x = cor, fill = col)) +
    geom_density(alpha = 0.4, color = NA) +
    xlab("Spearman's correlation") + ylab("Density") +
    ggpubr::theme_pubr() +
    theme(legend.position = "none") +
    scale_fill_manual(breaks = c("red", "blue"), values = c("red", "blue")) +
    geom_vline(xintercept = mean_cor, linetype = 2) +
    annotate("text", x = mean_cor, y = max_density * 1.05,
             label = paste("Mean =", sprintf("%.2f", mean_cor)),
             hjust = -0.05, size = 3.5) +
    annotate("text", x = min(x$cor, na.rm = TRUE) + 0.01, y = max_density * 1.05,
             label = paste0(n_pairs, " pairs\n",
                            sprintf("%.2f%%", 100 * positive_cor_ratio), " positive correlation\n",
                            sprintf("%.2f%%", 100 * sig_positive_cor_ratio),
                            " significant\npositive correlation\n(adjusted P <= ", alpha, ")"),
             vjust = 1, hjust = 0, size = 3.5)
}

# ── Plot: feature-wise correlation density (single color) ────────────────────

#' Plot feature-wise correlation as a density plot (single color)
#'
#' @param feature_cor_df Data frame returned by
#'   \code{\link{calculate_feature_wise_correlation_SE}}.
#' @param alpha Significance threshold for BH-adjusted p-values. Default is
#'   \code{0.01}.
#' @return A \code{ggplot} object.
#' @importFrom ggplot2 ggplot aes geom_density geom_vline annotate xlab ylab
#' @importFrom ggpubr theme_pubr
#' @export
plot_feature_correlation_distribution <- function(feature_cor_df, alpha = 0.01) {
  feature_cor_df$pvalue <- p.adjust(feature_cor_df$pvalue, method = "BH")
  mean_cor               <- mean(feature_cor_df$cor, na.rm = TRUE)
  n_pairs                <- nrow(feature_cor_df)
  positive_cor_ratio     <- sum(feature_cor_df$cor > 0, na.rm = TRUE) / n_pairs
  sig_positive_cor_ratio <- sum(feature_cor_df$pvalue <= alpha, na.rm = TRUE) / n_pairs
  max_density            <- max(density(feature_cor_df$cor, na.rm = TRUE)$y)

  ggplot(feature_cor_df, aes(x = cor)) +
    geom_density(fill = "skyblue", alpha = 0.5) +
    geom_vline(xintercept = mean_cor, linetype = 2, color = "darkred") +
    xlab("Spearman's correlation") + ylab("Density") +
    ggpubr::theme_pubr() +
    annotate("text", x = mean_cor, y = max_density * 1.05,
             label = paste("Mean =", sprintf("%.2f", mean_cor)),
             hjust = -0.05, size = 3.5, color = "darkred") +
    annotate("text", x = min(feature_cor_df$cor, na.rm = TRUE) + 0.01, y = max_density * 1.05,
             label = paste0(n_pairs, " pairs\n",
                            sprintf("%.2f%%", 100 * positive_cor_ratio), " positive correlation\n",
                            sprintf("%.2f%%", 100 * sig_positive_cor_ratio),
                            " significant positive correlation\n(adjusted P ≤ ", alpha, ")"),
             vjust = 1, hjust = 0, size = 3.5)
}

# ── Plot: feature-wise correlation density + histogram combo ─────────────────

#' Plot feature-wise correlation as a combined density and histogram
#'
#' @param feature_cor_df Data frame returned by
#'   \code{\link{calculate_feature_wise_correlation_SE}}.
#' @param alpha Significance threshold for BH-adjusted p-values. Default is
#'   \code{0.01}.
#' @return A \code{ggplot} object.
#' @importFrom ggplot2 ggplot aes geom_histogram geom_density geom_vline
#'   annotate xlab ylab
#' @importFrom ggpubr theme_pubr
#' @export
plot_feature_correlation_distribution_combo <- function(feature_cor_df, alpha = 0.01) {
  feature_cor_df$pvalue <- p.adjust(feature_cor_df$pvalue, method = "BH")
  median_cor             <- median(feature_cor_df$cor, na.rm = TRUE)
  n_pairs                <- nrow(feature_cor_df)
  positive_cor_ratio     <- sum(feature_cor_df$cor > 0, na.rm = TRUE) / n_pairs
  sig_positive_cor_ratio <- sum(feature_cor_df$pvalue <= alpha, na.rm = TRUE) / n_pairs
  max_density            <- max(density(feature_cor_df$cor, na.rm = TRUE)$y)

  ggplot(feature_cor_df, aes(x = cor)) +
    geom_histogram(aes(y = ..density..), binwidth = 0.02,
                   fill = "gray90", color = "gray60", alpha = 0.6) +
    geom_density(fill = "skyblue", alpha = 0.3, color = "blue", size = 1) +
    geom_vline(xintercept = median_cor, linetype = 2, color = "darkred") +
    xlab("Spearman's correlation") + ylab("Density") +
    ggpubr::theme_pubr() +
    annotate("text", x = median_cor, y = max_density * 1.05,
             label = paste("Median =", sprintf("%.2f", median_cor)),
             hjust = -0.05, size = 3.5, color = "darkred") +
    annotate("text", x = min(feature_cor_df$cor, na.rm = TRUE) + 0.01, y = max_density * 1.05,
             label = paste0(n_pairs, " pairs\n",
                            sprintf("%.2f%%", 100 * positive_cor_ratio), " positive correlation\n",
                            sprintf("%.2f%%", 100 * sig_positive_cor_ratio),
                            " significant positive correlation\n(adjusted P ≤ ", alpha, ")"),
             vjust = 1, hjust = 0, size = 3.5)
}

# ── Plot: sample-wise correlation dotplot ────────────────────────────────────

#' Plot sample-wise correlation as a dot plot
#'
#' @param sample_cor_df Data frame returned by
#'   \code{\link{calculate_sample_wise_correlation_SE}}.
#' @return A \code{ggplot} object.
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab ggtitle theme_bw theme
#'   element_text
#' @importFrom stats reorder
#' @importFrom dplyr arrange
#' @export
plot_sample_wise_correlation <- function(sample_cor_df) {
  cor_mean   <- mean(sample_cor_df$correlation, na.rm = TRUE)
  cor_median <- median(sample_cor_df$correlation, na.rm = TRUE)

  ggplot(sample_cor_df %>% arrange(correlation),
         aes(x = reorder(sample, correlation), y = correlation)) +
    geom_point(color = "darkblue", size = 2) +
    xlab("Sample (ordered)") +
    ylab("Sample-wise correlation") +
    ggtitle(paste("Mean:", sprintf("%.2f", cor_mean),
                  " | Median:", sprintf("%.2f", cor_median))) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ── Plot: Venn diagram of features across SE objects ─────────────────────────

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
#' @importFrom S4Vectors metadata
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
