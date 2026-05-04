#' Plot missing value heatmap
#'
#' Creates a binary heatmap visualization showing the pattern of missing values
#' across samples and proteins.
#'
#' @param se A \code{SummarizedExperiment} object containing proteomics data
#'   with missing values.
#'
#' @return A ComplexHeatmap object showing the missing value pattern where
#'   white cells represent missing values and black cells represent valid values.
#'   Only proteins with at least one missing value are displayed.
#'
#' @details
#' This visualization helps identify whether missing values are randomly
#' distributed or show systematic patterns (e.g., all missing in certain
#' samples or conditions), which can inform imputation strategy selection.
#'
#' @examples
#' \dontrun{
#' # Plot missing value pattern
#' plot_missval_heatmap(se)
#' }
#'
#' @seealso \code{\link{plotCumulativeMissingPercent}}, \code{\link{impute}}
#'
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom SummarizedExperiment assay
#' @importFrom grid gpar
#'
#' @export
plot_missval_heatmap <- function(se) {
  # Show error if input is not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))

  se_assay <- assay(se)
  # Show error if there are no missing values
  if (!any(is.na(se_assay))) {
    stop("No missing values in '", deparse(substitute(se)), "'",
      call. = FALSE
    )
  }

  # Make assay data binary (1 = valid value, 0 = missing value)
  df <- se_assay %>% data.frame(., check.names = F)

  missval <- df[apply(df, 1, function(x) any(is.na(x))), ]
  missval <- ifelse(is.na(missval), 0, 1)

  # Plot binary heatmap
  ht2 <- Heatmap(missval,
    col = c("white", "black"),
    column_names_side = "top",
    show_row_names = FALSE,
    show_column_names = TRUE,
    name = paste0("Missing values pattern (", dim(missval)[1], " proteins )"),
    column_names_gp = gpar(fontsize = 16),
    heatmap_legend_param = list(
      at = c(0, 1),
      labels = c("Missing value", "Valid value")
    )
  )
  return(draw(ht2, heatmap_legend_side = "top"))
}

#' Plot sample correlation heatmap
#'
#' Creates a heatmap visualization of Pearson correlation coefficients between
#' all sample pairs, useful for quality control and identifying outliers.
#'
#' @param dep A \code{SummarizedExperiment} object containing proteomics data.
#' @param significant Logical indicating whether to use only significant proteins
#'   for correlation calculation. Requires a "significant" column from
#'   \code{\link{add_rejections}}. Default is \code{FALSE}.
#' @param lower Numeric value (-1 to 1) for the lower bound of the color scale.
#'   Default is -1.
#' @param upper Numeric value (-1 to 1) for the upper bound of the color scale.
#'   Default is 1.
#' @param pal Character string specifying the RColorBrewer palette name.
#'   Must be a sequential or diverging palette (not qualitative).
#'   Default is "PRGn".
#' @param pal_rev Logical indicating whether to reverse the color palette.
#'   Default is \code{FALSE}.
#' @param indicate Character vector specifying colData columns to use for
#'   heatmap annotation. Default is \code{NULL}.
#' @param font_size Numeric value for sample name font size. Default is 10.
#' @param exp Character string for experiment type. If \code{NULL}, uses
#'   metadata value. Default is \code{NULL}.
#' @param use Character string specifying how to handle missing values in
#'   correlation calculation. Passed to \code{cor()}. Default is "complete.obs".
#' @param ... Additional arguments passed to \code{ComplexHeatmap::Heatmap()}.
#'
#' @return A ComplexHeatmap object showing the sample correlation matrix
#'   with hierarchical clustering.
#'
#' @details
#' High correlations between biological replicates and distinct clustering
#' of different conditions indicate good data quality. Samples that don't
#' cluster with their replicates may be outliers.
#'
#' @examples
#' \dontrun{
#' # Basic correlation heatmap
#' plot_correlation_heatmap(se)
#'
#' # With condition annotation
#' plot_correlation_heatmap(se, indicate = "condition")
#'
#' # Using only significant proteins
#' plot_correlation_heatmap(se_diff, significant = TRUE)
#' }
#'
#' @seealso \code{\link{plot_pca}}, \code{\link{get_cluster_heatmap}}
#'
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom SummarizedExperiment assay colData rowData metadata
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar
#' @importFrom dplyr select all_of
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr replace_na
#' @importFrom stats cor
#'
#' @export
plot_correlation_heatmap <- function(dep, significant = FALSE, lower = -1, upper = 1,
                                     pal = "PRGn", pal_rev = FALSE, indicate = NULL,
                                     font_size = 10, exp = NULL, use="complete.obs",...) {
  if (is.null(exp)) {
    exp <- metadata(dep)$exp
  }
  # Show error if inputs are not the required classes
  assertthat::assert_that(
    inherits(dep, "SummarizedExperiment"),
    is.logical(significant),
    length(significant) == 1,
    is.numeric(lower),
    length(lower) == 1,
    is.numeric(upper),
    length(upper) == 1,
    is.character(pal),
    length(pal) == 1,
    is.logical(pal_rev),
    length(pal_rev) == 1,
    is.numeric(font_size),
    length(font_size) == 1
  )

  # Check for valid lower and upper values
  if (!(lower >= -1 & upper >= -1 & lower <= 1 & upper <= 1)) {
    stop("'lower' and/or 'upper' arguments are not valid
         Run plot_pca() with 'lower' and 'upper' between -1 and 1",
         call. = FALSE
    )
  }

  # Check for valid pal
  pals <- brewer.pal.info %>%
    rownames_to_column() %>%
    filter(category != "qual")
  if (!pal %in% pals$rowname) {
    stop("'", pal, "' is not a valid color panel",
         " (qualitative panels also not allowed)\n",
         "Run plot_pca() with one of the following 'pal' options: ",
         paste(pals$rowname, collapse = "', '"), "'",
         call. = FALSE
    )
  }

  # if(any(is.na(assay(dep)))) {
  #   stop("Missing values in '", deparse(substitute(dep)), "'. Use plot_dist() instead")
  # }

  # Heatmap annotation
  if (!is.null(indicate)) {
    assertthat::assert_that(is.character(indicate))

    col_data <- colData(dep) %>%
      as.data.frame()
    columns <- colnames(col_data)
    if (any(!indicate %in% columns)) {
      stop("'",
           paste0(indicate, collapse = "' and/or '"),
           "' column(s) is/are not present in ",
           deparse(substitute(dep)),
           ".\nValid columns are: '",
           paste(columns, collapse = "', '"),
           "'.",
           call. = FALSE
      )
    }

    # Get annotation
    anno <- colData(dep) %>%
      data.frame() %>%
      select(all_of(indicate))

    # Annotation color
    names <- colnames(anno)
    anno_col <- vector(mode = "list", length = length(names))
    names(anno_col) <- names
    for (i in names) {
      var <- anno[[i]] %>%
        unique() %>%
        sort()
      if (length(var) == 1) {
        cols <- c("black")
      } else if (length(var) == 2) {
        cols <- c("orangered", "cornflowerblue")
      } else if (length(var) < 7 & length(var) > 2) {
        cols <- brewer.pal(length(var), "Pastel1")
      } else if (length(var) <= 12) {
        cols <- brewer.pal(length(var), "Set3")
      } else {
        cols <- colorRampPalette(brewer.pal(12, "Set3"))(length(var))
      }
      names(cols) <- var
      anno_col[[i]] <- cols
    }

    # HeatmapAnnotation object
    ha1 <- HeatmapAnnotation(
      df = anno,
      col = anno_col,
      show_annotation_name = TRUE
    )
  } else {
    ha1 <- NULL
  }

  # Filter for significant proteins
  if (significant) {
    # Check for significant column
    if (!"significant" %in% colnames(rowData(dep, use.names = FALSE))) {
      stop("'significant' column is not present in '",
           deparse(substitute(dep)),
           "'\nRun add_rejections() to obtain the required column",
           call. = FALSE
      )
    }

    dep <- dep[replace_na(rowData(dep, use.names = FALSE)$significant, FALSE), ]
  }

  # Calculate correlation matrix
  data <- assay(dep)
  cor_mat <- cor(data, use = use)
  lower <- min(cor_mat)
  upper <- max(cor_mat)

  # Plot heatmap
  ht1 <- Heatmap(cor_mat,
                 col = circlize::colorRamp2(
                   seq(lower, upper, ((upper - lower) / 7)),
                   if (pal_rev) {
                     rev(RColorBrewer::brewer.pal(8, pal))
                   } else {
                     RColorBrewer::brewer.pal(8, pal)
                   }
                 ),
                 heatmap_legend_param = list(
                   color_bar = "continuous",
                   legend_direction = "horizontal",
                   legend_width = unit(5, "cm"),
                   title_position = "topcenter"
                 ),
                 name = "Pearson correlation",
                 column_names_gp = gpar(fontsize = font_size),
                 row_names_gp = gpar(fontsize = font_size),
                 top_annotation = ha1,
                 ...
  )

  return(ht1)
}


#' Clustered heatmap of differentially expressed proteins
#'
#' Creates a heatmap visualization of differentially expressed proteins with
#' hierarchical or k-means clustering, showing either fold changes or
#' centered intensities.
#'
#' @param dep A \code{SummarizedExperiment} object containing differential
#'   expression results from \code{\link{test_limma}} or \code{\link{test_diff}}
#'   with significance annotations from \code{\link{add_rejections}}.
#' @param type Character string specifying the data to visualize:
#'   \itemize{
#'     \item \code{"contrast"}: Log2 fold changes for each contrast
#'     \item \code{"centered"}: Mean-centered log2 intensities per sample
#'   }
#' @param kmeans Logical indicating whether to perform k-means clustering
#'   on proteins. Default is \code{FALSE} (hierarchical clustering).
#' @param k Numeric value for the number of k-means clusters when
#'   \code{kmeans = TRUE}. Default is 6.
#' @param col_limit Numeric value for the color scale limits (-col_limit to
#'   +col_limit). Default is 6.
#' @param indicate Character vector specifying colData columns for heatmap
#'   annotation. Only applicable when \code{type = "centered"}.
#'   Default is \code{NULL}.
#' @param alpha Numeric value for the significance threshold on adjusted
#'   p-values. Default is 0.01.
#' @param lfc Numeric value for the log2 fold change threshold. Default is 1.
#' @param clustering_distance Character string for the distance metric used
#'   in hierarchical clustering. Options include "euclidean", "pearson",
#'   "spearman", "gower" (for data with missing values), etc.
#' @param row_font_size Numeric value for protein name font size. Default is 6.
#' @param col_font_size Numeric value for sample/contrast name font size.
#'   Default is 10.
#' @param plot Logical indicating whether to draw the heatmap. Default is \code{TRUE}.
#' @param ... Additional arguments passed to \code{ComplexHeatmap::Heatmap()}.
#'
#' @return A list containing:
#'   \itemize{
#'     \item The ComplexHeatmap object
#'     \item Row cluster order information
#'   }
#'   If no differentially expressed proteins are found, returns a ggplot
#'   object with a message.
#'
#' @details
#' Only proteins meeting both the \code{alpha} (adjusted p-value) and \code{lfc}
#' (log2 fold change) thresholds are included in the heatmap. When data contains
#' missing values, Gower distance is automatically used for clustering.
#'
#' K-means clusters are ordered by their average fold change magnitude to
#' facilitate interpretation.
#'
#' @examples
#' \dontrun{
#' # Basic cluster heatmap with fold changes
#' get_cluster_heatmap(se_diff, type = "contrast")
#'
#' # Centered intensities with k-means clustering
#' get_cluster_heatmap(se_diff, type = "centered", kmeans = TRUE, k = 4)
#'
#' # With condition annotation
#' get_cluster_heatmap(se_diff, type = "centered", indicate = "condition")
#'
#' # Stricter thresholds
#' get_cluster_heatmap(se_diff, alpha = 0.001, lfc = 2)
#' }
#'
#' @seealso \code{\link{plot_correlation_heatmap}}, \code{\link{test_limma}},
#'   \code{\link{add_rejections}}
#'
#' @importFrom ComplexHeatmap Heatmap draw row_order
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom circlize colorRamp2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom dplyr select mutate group_by summarize arrange pull ends_with
#' @importFrom tidyr gather
#' @importFrom tibble column_to_rownames
#' @importFrom ggplot2 ggplot annotate theme_void
#' @importFrom grid gpar
#' @importFrom stats kmeans
#' @importFrom cluster daisy
#'
#' @export
get_cluster_heatmap <- function(dep, type = c("contrast", "centered"),
                                kmeans = FALSE, k = 6,
                                col_limit = 6, indicate = NULL,
                                alpha = 0.01, lfc = 1,
                                clustering_distance = c("euclidean", "maximum", "manhattan", "canberra",
                                                        "binary", "minkowski", "pearson", "spearman", "kendall", "gower"),
                                row_font_size = 6, col_font_size = 10, plot = TRUE, ...) {

  # Show error if inputs are not the required classes
  if(is.integer(k)) k <- as.numeric(k)
  if(is.integer(col_limit)) col_limit <- as.numeric(col_limit)
  if(is.integer(row_font_size)) row_font_size <- as.numeric(row_font_size)
  if(is.integer(col_font_size)) col_font_size <- as.numeric(col_font_size)
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.character(type),
                          is.logical(kmeans),
                          is.numeric(k),
                          length(k) == 1,
                          is.numeric(col_limit),
                          length(col_limit) == 1,
                          is.numeric(row_font_size),
                          length(row_font_size) == 1,
                          is.numeric(col_font_size),
                          length(col_font_size) == 1,
                          is.logical(plot),
                          length(plot) == 1)

  # Show error if inputs do not contain required columns
  type <- match.arg(type)
  clustering_distance <- match.arg(clustering_distance)

  # Extract row and col data
  row_data <- rowData(dep)
  col_data <- colData(dep) %>%
    as.data.frame()

  # Show error if inputs do not contain required columns
  if(any(!c("label", "condition", "replicate") %in% colnames(col_data))) {
    stop(paste0("'label', 'condition' and/or 'replicate' columns are not present in '",
                deparse(substitute(dep)), "'"),
         call. = FALSE)
  }
  if(length(grep("_diff", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_diff' columns are not present in '",
                deparse(substitute(dep)),
                "'.\nRun test_diff() to obtain the required columns."),
         call. = FALSE)
  }
  if(!"significant" %in% colnames(row_data)) {
    stop(paste0("'significant' column is not present in '",
                deparse(substitute(dep)),
                "'.\nRun add_rejections() to obtain the required column."),
         call. = FALSE)
  }

  # Heatmap annotation
  if(!is.null(indicate) & type == "contrast") {
    warning("Heatmap annotation only applicable for type = 'centered'",
            call. = FALSE)
  }
  if(!is.null(indicate) & type == "centered") {
    ha1 <- get_annotation(dep, indicate)
  } else {
    ha1 <- NULL
  }

  conditions <- gsub("_diff", "",colnames(row_data)[grepl("_diff", colnames(row_data))])
  cols_p <- paste0(conditions, "_p.adj")
  cols_lfc <- paste0(conditions, "_diff")
  p <- as.matrix(row_data[,cols_p]) <= alpha
  lfc <- abs(as.matrix(row_data[,cols_lfc])) >= lfc
  p[is.na(p)] <- FALSE
  lfc[is.na(p)] <- FALSE
  colnames(p) <- conditions
  colnames(lfc) <- conditions
  combined_rejections <- p
  colnames(combined_rejections) <- conditions
  combined_rejections <- p & lfc
  filtered <- dep[apply(combined_rejections, 1, any), ]

  if (nrow(filtered) == 0){
    return(ggplot() +
             annotate("text", x = 4, y = 25, size=8, label = "No differential expressed genes available for the heatmap") +
             theme_void())
  }

  # Check for missing values
  if(any(is.na(assay(filtered)))) {
    warning("Missing values in '", deparse(substitute(dep)), "'. ",
            "Using clustering_distance = 'gower'",
            call. = FALSE)
    clustering_distance <- "gower"
    obs_NA <- TRUE
  } else {
    obs_NA <- FALSE
  }

  # Get centered intensity values ('centered')
  if(type == "centered") {
    rowData(filtered)$mean <- rowMeans(assay(filtered), na.rm = TRUE)
    df <- assay(filtered) - rowData(filtered)$mean
  }
  # Get contrast fold changes ('contrast')
  if(type == "contrast") {
    df <- rowData(filtered) %>%
      data.frame(.) %>%
      column_to_rownames(var = "name") %>%
      select(dplyr::ends_with("_diff"))
    colnames(df) <-
      gsub("_diff", "", colnames(df)) %>%
      gsub("_vs_", " vs ", .)
  }

  # Facultative kmeans clustering
  if(kmeans & obs_NA) {
    warning("Cannot perform kmeans clustering with missing values",
            call. = FALSE)
    kmeans <- FALSE
  }
  if(kmeans & !obs_NA) {
    set.seed(1)
    df_kmeans <- kmeans(df, k)
    if(type == "centered") {
      # Order the k-means clusters according to the maximum fold change
      # in all samples averaged over the proteins in the cluster
      order <- data.frame(df) %>%
        cbind(., cluster = df_kmeans$cluster) %>%
        dplyr::mutate(row = apply(.[, seq_len(ncol(.) - 1)], 1, function(x) max(x))) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarize(index = sum(row)/n()) %>%
        dplyr::arrange(desc(index)) %>%
        dplyr::pull(cluster) %>%
        match(seq_len(k), .)
      df_kmeans$cluster <- order[df_kmeans$cluster]
    }
    if(type == "contrast") {
      # Order the k-means clusters according to their average fold change
      order <- cbind(df, cluster = df_kmeans$cluster) %>%
        dplyr::gather(condition, diff, -cluster) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarize(row = mean(diff)) %>%
        dplyr::arrange(desc(row)) %>%
        dplyr::pull(cluster) %>%
        match(seq_len(k), .)
      df_kmeans$cluster <- order[df_kmeans$cluster]
    }
  }

  if(ncol(df) == 1) {
    col_clust = FALSE
  } else {
    col_clust = TRUE
  }
  if(nrow(df) == 1) {
    row_clust = FALSE
  } else {
    row_clust = TRUE
  }
  if(clustering_distance == "gower") {
    clustering_distance <- function(x) {
      dist <- cluster::daisy(x, metric = "gower")
      dist[is.na(dist)] <- max(dist, na.rm = TRUE)
      return(dist)
    }
  }

  # Legend info
  legend <- ifelse(type == "contrast",
                   "log2 Fold change",
                   "log2 Centered intensity")

  # use sample name for the heatmap
  temp <- colData(filtered)
  rownames(temp) <- temp$label
  colnames(df) <- temp[colnames(df), "sample_name"]

  # Heatmap
  ht1 = Heatmap(df,
                col = circlize::colorRamp2(
                  seq(-col_limit, col_limit, (col_limit/5)),
                  rev(RColorBrewer::brewer.pal(11, "RdBu"))),
                split = if(kmeans) {df_kmeans$cluster} else {NULL},
                cluster_rows = col_clust,
                cluster_columns = row_clust,
                row_names_side = "left",
                column_names_side = "top",
                clustering_distance_rows = clustering_distance,
                clustering_distance_columns = clustering_distance,
                heatmap_legend_param = list(color_bar = "continuous",
                                            legend_direction = "horizontal",
                                            legend_width = unit(5, "cm"),
                                            title_position = "lefttop"),
                name = legend,
                row_names_gp = gpar(fontsize = row_font_size),
                column_names_gp = gpar(fontsize = col_font_size),
                top_annotation = ha1,
                ...)
  # return (row_order(ht1))
  # Return data.frame
  p <- draw(ht1, heatmap_legend_side = "top")
  row_clusters <- row_order(ht1)
  #mat<-as.matrix(df)

  # for (i in 1:length(row_clusters)){
  #   if (i==1){
  #     clu <-t(t(row.names(ht1[row_clusters[[i]],])))
  #     out <-cbind (clu, paste("cluster", i, sep=""))
  #     colnames(out)<- c("ProteinID", "Cluster")
  #   }
  #   else{
  #     clu <- t(t(row.names(ht1[row_clusters[[i]],])))
  #     clu <- cbind(clu, paste("cluster", i, sep = ""))
  #     out <- cbind(out, clu)
  #   }
  # }
  heatmap_list <- list(ht1, row_clusters)
  return(heatmap_list)
}

# Internal function to get ComplexHeatmap::HeatmapAnnotation object
get_annotation <- function(dep, indicate) {
  assertthat::assert_that(
    inherits(dep, "SummarizedExperiment"),
    is.character(indicate))

  # Check indicate columns
  col_data <- colData(dep) %>%
    as.data.frame()
  columns <- colnames(col_data)
  if(all(!indicate %in% columns)) {
    stop("'",
         paste0(indicate, collapse = "' and/or '"),
         "' column(s) is/are not present in ",
         deparse(substitute(dep)),
         ".\nValid columns are: '",
         paste(columns, collapse = "', '"),
         "'.",
         call. = FALSE)
  }
  if(any(!indicate %in% columns)) {
    indicate <- indicate[indicate %in% columns]
    warning("Only used the following indicate column(s): '",
            paste0(indicate, collapse = "', '"),
            "'")
  }

  # Get annotation
  anno <- dplyr::select(col_data, indicate)

  # Annotation color
  names <- colnames(anno)
  anno_col <- vector(mode="list", length=length(names))
  names(anno_col) <- names
  for(i in names) {
    var = anno[[i]] %>% unique() %>% sort()
    if(length(var) == 1) {
      cols <- c("black")
    } else if(length(var) == 2) {
      cols <- c("orangered", "cornflowerblue")
    } else if(length(var) < 7 & length(var) > 2) {
      cols <- RColorBrewer::brewer.pal(length(var), "Pastel1")
    } else if(length(var) <= 12) {
      cols <- RColorBrewer::brewer.pal(length(var), "Set3")
    } else {
      cols <- colorRampPalette(brewer.pal(12, "Set3"))(length(var))
    }
    names(cols) <- var
    anno_col[[i]] <-  cols
  }

  # HeatmapAnnotation object
  ComplexHeatmap::HeatmapAnnotation(df = anno,
                                    col = anno_col,
                                    show_annotation_name = TRUE)
}

#' Plot replicate correlation half-heatmap
#'
#' Creates a triangular heatmap visualization showing Pearson correlation
#' coefficients between specified replicate samples, with correlation values
#' displayed in each cell.
#'
#' @param se A \code{SummarizedExperiment} object containing proteomics data.
#' @param rep_samples Character vector of sample names (column names from assay)
#'   to include in the correlation matrix. Typically technical or biological
#'   replicates that should show high correlation.
#'
#' @return A ggplot object (converted from ComplexHeatmap) showing a triangular
#'   correlation matrix with:
#'   \itemize{
#'     \item Color scale from blue (-1) through white (0) to red (+1)
#'     \item Correlation coefficients displayed in lower triangle
#'     \item No clustering (samples in specified order)
#'   }
#'
#' @details
#' This visualization is particularly useful for quality control of technical
#' replicates, where correlations should typically be > 0.95. The triangular
#' format reduces visual redundancy compared to full correlation matrices.
#'
#' @examples
#' \dontrun{
#' # Plot correlation between specific replicates
#' plot_replicate_heatmap(se, rep_samples = c("Sample1_Rep1", "Sample1_Rep2",
#'                                            "Sample1_Rep3"))
#'
#' # All samples from one condition
#' condition_samples <- colData(se)$sample_name[colData(se)$condition == "Control"]
#' plot_replicate_heatmap(se, rep_samples = condition_samples)
#' }
#'
#' @seealso \code{\link{plot_correlation_heatmap}}
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom SummarizedExperiment assay
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar grid.rect grid.text
#' @importFrom ggplotify as.ggplot
#' @importFrom stats cor
#'
#' @export
plot_replicate_heatmap <- function(se,
                                   rep_samples){
  d <- assay(se)
  samples <- rep_samples
  corr_mat <- round(cor(d[,samples], use = "pairwise.complete.obs"), 2)
  x_samples <- rep_samples
  y_samples <- rep_samples
  corr_mat <- corr_mat[x_samples, y_samples]

  col_fun = circlize::colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))
  ht <- Heatmap(
    corr_mat,
    col = col_fun,
    cluster_rows = F,
    cluster_columns = F,
    rect_gp = grid::gpar(type = "none"),
    show_heatmap_legend = T,
    show_column_names = T,
    show_row_names = T,
    row_dend_side = "right",
    row_names_side = "left",
    cell_fun = function(j, i, x, y, w, h, fill) {
      if (as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
        grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = fill, col = fill))
      }
      if (i >= j) {
        grid::grid.text(sprintf("%.2f", corr_mat[i, j]), x, y, gp = grid::gpar(fontsize = 6.5))
      }
    },
    name = "correlation"
  )
  gg_heatmap <- ggplotify::as.ggplot(ht)
  return(gg_heatmap)
}

#' ORA Pathway Heatmap across all contrasts
#'
#' Takes the data frame returned by \code{or_test()} and draws a
#' \code{ComplexHeatmap} where rows = pathway terms and columns = contrasts.
#'
#' @param ora_results data.frame from \code{or_test()}.
#' @param value_type \code{"log2OR"} (log2 odds ratio), \code{"pvalue"}
#'   (-log10 hypergeometric p), or \code{"size"} (hit count).
#' @param alpha Significance cutoff — only terms significant in at least one
#'   contrast are shown.
#' @param use_adjp Use adjusted p-value for the significance filter?
#' @return A \code{ComplexHeatmap} object.
#' @seealso \code{\link{or_test}}
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar unit
#' @importFrom tidyr pivot_wider
#' @export
plot_ora_heatmap <- function(ora_results, value_type = "log2OR",
                             alpha = 0.05, use_adjp = TRUE) {
  required_cols <- c("Term", "contrast", "log_odds", "p_hyper", "IN")
  missing <- setdiff(required_cols, colnames(ora_results))
  if (length(missing) > 0)
    stop("ORA results missing columns: ", paste(missing, collapse = ", "))

  has_direction <- "direction" %in% colnames(ora_results)
  if (has_direction)
    ora_results <- ora_results %>%
      dplyr::mutate(Term = paste0(Term, " (", direction, ")"))

  p_col <- if (use_adjp && "p.adjust_hyper" %in% colnames(ora_results))
    "p.adjust_hyper" else "p_hyper"

  sig_terms <- ora_results %>%
    dplyr::filter(.data[[p_col]] < alpha) %>%
    dplyr::pull(Term) %>%
    unique()

  if (length(sig_terms) == 0)
    stop("No ORA terms with ", if (use_adjp) "adj. " else "", "p < ", alpha,
         " in any contrast.\nTry raising the cutoff or choosing a different database.")

  df_sub <- ora_results %>% dplyr::filter(Term %in% sig_terms)

  if (value_type == "log2OR") {
    fill_label <- "log2 OR"; fill_na <- 0
    df_wide    <- df_sub %>%
      dplyr::select(Term, contrast, val = log_odds) %>%
      tidyr::pivot_wider(names_from = contrast, values_from = val,
                         values_fill = fill_na)
  } else if (value_type == "pvalue") {
    fill_label <- "-log10(p)"; fill_na <- 0
    df_wide    <- df_sub %>%
      dplyr::mutate(val = -log10(pmax(p_hyper, 1e-300))) %>%
      dplyr::select(Term, contrast, val) %>%
      tidyr::pivot_wider(names_from = contrast, values_from = val,
                         values_fill = fill_na)
  } else {
    fill_label <- "Hits (size)"; fill_na <- 0L
    df_wide    <- df_sub %>%
      dplyr::select(Term, contrast, val = IN) %>%
      tidyr::pivot_wider(names_from = contrast, values_from = val,
                         values_fill = fill_na)
  }

  mat           <- as.matrix(df_wide[, -1, drop = FALSE])
  rownames(mat) <- df_wide$Term
  rownames(mat) <- gsub(
    "^HALLMARK_|^KEGG_|^REACTOME_|^WP_|^GOBP_|^GOMF_|^GOCC_",
    "", rownames(mat))
  rownames(mat) <- gsub("_", " ", rownames(mat))

  if (value_type == "log2OR") {
    max_val <- max(abs(mat), na.rm = TRUE)
    if (!is.finite(max_val) || max_val == 0) max_val <- 1
    col_fun <- circlize::colorRamp2(
      c(-max_val, 0, max_val), c("#3498db", "white", "#e74c3c"))
  } else {
    max_val <- max(mat, na.rm = TRUE)
    if (!is.finite(max_val) || max_val == 0) max_val <- 1
    col_fun <- circlize::colorRamp2(c(0, max_val), c("white", "#e74c3c"))
  }

  ComplexHeatmap::Heatmap(
    mat,
    name              = fill_label,
    col               = col_fun,
    cluster_rows      = TRUE,
    cluster_columns   = FALSE,
    show_column_names = TRUE,
    row_names_gp      = grid::gpar(fontsize = 8),
    column_names_gp   = grid::gpar(fontsize = 9),
    column_names_rot  = 45,
    row_names_max_width = grid::unit(12, "cm"),
    na_col            = "grey90",
    rect_gp           = grid::gpar(col = "white", lwd = 0.5),
    column_title      = paste0("ORA Heatmap – ", fill_label),
    column_title_gp   = grid::gpar(fontsize = 11, fontface = "bold")
  )
}

#' GSVA Pathway Heatmap
#'
#' Visualises a GSVA enrichment score matrix as a \code{ComplexHeatmap}.
#' Pass the matrix returned by \code{GSVA::gsva()} directly.
#' Rows = gene sets, columns = samples.
#'
#' @param gsva_mat Numeric matrix of GSVA enrichment scores, as returned by
#'   \code{GSVA::gsva()}.  Rows are gene sets, columns are samples.
#' @param condition Optional character vector of condition labels, one per
#'   sample column.  Used for the top colour-bar annotation.  If \code{NULL},
#'   no annotation is drawn.
#' @param top_n Integer. Number of most-variable gene sets to display.
#'   Default 25.
#' @param order_by_condition Logical. Sort columns by condition label
#'   (\code{TRUE}) or leave in the original order / cluster
#'   (\code{FALSE})?  Default \code{TRUE}.
#' @return A \code{ComplexHeatmap} object.
#' @examples
#' \dontrun{
#' param    <- GSVA::gsvaParam(exprData = expr_mat, geneSets = gene_sets)
#' gsva_mat <- GSVA::gsva(param, verbose = FALSE)
#' plot_gsva_heatmap(gsva_mat, condition = colData(se)$condition)
#' }
#' @seealso \code{\link{plot_ora_heatmap}}
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom circlize colorRamp2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom grid gpar unit
#' @export
plot_gsva_heatmap <- function(gsva_mat, condition = NULL,
                              top_n = 25L, order_by_condition = TRUE) {
  if (!is.matrix(gsva_mat) || !is.numeric(gsva_mat))
    stop("gsva_mat must be a numeric matrix (output of GSVA::gsva()).")

  row_vars <- apply(gsva_mat, 1, stats::var, na.rm = TRUE)
  top_idx  <- order(-row_vars)[seq_len(min(top_n, nrow(gsva_mat)))]
  gsva_mat <- gsva_mat[top_idx, , drop = FALSE]

  rownames(gsva_mat) <- gsub(
    "^HALLMARK_|^KEGG_|^REACTOME_|^WP_|^GOBP_|^GOMF_|^GOCC_",
    "", rownames(gsva_mat))
  rownames(gsva_mat) <- gsub("_", " ", rownames(gsva_mat))

  if (is.null(condition)) condition <- rep("Sample", ncol(gsva_mat))

  if (order_by_condition) {
    col_order <- order(condition)
    gsva_mat  <- gsva_mat[, col_order, drop = FALSE]
    condition <- condition[col_order]
  }

  cond_levels <- unique(condition)
  n_cond      <- length(cond_levels)
  palette     <- if (n_cond <= 8) {
    RColorBrewer::brewer.pal(max(3L, n_cond), "Set2")[seq_len(n_cond)]
  } else {
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_cond)
  }
  cond_colors <- stats::setNames(palette, cond_levels)

  top_anno <- ComplexHeatmap::HeatmapAnnotation(
    Condition = condition,
    col       = list(Condition = cond_colors),
    annotation_name_side = "left",
    annotation_name_gp   = grid::gpar(fontsize = 9)
  )

  max_val <- max(abs(gsva_mat), na.rm = TRUE)
  if (max_val == 0 || is.na(max_val)) max_val <- 1
  col_fun <- circlize::colorRamp2(
    c(-max_val, 0, max_val), c("#3498db", "white", "#e74c3c"))

  ComplexHeatmap::Heatmap(
    gsva_mat,
    name              = "GSVA score",
    col               = col_fun,
    top_annotation    = top_anno,
    cluster_rows      = TRUE,
    cluster_columns   = TRUE,
    show_column_names = TRUE,
    row_names_gp      = grid::gpar(fontsize = 8),
    column_names_gp   = grid::gpar(fontsize = 9),
    column_names_rot  = 45,
    row_names_max_width = grid::unit(12, "cm"),
    na_col            = "grey90",
    rect_gp           = grid::gpar(col = "white", lwd = 0.5),
    column_title_gp   = grid::gpar(fontsize = 11, fontface = "bold")
  )
}
