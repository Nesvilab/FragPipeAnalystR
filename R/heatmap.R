#' Plot missing value heatmap
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

#' Plot correlation heatmap
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


#' Plot heatmap using differentially expressed features
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
