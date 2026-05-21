#' Principal Component Analysis (PCA) plot
#'
#' Creates a PCA visualization of proteomics data, with options for static
#' (ggplot2) or interactive (plotly) output.
#'
#' @param dep A \code{SummarizedExperiment} object containing proteomics data.
#' @param x Numeric value specifying which PC to plot on x-axis. Default is 1.
#' @param y Numeric value specifying which PC to plot on y-axis. Default is 2.
#' @param indicate Character vector specifying colData columns to use for
#'   coloring/shaping points. Maximum of 3 features allowed.
#'   Default is \code{c("condition", "replicate")}.
#' @param label Logical indicating whether to add sample labels to points.
#'   Default is \code{FALSE}.
#' @param n Numeric value specifying the number of top variable features to use
#'   for PCA calculation. Set to 0 to use all features. Default is 500.
#' @param point_size Numeric value specifying the size of points. Default is 8.
#' @param label_size Numeric value specifying the size of labels when
#'   \code{label = TRUE}. Default is 3.
#' @param plot Logical indicating whether to return a plot (\code{TRUE}) or
#'   a data frame of PCA coordinates (\code{FALSE}). Default is \code{TRUE}.
#' @param ID_col Character string specifying the column in colData to use for
#'   sample identification. Default is "sample_name".
#' @param exp Character string specifying the experiment type. If \code{NULL},
#'   uses the value from metadata. Default is \code{NULL}.
#' @param scale Logical indicating whether to scale features before PCA.
#'   Default is \code{FALSE}.
#' @param interactive Logical indicating whether to create an interactive
#'   plotly visualization (\code{TRUE}) or static ggplot2 plot (\code{FALSE}).
#'   Default is \code{FALSE}.
#'
#' @return If \code{plot = TRUE}, returns either a \code{ggplot} object (when
#'   \code{interactive = FALSE}) or a \code{plotly} object (when
#'   \code{interactive = TRUE}).
#'   If \code{plot = FALSE}, returns a data frame with PCA coordinates and
#'   indicated features.
#'
#' @details
#' The function calculates PCA using the top n most variable features.
#' For static plots:
#' \itemize{
#'   \item 1 indicate feature: mapped to color
#'   \item 2 indicate features: mapped to color and shape
#'   \item 3 indicate features: mapped to color, shape, and facets
#' }
#'
#' For interactive TMT plots with 2 indicate features, a dropdown menu allows
#' switching between condition-colored and plex-colored views.
#'
#' @examples
#' \dontrun{
#' # Basic PCA plot
#' plot_pca(se, indicate = "condition")
#'
#' # Interactive PCA with top 1000 features
#' plot_pca(se, n = 1000, interactive = TRUE)
#'
#' # PCA coordinates as data frame
#' pca_coords <- plot_pca(se, plot = FALSE)
#'
#' # Plot PC1 vs PC3
#' plot_pca(se, x = 1, y = 3)
#' }
#'
#' @seealso \code{\link{plot_correlation_heatmap}}
#'
#' @importFrom ggplot2 ggplot aes geom_point labs theme_bw theme element_blank
#'   element_line facet_wrap geom_text scale_shape_manual
#' @importFrom plotly plot_ly add_trace layout
#' @importFrom SummarizedExperiment assay colData
#' @importFrom S4Vectors metadata
#' @importFrom dplyr left_join select
#' @importFrom tibble rownames_to_column
#' @importFrom stats prcomp complete.cases
#'
#' @export
plot_pca <- function(dep, x = 1, y = 2, indicate = c("condition", "replicate"),
                     label = FALSE, n = 500, point_size = 8, label_size = 3, plot = TRUE, ID_col = "sample_name", exp = NULL, scale=F, interactive = F) {
  if (is.null(exp)) {
    exp <- metadata(dep)$exp
  }

  if (is.integer(x)) x <- as.numeric(x)
  if (is.integer(y)) y <- as.numeric(y)
  if (is.integer(n)) n <- as.numeric(n)
  if (is.integer(point_size)) point_size <- as.numeric(point_size)
  if (is.integer(label_size)) label_size <- as.numeric(label_size)
  # Show error if inputs are not the required classes
  assertthat::assert_that(
    inherits(dep, "SummarizedExperiment"),
    is.numeric(x),
    length(x) == 1,
    is.numeric(y),
    length(y) == 1,
    is.numeric(n),
    length(n) == 1,
    is.character(indicate),
    is.logical(label),
    is.numeric(point_size),
    length(point_size) == 1,
    is.numeric(label_size),
    length(label_size) == 1,
    is.logical(plot),
    length(plot) == 1
  )

  # Check for valid x and y values
  if (x > ncol(dep) | y > ncol(dep)) {
    stop(
      paste0(
        "'x' and/or 'y' arguments are not valid\n",
        "Run plot_pca() with 'x' and 'y' <= ",
        ncol(dep), "."
      ),
      call. = FALSE
    )
  }


  # Check for valid 'indicate'
  columns <- colnames(colData(dep))
  if (!is.null(indicate)) {
    if (length(indicate) > 3) {
      stop("Too many features in 'indicate'
        Run plot_pca() with a maximum of 3 indicate features")
    }
    if (any(!indicate %in% columns)) {
      stop(
        paste0(
          "'",
          paste0(indicate, collapse = "' and/or '"),
          "' column(s) is/are not present in ",
          deparse(substitute(dep)),
          ".\nValid columns are: '",
          paste(columns, collapse = "', '"),
          "'."
        ),
        call. = FALSE
      )
    }
  }

  data <- assay(dep)
  data <- data[complete.cases(data), ]
  # Get the variance per protein and take the top n variable proteins
  var <- apply(data, 1, sd)
  # Check for valid 'n' value
  if (n == 0) {
    df <- data
    n <- nrow(data)
  } else if (n > nrow(data)) {
    message(paste(
      "'n' argument is larger than number of features availble(",
      nrow(data), ").", nrow(data), "features will be used for PCA calculation."
    ))
    df <- data
    n <- nrow(data)
  } else {
    df <- data[order(var, decreasing = TRUE)[seq_len(n)], ]
  }

  # Calculate PCA
  pca <- prcomp(t(df), scale = scale)
  pca_df <- pca$x %>%
    data.frame() %>%
    rownames_to_column() %>%
    left_join(., data.frame(colData(dep)), by = c("rowname" = ID_col))

  # Calculate the percentage of variance explained
  percent <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  # Make factors of indicate features
  for (feat in indicate) {
    pca_df[[feat]] <- as.factor(pca_df[[feat]])
  }

  if (interactive) {
    if (length(indicate) == 1) {
      p <- plot_ly(data = pca_df, type = "scatter", mode = "markers", marker = list(size = point_size)) %>%
        plotly::layout(
            title = paste0("PCA plot (", n, " features used)"),
            xaxis = list(title = paste0("PC", x, ": ", percent[x], "%")),
            yaxis = list(title = paste0("PC", y, ": ", percent[y], "%"))
            ) %>%
        add_trace(
          type = "scatter",
          x = ~PC1,
          y = ~PC2,
          color = as.formula(paste0("~", indicate[1])),
          mode = "markers",
          legendgroup = indicate[1],
          legendgrouptitle_text = indicate[1]
          )
      } else if (length(indicate) == 2) {
      if (exp == "TMT") {
        pca_df$plex <- as.factor(pca_df$plex)
        p <- plot_ly() %>%
          # Overlay color for gears
          add_trace(
            data = pca_df, type = "scatter",
            x = ~PC1,
            y = ~PC2,
            symbol = as.formula(paste0("~", indicate[2])),
            marker = list(color = "grey", size = point_size + 3),
            text = pca_df$sample_name,
            hoverinfo = "text",
            mode = "markers",
            legendgroup = indicate[2],
            legendgrouptitle_text = indicate[2]
          ) %>%
          add_trace(
            data = pca_df, type = "scatter",
            x = ~PC1,
            y = ~PC2,
            marker = list(size = point_size),
            color = as.formula(paste0("~", indicate[1])),
            mode = "markers",
            text = pca_df$sample_name,
            hoverinfo = "text",
            legendgroup = indicate[1],
            legendgrouptitle_text = indicate[1]
          ) %>%
          add_trace(
            data = pca_df, type = "scatter",
            x = ~PC1,
            y = ~PC2,
            marker = list(size = point_size),
            color = as.formula(paste0("~", "plex")),
            text = pca_df$sample_name,
            hoverinfo = "text",
            mode = "markers",
            legendgroup = "plex",
            legendgrouptitle_text = "plex",
            xaxis = "x2", yaxis = "y2", visible = F
          ) %>%
          plotly::layout(
            title = paste0("PCA plot (", n, " features used)"),
            xaxis = list(title = paste0("PC", x, ": ", percent[x], "%")),
            xaxis2 = list(title = paste0("PC", x, ": ", percent[x], "%"), overlaying = "x", visible = F),
            yaxis = list(title = paste0("PC", y, ": ", percent[y], "%")),
            yaxis2 = list(title = paste0("PC", y, ": ", percent[y], "%"), overlaying = "y", visible = F),
            legend = list(
              itemclick = FALSE,
              itemdoubleclick = FALSE,
              groupclick = FALSE
            ),
            updatemenus = list(
              list(
                y = 0.8,
                buttons = list(
                  list(
                    method = "update",
                    args = list(
                      list(visible = unlist(Map(rep, x = c(T, T, F), each = c(
                        length(unique(pca_df$condition)),
                        length(unique(pca_df$replicate)),
                        length(unique(pca_df$plex))
                      )))),
                      list(
                        xaxis = list(
                          title = paste0("PC", x, ": ", percent[x], "%"),
                          visible = TRUE
                        ),
                        xaxis2 = list(overlaying = "x", visible = FALSE),
                        yaxis = list(
                          title = paste0("PC", y, ": ", percent[y], "%"),
                          visible = TRUE
                        ),
                        yaxis2 = list(overlaying = "y", visible = FALSE)
                      )
                    ),
                    label = "by condition"
                  ),
                  list(
                    method = "update",
                    args = list(
                      list(visible = unlist(Map(rep, x = c(F, F, T), each = c(
                        length(unique(pca_df$condition)),
                        length(unique(pca_df$replicate)),
                        length(unique(pca_df$plex))
                      )))),
                      list(
                        xaxis = list(visible = F),
                        xaxis2 = list(
                          title = paste0("PC", x, ": ", percent[x], "%"),
                          overlaying = "x", visible = T
                        ),
                        yaxis = list(visible = F),
                        yaxis2 = list(
                          title = paste0("PC", y, ": ", percent[y], "%"),
                          overlaying = "y", visible = T
                        )
                      )
                    ),
                    label = "by plex"
                  )
                )
              )
            )
          )
      } else {
        p <- plot_ly(
          data = pca_df, type = "scatter",
          mode = "markers", marker = list(size = point_size), text = ~rowname
        ) %>%
          # Overlay color for gears
          add_trace(
            type = "scatter",
            x = ~PC1,
            y = ~PC2,
            symbol = as.formula(paste0("~", indicate[2])),
            marker = list(color = "grey", size = point_size + 3),
            mode = "markers",
            text = pca_df$sample_name,
            hoverinfo = "text",
            legendgroup = indicate[2],
            legendgrouptitle_text = indicate[2]
          ) %>%
          add_trace(
            type = "scatter",
            x = ~PC1,
            y = ~PC2,
            color = as.formula(paste0("~", indicate[1])),
            mode = "markers",
            text = pca_df$sample_name,
            hoverinfo = "text",
            legendgroup = indicate[1],
            legendgrouptitle_text = indicate[1]
          ) %>%
          plotly::layout(
            title = paste0("PCA plot (", n, " features used)"),
            xaxis = list(title = paste0("PC", x, ": ", percent[x], "%")),
            yaxis = list(title = paste0("PC", y, ": ", percent[y], "%")),
            legend = list(
              itemclick = FALSE,
              itemdoubleclick = FALSE,
              groupclick = FALSE
            )
          )
      }
    }
  } else { # static plot by ggplot2
    shape_palette <- c(16, 17, 15, 3, 7, 8, 0, 1, 2, 4, 5, 6, 9, 10, 11, 12, 13, 14, 18, 20, 21, 22, 23, 24, 25, 19)

    p <- ggplot(pca_df, aes(get(paste0("PC", x)), get(paste0("PC", y)))) +
      labs(
        title = paste0("PCA plot - top ", n, " variable features"),
        x = paste0("PC", x, ": ", percent[x], "%"),
        y = paste0("PC", y, ": ", percent[y], "%")
      ) +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

    if (length(indicate) == 0) {
      p <- p + geom_point(size = point_size)
    }
    if (length(indicate) == 1) {
      p <- p + geom_point(aes(col = pca_df[[indicate[1]]]),
        size = point_size
      ) +
        labs(col = indicate[1])
    }
    if (length(indicate) == 2) {
      n_shapes <- nlevels(pca_df[[indicate[2]]])
      p <- p + geom_point(
        aes(
          col = pca_df[[indicate[1]]],
          shape = pca_df[[indicate[2]]]
        ),
        size = point_size
      ) +
        scale_shape_manual(values = shape_palette[seq_len(n_shapes)]) +
        labs(
          col = indicate[1],
          shape = indicate[2]
        )
    }
    if (length(indicate) == 3) {
      n_shapes <- nlevels(pca_df[[indicate[2]]])
      p <- p + geom_point(
        aes(
          col = pca_df[[indicate[1]]],
          shape = pca_df[[indicate[2]]]
        ),
        size = point_size
      ) +
        scale_shape_manual(values = shape_palette[seq_len(n_shapes)]) +
        facet_wrap(~ pca_df[[indicate[3]]])
      labs(
        col = indicate[1],
        shape = indicate[2]
      )
    }
    if (label) {
      p <- p + geom_text(aes(label = rowname), size = label_size)
    }
  }

  if (plot) {
    return(p)
  } else {
    df <- pca_df %>%
      select(rowname, paste0("PC", c(x, y)), match(indicate, colnames(pca_df)))
    colnames(df)[1] <- "sample"
    return(df)
  }
}
