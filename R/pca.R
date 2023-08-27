#' @export
plot_pca <- function(dep, x = 1, y = 2, indicate = c("condition", "replicate"),
                     label = FALSE, n = 500, point_size = 8, label_size = 3, plot = TRUE, ID_col = "ID", exp = "LFQ", scale=F, interactive = F) {
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
    p <- ggplot(pca_df, aes(get(paste0("PC", x)), get(paste0("PC", y)))) +
      labs(
        title = paste0("PCA plot - top ", n, " variable features"),
        x = paste0("PC", x, ": ", percent[x], "%"),
        y = paste0("PC", y, ": ", percent[y], "%")
      ) +
      coord_fixed() + theme_bw()

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
      p <- p + geom_point(
        aes(
          col = pca_df[[indicate[1]]],
          shape = pca_df[[indicate[2]]]
        ),
        size = point_size
      ) +
        labs(
          col = indicate[1],
          shape = indicate[2]
        )
    }
    if (length(indicate) == 3) {
      p <- p + geom_point(
        aes(
          col = pca_df[[indicate[1]]],
          shape = pca_df[[indicate[2]]]
        ),
        size = point_size
      ) +
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
