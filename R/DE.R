#' @export
test_limma <- function(se, type = c("control", "all", "manual"),
                       control = NULL, test = NULL,
                       design_formula = formula(~ 0 + condition),
                       paired = FALSE) {
  # require("dplyr", "tidyr", "purrr")

  # Show error if inputs are not the required classes
  assertthat::assert_that(
    inherits(se, "SummarizedExperiment"),
    is.character(type),
    class(design_formula) == "formula"
  )
  if (paired == FALSE) {
    design_formula <- design_formula
  } else {
    design_formula <- formula(~ 0 + condition + replicate)
  }


  # Show error if inputs do not contain required columns
  type <- match.arg(type)

  col_data <- colData(se)
  raw <- assay(se)

  if (any(!c("name") %in% colnames(rowData(se)))) {
    stop("'name' column is not present in '",
      deparse(substitute(se)),
      "'\nRun make_unique() and make_se() to obtain the required columns",
      call. = FALSE
    )
  }
  if (any(!c("label", "condition", "replicate") %in% colnames(col_data))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
      deparse(substitute(se)),
      "'\nRun make_se() or make_se_parse() to obtain the required columns",
      call. = FALSE
    )
  }
  if (any(is.na(raw))) {
    warning("Missing values in '", deparse(substitute(se)), "'")
  }

  if (!is.null(control)) {
    # Show error if control input is not valid
    assertthat::assert_that(
      is.character(control),
      length(control) == 1
    )
    if (!control %in% unique(col_data$condition)) {
      stop("run test_diff() with a valid control.\nValid controls are: '",
        paste0(unique(col_data$condition), collapse = "', '"), "'",
        call. = FALSE
      )
    }
  }

  # variables in formula
  variables <- terms.formula(design_formula) %>%
    attr(., "variables") %>%
    as.character() %>%
    .[-1]

  # Throw error if variables are not col_data columns
  if (any(!variables %in% colnames(col_data))) {
    stop("run make_diff() with an appropriate 'design_formula'")
  }
  if (variables[1] != "condition") {
    stop("first factor of 'design_formula' should be 'condition'")
  }

  # Obtain variable factors
  for (var in variables) {
    temp <- factor(col_data[[var]])
    assign(var, temp)
  }

  # Make an appropriate design matrix
  design <- model.matrix(design_formula, data = environment())
  colnames(design) <- gsub("condition", "", colnames(design))

  # Generate contrasts to be tested
  # Either make all possible combinations ("all"),
  # only the contrasts versus the control sample ("control") or
  # use manual contrasts
  conditions <- as.character(unique(condition))
  if (type == "all") {
    # All possible combinations
    cntrst <- apply(utils::combn(conditions, 2), 2, paste, collapse = " - ")

    if (!is.null(control)) {
      # Make sure that contrast containing
      # the control sample have the control as denominator
      flip <- grep(paste("^", control, sep = ""), cntrst)
      if (length(flip) >= 1) {
        cntrst[flip] <- cntrst[flip] %>%
          gsub(paste(control, "- ", sep = " "), "", .) %>%
          paste(" - ", control, sep = "")
      }
    }
  }
  if (type == "control") {
    # Throw error if no control argument is present
    if (is.null(control)) {
      stop("run test_diff(type = 'control') with a 'control' argument")
    }

    # Make contrasts
    cntrst <- paste(conditions[!conditions %in% control],
      control,
      sep = " - "
    )
  }
  if (type == "manual") {
    # Throw error if no test argument is present
    if (is.null(test)) {
      stop("run test_diff(type = 'manual') with a 'test' argument")
    }
    assertthat::assert_that(is.character(test))

    if (any(!unlist(strsplit(test, "_vs_")) %in% conditions)) {
      stop("run test_diff() with valid contrasts in 'test'",
        ".\nValid contrasts should contain combinations of: '",
        paste0(conditions, collapse = "', '"),
        "', for example '", paste0(conditions[1], "_vs_", conditions[2]),
        "'.",
        call. = FALSE
      )
    }

    cntrst <- gsub("_vs_", " - ", test)
  }
  # Print tested contrasts
  message(
    "Tested contrasts: ",
    paste(gsub(" - ", "_vs_", cntrst), collapse = ", ")
  )

  # Test for differential expression by empirical Bayes moderation
  # of a linear model on the predefined contrasts
  fit <- lmFit(raw, design = design)
  made_contrasts <- makeContrasts(contrasts = cntrst, levels = design)
  contrast_fit <- contrasts.fit(fit, made_contrasts)

  if (any(is.na(raw))) {
    for (i in cntrst) {
      covariates <- strsplit(i, " - ") %>% unlist()
      single_contrast <- makeContrasts(contrasts = i, levels = design[, covariates])
      single_contrast_fit <- contrasts.fit(fit[, covariates], single_contrast)
      contrast_fit$coefficients[, i] <- single_contrast_fit$coefficients[, 1]
      contrast_fit$stdev.unscaled[, i] <- single_contrast_fit$stdev.unscaled[, 1]
    }
  }

  eB_fit <- eBayes(contrast_fit)

  # function to retrieve the results of
  # the differential expression test using 'fdrtool'
  retrieve_fun <- function(comp, fit = eB_fit) {
    res <- topTable(fit,
      sort.by = "t", adjust.method = "BH", coef = comp,
      number = Inf, confint = TRUE
    )
    # res <- res[!is.na(res$t),]
    # fdr_res <- fdrtool::fdrtool(res$t, plot = FALSE, verbose = FALSE)
    # res$qval <- res$adj.P.Value
    # res$lfdr <- fdr_res$lfdr
    res$comparison <- rep(comp, dim(res)[1])
    res <- tibble::rownames_to_column(res)
    return(res)
  }

  # limma_res<- topTable(eB_fit, sort.by = 'B', adjust.method="BH", coef = cntrst, number = Inf, confint = T )
  # limma_res$comparison <- rep(cntrst, dim(limma_res)[1])
  # limma_res <- rownames_to_column(limma_res)
  # Retrieve the differential expression test results
  limma_res <- purrr::map_df(cntrst, retrieve_fun)

  # Select the logFC, CI and qval variables
  table <- limma_res %>%
    dplyr::select(rowname, logFC, CI.L, CI.R, P.Value, adj.P.Val, comparison) %>%
    dplyr::mutate(comparison = gsub(" - ", "_vs_", comparison)) %>%
    tidyr::gather(variable, value, -c(rowname, comparison)) %>%
    dplyr::mutate(variable = dplyr::recode(variable, logFC = "diff", P.Value = "p.val", adj.P.Val = "p.adj")) %>%
    tidyr::unite(temp, comparison, variable) %>%
    tidyr::spread(temp, value)
  rowData(se) <- as.data.frame(left_join(as.data.frame(rowData(se)), table,
    by = c("name" = "rowname")
  ))
  return(se)
}

#' @export
test_diff <- function(se, type = c("control", "all", "manual"),
                      control = NULL, test = NULL,
                      design_formula = formula(~ 0 + condition)) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(
    inherits(se, "SummarizedExperiment"),
    is.character(type),
    class(design_formula) == "formula"
  )

  # Show error if inputs do not contain required columns
  type <- match.arg(type)

  col_data <- colData(se)
  raw <- assay(se)

  if (any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '",
      deparse(substitute(se)),
      "'\nRun make_unique() and make_se() to obtain the required columns",
      call. = FALSE
    )
  }
  if (any(!c("label", "condition", "replicate") %in% colnames(col_data))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
      deparse(substitute(se)),
      "'\nRun make_se() or make_se_parse() to obtain the required columns",
      call. = FALSE
    )
  }
  if (any(is.na(raw))) {
    warning("Missing values in '", deparse(substitute(se)), "'")
  }

  if (!is.null(control)) {
    # Show error if control input is not valid
    assertthat::assert_that(
      is.character(control),
      length(control) == 1
    )
    if (!control %in% unique(col_data$condition)) {
      stop("run test_diff() with a valid control.\nValid controls are: '",
        paste0(unique(col_data$condition), collapse = "', '"), "'",
        call. = FALSE
      )
    }
  }

  # variables in formula
  variables <- terms.formula(design_formula) %>%
    attr(., "variables") %>%
    as.character() %>%
    .[-1]

  # Throw error if variables are not col_data columns
  if (any(!variables %in% colnames(col_data))) {
    stop("run make_diff() with an appropriate 'design_formula'")
  }
  # if(variables[1] != "condition") {
  #   stop("first factor of 'design_formula' should be 'condition'")
  # }

  # Obtain variable factors
  for (var in variables) {
    temp <- factor(col_data[[var]])
    assign(var, temp)
  }

  # Make an appropriate design matrix
  design <- model.matrix(design_formula, data = environment())
  colnames(design) <- gsub("condition", "", colnames(design))

  # Generate contrasts to be tested
  # Either make all possible combinations ("all"),
  # only the contrasts versus the control sample ("control") or
  # use manual contrasts
  conditions <- as.character(unique(col_data$condition))
  if (type == "all") {
    # All possible combinations
    cntrst <- apply(utils::combn(conditions, 2), 2, paste, collapse = " - ")

    if (!is.null(control)) {
      # Make sure that contrast containing
      # the control sample have the control as denominator
      flip <- grep(paste("^", control, sep = ""), cntrst)
      if (length(flip) >= 1) {
        cntrst[flip] <- cntrst[flip] %>%
          gsub(paste(control, "- ", sep = " "), "", .) %>%
          paste(" - ", control, sep = "")
      }
    }
  }
  if (type == "control") {
    # Throw error if no control argument is present
    if (is.null(control)) {
      stop("run test_diff(type = 'control') with a 'control' argument")
    }

    # Make contrasts
    cntrst <- paste(conditions[!conditions %in% control],
      control,
      sep = " - "
    )
  }
  if (type == "manual") {
    # Throw error if no test argument is present
    if (is.null(test)) {
      stop("run test_diff(type = 'manual') with a 'test' argument")
    }
    assertthat::assert_that(is.character(test))

    # if(any(!unlist(strsplit(test, "_vs_")) %in% conditions)) {
    #   stop("run test_diff() with valid contrasts in 'test'",
    #        ".\nValid contrasts should contain combinations of: '",
    #        paste0(conditions, collapse = "', '"),
    #        "', for example '", paste0(conditions[1], "_vs_", conditions[2]),
    #        "'.", call. = FALSE)
    # }

    cntrst <- gsub("_vs_", " - ", test)
  }
  # Print tested contrasts
  message(
    "Tested contrasts: ",
    paste(gsub(" - ", "_vs_", cntrst), collapse = ", ")
  )

  # Test for differential expression by empirical Bayes moderation
  # of a linear model on the predefined contrasts
  # print(design)
  # print(cntrst)
  fit <- lmFit(raw, design = design)
  made_contrasts <- makeContrasts(contrasts = cntrst, levels = design)
  # print(made_contrasts)
  contrast_fit <- contrasts.fit(fit, made_contrasts)

  if (type != "manual") {
    if (any(is.na(raw))) {
      for (i in cntrst) {
        covariates <- strsplit(i, " - ") %>% unlist()
        single_contrast <- makeContrasts(contrasts = i, levels = design[, covariates])
        single_contrast_fit <- contrasts.fit(fit[, covariates], single_contrast)
        contrast_fit$coefficients[, i] <- single_contrast_fit$coefficients[, 1]
        contrast_fit$stdev.unscaled[, i] <- single_contrast_fit$stdev.unscaled[, 1]
      }
    }
  }
  eB_fit <- eBayes(contrast_fit)

  # function to retrieve the results of
  # the differential expression test using 'fdrtool'
  retrieve_fun <- function(comp, fit = eB_fit) {
    res <- topTable(fit,
      sort.by = "t", coef = comp,
      number = Inf, confint = TRUE
    )
    res <- res[!is.na(res$t), ]
    fdr_res <- fdrtool::fdrtool(res$t, plot = FALSE, verbose = FALSE)
    res$qval <- fdr_res$qval
    res$lfdr <- fdr_res$lfdr
    res$comparison <- rep(comp, dim(res)[1])
    res <- rownames_to_column(res)
    return(res)
  }

  # Retrieve the differential expression test results
  limma_res <- map_df(cntrst, retrieve_fun)

  # Select the logFC, CI and qval variables
  table <- limma_res %>%
    select(rowname, logFC, CI.L, CI.R, P.Value, qval, comparison) %>%
    mutate(comparison = gsub(" - ", "_vs_", comparison)) %>%
    gather(variable, value, -c(rowname, comparison)) %>%
    mutate(variable = recode(variable, logFC = "diff", P.Value = "p.val", qval = "p.adj")) %>%
    unite(temp, comparison, variable) %>%
    spread(temp, value)
  rowData(se) <- as.data.frame(left_join(as.data.frame(rowData(se)), table,
    by = c("name" = "rowname")
  ))
  return(se)
}

# customized from DEP's add_rejections: https://github.com/arnesmits/DEP/blob/b425d8d0db67b15df4b8bcf87729ef0bf5800256/R/functions.R#L1032
#' Mark significant proteins
#'
#' \code{add_rejections} marks significant proteins based on defined cutoffs.
#'
#' @param diff SummarizedExperiment,
#' Proteomics dataset on which differential enrichment analysis
#' has been performed (output from \code{\link{test_diff}()}).
#' @param alpha Numeric(1),
#' Sets the threshold for the adjusted P value.
#' @param lfc Numeric(1),
#' Sets the threshold for the log2 fold change.
#' @return A SummarizedExperiment object
#' annotated with logical columns indicating significant proteins.
#' @export
add_rejections <- function(diff, alpha = 0.05, lfc = 1) {
  # Show error if inputs are not the required classes
  if (is.integer(alpha)) alpha <- as.numeric(alpha)
  if (is.integer(lfc)) lfc <- as.numeric(lfc)
  assertthat::assert_that(
    inherits(diff, "SummarizedExperiment"),
    is.numeric(alpha),
    length(alpha) == 1,
    is.numeric(lfc),
    length(lfc) == 1
  )

  row_data <- rowData(diff, use.names = FALSE) %>%
    as.data.frame()
  # Show error if inputs do not contain required columns
  if (any(!c("name", "ID") %in% colnames(row_data))) {
    stop("'name' and/or 'ID' columns are not present in '",
      deparse(substitute(diff)),
      "'\nRun make_unique() and make_se() to obtain the required columns",
      call. = FALSE
    )
  }
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in '",
      deparse(substitute(diff)),
      "'\nRun test_diff() to obtain the required columns",
      call. = FALSE
    )
  }

  # get all columns with adjusted p-values and log2 fold changes
  cols_p <- grep("_p.adj", colnames(row_data))
  cols_diff <- grep("_diff", colnames(row_data))

  # Mark differential expressed proteins by
  # applying alpha and log2FC parameters per protein
  if (length(cols_p) == 1) {
    rowData(diff)$significant <-
      row_data[, cols_p] <= alpha & abs(row_data[, cols_diff]) >= lfc
    rowData(diff)$contrast_significant <-
      rowData(diff, use.names = FALSE)$significant
    colnames(rowData(diff))[ncol(rowData(diff, use.names = FALSE))] <-
      gsub("p.adj", "significant", colnames(row_data)[cols_p])
  }
  if (length(cols_p) > 1) {
    p_reject <- row_data[, cols_p] <= alpha
    p_reject[is.na(p_reject)] <- FALSE
    diff_reject <- abs(row_data[, cols_diff]) >= lfc
    diff_reject[is.na(diff_reject)] <- FALSE
    sign_df <- p_reject & diff_reject
    sign_df <- cbind(sign_df,
      significant = apply(sign_df, 1, function(x) any(x))
    )
    colnames(sign_df) <- gsub("_p.adj", "_significant", colnames(sign_df))
    sign_df <- cbind(name = row_data$name, as.data.frame(sign_df))
    rowData(diff) <- as.data.frame(left_join(as.data.frame(rowData(diff)), sign_df,
      by = c("name" = "name")
    ))
  }
  return(diff)
}

# plot_volcano_new from FragPipe-Analyst
#' @export
plot_volcano <- function(dep, contrast, label_size = 3,
                         add_names = TRUE, adjusted = T, plot = TRUE, alpha = 0.05, lfc = 1) {
  # Show error if inputs are not the required classes
  if (is.integer(label_size)) label_size <- as.numeric(label_size)
  assertthat::assert_that(
    inherits(dep, "SummarizedExperiment"),
    is.character(contrast),
    length(contrast) == 1,
    is.numeric(label_size),
    length(label_size) == 1,
    is.logical(add_names),
    length(add_names) == 1,
    is.logical(adjusted),
    length(adjusted) == 1,
    is.logical(plot),
    length(plot) == 1
  )

  row_data <- rowData(dep, use.names = FALSE)

  # Show error if inputs do not contain required columns
  if (any(!c("name", "ID") %in% colnames(row_data))) {
    stop(
      paste0(
        "'name' and/or 'ID' columns are not present in '",
        deparse(substitute(dep)),
        "'.\nRun make_unique() to obtain required columns."
      ),
      call. = FALSE
    )
  }
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop(
      paste0(
        "'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '",
        deparse(substitute(dep)),
        "'.\nRun test_diff() to obtain the required columns."
      ),
      call. = FALSE
    )
  }
  if (length(grep("_significant", colnames(row_data))) < 1) {
    stop(
      paste0(
        "'[contrast]_significant' columns are not present in '",
        deparse(substitute(dep)),
        "'.\nRun add_rejections() to obtain the required columns."
      ),
      call. = FALSE
    )
  }

  # Show error if an unvalid contrast is given
  if (length(grep(
    paste("^", contrast, "_diff", sep = ""),
    colnames(row_data)
  )) == 0) {
    valid_cntrsts <- row_data %>%
      data.frame() %>%
      select(ends_with("_diff")) %>%
      colnames(.) %>%
      gsub("_diff", "", .)
    valid_cntrsts_msg <- paste0(
      "Valid contrasts are: '",
      paste0(valid_cntrsts, collapse = "', '"),
      "'"
    )
    stop("Not a valid contrast, please run `plot_volcano()` with a valid contrast as argument\n",
      valid_cntrsts_msg,
      call. = FALSE
    )
  }

  # Generate a data.frame containing all info for the volcano plot
  diff <- grep(
    paste("^", contrast, "_diff", sep = ""),
    colnames(row_data)
  )
  if (adjusted) {
    p_values <- grep(
      paste("^", contrast, "_p.adj", sep = ""),
      colnames(row_data)
    )
  } else {
    p_values <- grep(
      paste("^", contrast, "_p.val", sep = ""),
      colnames(row_data)
    )
  }
  signif <- abs(row_data[, diff]) >= lfc & row_data[, p_values] <= alpha

  df_tmp <- data.frame(
    diff = row_data[, diff],
    p_values = -log10(row_data[, p_values]),
    signif = signif,
    name = row_data$name
  )
  df <- df_tmp %>%
    data.frame() %>%
    filter(!is.na(signif)) %>%
    arrange(signif)
  name1 <- gsub("_vs_.*", "", contrast)
  name2 <- gsub(".*_vs_", "", contrast)
  # return(df)
  # Plot volcano with or without labels

  p <- ggplot(df, aes(diff, p_values)) +
    geom_vline(xintercept = 0) +
    geom_point(aes(col = signif)) +
    geom_text(data = data.frame(), aes(
      x = c(Inf, -Inf),
      y = c(-Inf, -Inf),
      hjust = c(1, 0),
      vjust = c(-1, -1),
      label = c(name1, name2),
      size = 5,
      fontface = "bold"
    )) +
    labs(
      title = contrast,
      x = expression(log[2] ~ "Fold change")
    ) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey"))
  if (add_names) {
    p <- p + ggrepel::geom_text_repel(
      data = filter(df, signif),
      aes(label = name),
      size = label_size,
      box.padding = unit(0.1, "lines"),
      point.padding = unit(0.1, "lines"),
      segment.size = 0.5
    )
  }
  if (adjusted) {
    p <- p + labs(y = expression(-log[10] ~ "Adjusted p-value"))
  } else {
    p <- p + labs(y = expression(-log[10] ~ "P-value"))
  }
  if (plot) {
    # return(list(p, df))
    # return(df)
    return(p)
  } else {
    df <- df %>%
      select(name, diff, p_values, signif)
    colnames(df)[c(1, 2, 3)] <- c("protein", "log2_fold_change", "p_value_-log10")
    if (adjusted) {
      colnames(df)[3] <- "adjusted_p_value_-log10"
    }
    return(df)
  }
}

# From DEP
theme_DEP1 <- function() {
  # Use theme_bw() as default
  basesize <- 12
  theme <- ggplot2::theme_bw(base_size = basesize)

  # Change plot title appearance
  theme$plot.title$face <- "bold"
  theme$plot.title$size <- basesize + 2
  theme$plot.title$hjust <- 0.5

  # Change axis title appearance
  theme$axis.title.x$size <- basesize + 2

  theme$axis.title.y$size <- basesize + 2

  # Change axis text appearance
  theme$axis.text$size <- basesize
  theme$axis.text$colour <- "black"

  # Change legend title appearance
  theme$legend.title$size <- basesize + 2

  # Change legend text appearance
  theme$legend.text$size <- basesize

  # Change strip text (facet headers) appearance
  theme$strip.text$face <- "bold"
  theme$strip.text$size <- basesize + 2
  theme$strip.text$colour <- "black"

  return(theme)
}
