#' Differential expression analysis using limma
#'
#' Performs differential expression analysis on a SummarizedExperiment object
#' using the limma package with empirical Bayes moderation.
#'
#' @param se A \code{SummarizedExperiment} object containing log2-transformed
#'   proteomics data. Must contain 'name' column in rowData and 'label',
#'   'condition', 'replicate' columns in colData.
#' @param type Character string specifying the type of contrasts to test.
#'   Options are:
#'   \itemize{
#'     \item \code{"control"}: Test all conditions against a specified control.
#'     \item \code{"all"}: Test all pairwise combinations of conditions.
#'     \item \code{"others"}: Test each condition against all other conditions combined.
#'     \item \code{"manual"}: Use manually specified contrasts via the \code{test} parameter.
#'   }
#' @param control Character vector specifying the control condition(s).
#'   Required when \code{type = "control"}. Can specify multiple controls.
#' @param test Character vector of manual contrasts to test when
#'   \code{type = "manual"}. Contrasts should be in the format
#'   "condition1_vs_condition2".
#' @param design_formula A formula specifying the design matrix.
#'   Default is \code{~ 0 + condition}. The first term should be 'condition'.
#' @param paired Logical indicating whether to perform paired analysis.
#'   If \code{TRUE}, adds replicate to the design formula.
#' @param numeric_var Character vector of variable names in colData that should
#'   be treated as numeric covariates rather than factors.
#'
#' @return A \code{SummarizedExperiment} object with differential expression
#'   results added to rowData. For each contrast, the following columns are added:
#'   \itemize{
#'     \item \code{[contrast]_diff}: Log2 fold change
#'     \item \code{[contrast]_CI.L}: Lower confidence interval
#'     \item \code{[contrast]_CI.R}: Upper confidence interval
#'     \item \code{[contrast]_p.val}: Raw p-value
#'     \item \code{[contrast]_p.adj}: BH-adjusted p-value
#'   }
#'
#' @examples
#' \dontrun{
#' # Test all conditions against control
#' se_diff <- test_limma(se, type = "control", control = "Ctrl")
#'
#' # Test all pairwise comparisons
#' se_diff <- test_limma(se, type = "all")
#'
#' # Manual contrast
#' se_diff <- test_limma(se, type = "manual", test = "Treatment_vs_Control")
#'
#' # Paired analysis
#' se_diff <- test_limma(se, type = "control", control = "Ctrl", paired = TRUE)
#' }
#'
#' @seealso \code{\link{test_diff}}, \code{\link{add_rejections}}, \code{\link{plot_volcano}}
#'
#' @export
test_limma <- function(se, type = c("control", "all", "others", "manual"),
                       control = NULL, test = NULL,
                       design_formula = formula(~ 0 + condition),
                       paired = FALSE, numeric_var = NULL) {

  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.character(type),
                          class(design_formula) == "formula")
  if (paired == FALSE){
    design_formula <- design_formula
  } else {
    design_formula <- formula(~ 0 + condition + replicate)
  }


  # Show error if inputs do not contain required columns
  type <- match.arg(type)

  col_data <- colData(se)
  raw <- assay(se)

  if(any(!c("name") %in% colnames(rowData(se)))) {
    stop("'name' column is not present in '",
         deparse(substitute(se)),
         "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if(any(!c("label", "condition", "replicate") %in% colnames(col_data))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
         deparse(substitute(se)),
         "'\nRun make_se() or make_se_parse() to obtain the required columns",
         call. = FALSE)
  }
  if(any(is.na(raw))) {
    warning("Missing values in '", deparse(substitute(se)), "'")
  }

  if(!is.null(control)) {
    # Show error if control input is not valid
    if(any(!control %in% unique(col_data$condition))) {
      stop("run test_limma() with a valid control.\nValid controls are: '",
           paste0(unique(col_data$condition), collapse = "', '"), "'",
           call. = FALSE)
    }
  }

  # variables in formula
  variables <- terms.formula(design_formula) %>%
    attr(., "variables") %>%
    as.character() %>%
    .[-1]

  # Throw error if variables are not col_data columns
  if(any(!variables %in% colnames(col_data))) {
    stop("run make_diff() with an appropriate 'design_formula'")
  }
  if(variables[1] != "condition") {
    stop("first factor of 'design_formula' should be 'condition'")
  }

  # Obtain variable factors
  for(var in variables) {
    val <- col_data[[var]]
    if (var %in% numeric_var) {
      temp <- as.numeric(val)
    } else {
      temp <- factor(val)
    }
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
  if(type == "all") {
    # All possible combinations
    cntrst <- apply(utils::combn(conditions, 2), 2, paste, collapse = " - ")

    if(!is.null(control)) {
      # Make sure that contrast containing
      # the control sample have the control as denominator
      flip <- grep(paste("^", control, sep = ""), cntrst)
      if(length(flip) >= 1) {
        cntrst[flip] <- cntrst[flip] %>%
          gsub(paste(control, "- ", sep = " "), "", .) %>%
          paste(" - ", control, sep = "")
      }
    }

  } else if(type == "control") {
    if(is.null(control)) { # Throw error if no control argument is present
      stop("Please sepecify 'control' condition")
    } else { # Make contrasts
      conditions_other_than_control <- conditions[!conditions %in% control]
      cntrst <- c()
      for (i in 1:length(control)) {
        cntrst <- c(cntrst,
                    paste(conditions_other_than_control,
                          control[i],sep = " - "))
      }
    }
  } else if(type == "manual") {
    if(is.null(test)) { # Throw error if no test argument is present
      stop("run test_diff(type = 'manual') with a 'test' argument")
    }
    assertthat::assert_that(is.character(test))

    if(any(!unlist(strsplit(test, "_vs_")) %in% conditions)) {
      stop("run test_diff() with valid contrasts in 'test'",
           ".\nValid contrasts should contain combinations of: '",
           paste0(conditions, collapse = "', '"),
           "', for example '", paste0(conditions[1], "_vs_", conditions[2]),
           "'.", call. = FALSE)
    }
    cntrst <- gsub("_vs_", " - ", test)
  } else if (type == "others") {
    for (i in 1:length(conditions)) {
      design <- cbind(design, ifelse(design[, conditions[i]], 0, 1))
      colnames(design)[ncol(design)] <- paste0("NOT_", conditions[i])
    }
  }

  # Print tested contrasts
  if (type == "others") {
    message("Tested contrasts: ",
            paste(paste0(conditions, "_vs_others"), collapse = ", "))
    limma_res <- data.frame()
    for (c in conditions) {
      sub_design <- design[,c(c, paste0("NOT_",c))]
      fit <- lmFit(raw, design = sub_design)
      made_contrasts <- makeContrasts(contrasts = paste0(c, "-", "NOT_", c), levels = sub_design)
      contrast_fit <- contrasts.fit(fit, made_contrasts)
      eB_fit <- eBayes(contrast_fit)
      temp <- topTable(eB_fit, sort.by = 't', adjust.method="BH", coef = paste0(c, "-", "NOT_", c), number = Inf, confint = T)
      temp <- temp[,c("logFC", "CI.L", "CI.R", "P.Value", "adj.P.Val")]
      colnames(temp) <- c("diff", "CI.L", "CI.R", "p.val", "p.adj")
      colnames(temp) <- paste0(c, "_vs_others_", colnames(temp))
      if (nrow(limma_res) == 0) {
        limma_res <- temp
      } else {
        limma_res <- merge(limma_res, temp, by="row.names")
      }
    }
    rowData(se) <- as.data.frame(left_join(as.data.frame(rowData(se)), limma_res, by=c("ID"="Row.names")))
  } else {
    message("Tested contrasts: ",
            paste(gsub(" - ", "_vs_", cntrst), collapse = ", "))
    # Test for differential expression by empirical Bayes moderation
    # of a linear model on the predefined contrasts
    fit <- lmFit(raw, design = design)
    made_contrasts <- makeContrasts(contrasts = cntrst, levels = design)
    contrast_fit <- contrasts.fit(fit, made_contrasts)

    if (!type %in% c("others")) {
      if(any(is.na(raw))) {
        for(i in cntrst) {
          covariates <- strsplit(i, " - ") %>% unlist
          single_contrast <- makeContrasts(contrasts = i, levels = design[, covariates])
          single_contrast_fit <- contrasts.fit(fit[, covariates], single_contrast)
          contrast_fit$coefficients[, i] <- single_contrast_fit$coefficients[, 1]
          contrast_fit$stdev.unscaled[, i] <- single_contrast_fit$stdev.unscaled[, 1]
        }
      }
    } else {
      if(any(is.na(raw))) {
        for(i in cntrst) {
          covariates <- cntrst
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
    retrieve_fun <- function(comp, fit = eB_fit){
      res <- topTable(fit, sort.by = "t", adjust.method="BH", coef = comp,
                      number = Inf, confint = TRUE)
      res$comparison <- rep(comp, dim(res)[1])
      res <- tibble::rownames_to_column(res)
      return(res)
    }

    #limma_res<- topTable(eB_fit, sort.by = 'B', adjust.method="BH", coef = cntrst, number = Inf, confint = T )
    # limma_res$comparison <- rep(cntrst, dim(limma_res)[1])
    #limma_res <- rownames_to_column(limma_res)
    # Retrieve the differential expression test results
    limma_res <- purrr::map_df(cntrst, retrieve_fun)

    # Select the logFC, CI and qval variables
    table <- limma_res %>%
      dplyr::select(rowname, logFC, CI.L, CI.R, P.Value, adj.P.Val, comparison) %>%
      dplyr::mutate(comparison = gsub(" - ", "_vs_", comparison)) %>%
      tidyr::gather(variable, value, -c(rowname,comparison)) %>%
      dplyr::mutate(variable = dplyr::recode(variable, logFC = "diff", P.Value = "p.val", adj.P.Val = "p.adj")) %>%
      tidyr::unite(temp, comparison, variable) %>%
      tidyr::spread(temp, value)
    rowData(se) <- as.data.frame(left_join(as.data.frame(rowData(se)), table,
                                           by=c("ID"="rowname")))
  }
  return(se)
}

#' Differential expression analysis using limma with fdrtool
#'
#' Performs differential expression analysis on a SummarizedExperiment object
#' using the limma package with fdrtool for FDR estimation. This function uses
#' an alternative FDR calculation method compared to \code{\link{test_limma}}.
#'
#' @param se A \code{SummarizedExperiment} object containing log2-transformed
#'   proteomics data. Must contain 'name' and 'ID' columns in rowData and
#'   'label', 'condition', 'replicate' columns in colData.
#' @param type Character string specifying the type of contrasts to test.
#'   Options are:
#'   \itemize{
#'     \item \code{"control"}: Test all conditions against a specified control.
#'     \item \code{"all"}: Test all pairwise combinations of conditions.
#'     \item \code{"manual"}: Use manually specified contrasts via the \code{test} parameter.
#'   }
#' @param control Character string specifying the control condition.
#'   Required when \code{type = "control"}.
#' @param test Character vector of manual contrasts to test when
#'   \code{type = "manual"}. Contrasts should be in the format
#'   "condition1_vs_condition2".
#' @param design_formula A formula specifying the design matrix.
#'   Default is \code{~ 0 + condition}.
#'
#' @return A \code{SummarizedExperiment} object with differential expression
#'   results added to rowData. For each contrast, the following columns are added:
#'   \itemize{
#'     \item \code{[contrast]_diff}: Log2 fold change
#'     \item \code{[contrast]_CI.L}: Lower confidence interval
#'     \item \code{[contrast]_CI.R}: Upper confidence interval
#'     \item \code{[contrast]_p.val}: Raw p-value
#'     \item \code{[contrast]_p.adj}: Adjusted p-value (q-value from fdrtool)
#'   }
#'
#' @details
#' This function differs from \code{\link{test_limma}} in that it uses the
#' fdrtool package to estimate local false discovery rates, which can provide
#' more accurate FDR estimates in some cases. The fdrtool approach is based on
#' empirical null modeling.
#'
#' @examples
#' \dontrun{
#' # Test all conditions against control
#' se_diff <- test_diff(se, type = "control", control = "Ctrl")
#'
#' # Test all pairwise comparisons
#' se_diff <- test_diff(se, type = "all")
#'
#' # Manual contrast
#' se_diff <- test_diff(se, type = "manual", test = "Treatment_vs_Control")
#' }
#'
#' @seealso \code{\link{test_limma}}, \code{\link{add_rejections}}, \code{\link{plot_volcano}}
#'
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
    sign_df <- cbind(ID = row_data$ID, as.data.frame(sign_df))
    rowData(diff) <- as.data.frame(left_join(as.data.frame(rowData(diff)), sign_df,
      by = c("ID" = "ID")
    ))
  }
  return(diff)
}

#' Create a volcano plot for differential expression results
#'
#' Generates a volcano plot showing log2 fold changes vs. -log10 p-values
#' for a specified contrast from differential expression analysis.
#'
#' @param dep A \code{SummarizedExperiment} object containing differential
#'   expression results from \code{\link{test_limma}} or \code{\link{test_diff}},
#'   with significance annotations from \code{\link{add_rejections}}.
#' @param contrast Character string specifying which contrast to plot.
#'   Must match a contrast name from the differential expression analysis
#'   (e.g., "Treatment_vs_Control").
#' @param label_size Numeric value specifying the size of protein labels.
#'   Default is 3.
#' @param name_col Character string specifying which column from rowData to use
#'   for labeling points. Default is \code{NULL}, which uses "ID".
#' @param add_names Logical indicating whether to add labels for significant
#'   proteins. Default is \code{TRUE}.
#' @param adjusted Logical indicating whether to use adjusted p-values
#'   (\code{TRUE}) or raw p-values (\code{FALSE}). Default is \code{TRUE}.
#' @param plot Logical indicating whether to return a plot (\code{TRUE}) or
#'   a data frame (\code{FALSE}). Default is \code{TRUE}.
#' @param alpha Numeric value for the significance threshold on adjusted/raw
#'   p-values. Default is 0.05.
#' @param lfc Numeric value for the log2 fold change threshold. Default is 1.
#'
#' @return If \code{plot = TRUE}, returns a \code{ggplot} object.
#'   If \code{plot = FALSE}, returns a data frame with columns:
#'   \itemize{
#'     \item \code{protein}: Protein name
#'     \item \code{log2_fold_change}: Log2 fold change
#'     \item \code{p_value_-log10} or \code{adjusted_p_value_-log10}: -log10 transformed p-value
#'     \item \code{signif}: Logical indicating significance
#'   }
#'
#' @examples
#' \dontrun{
#' # Basic volcano plot
#' plot_volcano(se_diff, contrast = "Treatment_vs_Control")
#'
#' # Without labels
#' plot_volcano(se_diff, contrast = "Treatment_vs_Control", add_names = FALSE)
#'
#' # Custom significance thresholds
#' plot_volcano(se_diff, contrast = "Treatment_vs_Control", alpha = 0.01, lfc = 2)
#'
#' # Get data instead of plot
#' volcano_data <- plot_volcano(se_diff, contrast = "Treatment_vs_Control", plot = FALSE)
#' }
#'
#' @seealso \code{\link{test_limma}}, \code{\link{test_diff}}, \code{\link{add_rejections}},
#'   \code{\link{plot_peptide_volcano}}
#'
#' @importFrom ggplot2 ggplot aes geom_vline geom_point geom_text labs theme_bw
#'   theme element_blank element_line scale_color_manual
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr filter arrange select
#'
#' @export
plot_volcano <- function(dep, contrast, label_size = 3, name_col = NULL,
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
  if (is.null(name_col)) {
    name_col <- "ID"
  }
  if (any(!c("name", "ID", name_col) %in% colnames(row_data))) {
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
    name = row_data$name,
    ID = row_data$ID,
    label = row_data[,name_col]
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
      aes(label = label),
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

#' Create a volcano plot for peptide-level differential expression
#'
#' Generates a volcano plot for peptide or site-level data, with options to
#' highlight specific peptides and show related peptides from the same protein.
#'
#' @param dep A \code{SummarizedExperiment} object containing peptide or site-level
#'   differential expression results. Must have metadata level set to "peptide" or "site".
#' @param contrast Character string specifying which contrast to plot.
#' @param peptides Character vector of peptide IDs to highlight in the plot.
#'   Default is \code{NA} (no specific highlighting).
#' @param show_other_peptides Logical indicating whether to show other peptides
#'   from the same protein(s) as the highlighted peptides in blue. Default is \code{TRUE}.
#' @param show_gene Logical indicating whether to display gene names in labels
#'   instead of full peptide IDs. Default is \code{FALSE}.
#' @param label_size Numeric value specifying the size of peptide labels.
#'   Default is 3.
#' @param name_col Character string specifying which column from rowData to use
#'   for labeling points. Default is \code{NULL}, which uses "ID".
#' @param add_names Logical indicating whether to add labels for significant
#'   peptides. Default is \code{TRUE}.
#' @param adjusted Logical indicating whether to use adjusted p-values
#'   (\code{TRUE}) or raw p-values (\code{FALSE}). Default is \code{TRUE}.
#' @param alpha Numeric value for the significance threshold. Default is 0.05.
#' @param lfc Numeric value for the log2 fold change threshold. Default is 1.
#'
#' @return A \code{ggplot} object showing the peptide volcano plot with:
#'   \itemize{
#'     \item Grey points for non-significant peptides
#'     \item Black points for significant peptides
#'     \item Maroon points for specified peptides of interest
#'     \item Blue points for other peptides from the same protein (if \code{show_other_peptides = TRUE})
#'   }
#'
#' @examples
#' \dontrun{
#' # Basic peptide volcano plot
#' plot_peptide_volcano(se_peptide, contrast = "Treatment_vs_Control")
#'
#' # Highlight specific peptides
#' plot_peptide_volcano(se_peptide, contrast = "Treatment_vs_Control",
#'                      peptides = c("PEPTIDE_1", "PEPTIDE_2"))
#'
#' # Show gene names in labels
#' plot_peptide_volcano(se_peptide, contrast = "Treatment_vs_Control",
#'                      peptides = c("PEPTIDE_1"), show_gene = TRUE)
#' }
#'
#' @seealso \code{\link{plot_volcano}}, \code{\link{test_limma}}, \code{\link{add_rejections}}
#'
#' @importFrom ggplot2 ggplot aes geom_vline geom_point geom_text labs theme_bw
#'   theme element_blank element_line scale_color_manual
#' @importFrom ggrepel geom_text_repel
#'
#' @export
plot_peptide_volcano <- function(dep, contrast, peptides=NA, show_other_peptides=T, show_gene=F,
                                 label_size = 3, name_col = NULL,
                                 add_names = TRUE, adjusted = T, alpha = 0.05, lfc = 1) {
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
    metadata(dep)$level %in% c("peptide", "site")
  )

  row_data <- rowData(dep, use.names = FALSE)

  # Show error if inputs do not contain required columns
  if (is.null(name_col)) {
    name_col <- "ID"
  }
  if (any(!c("name", "ID", name_col) %in% colnames(row_data))) {
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
    name = row_data$name,
    ID = row_data$ID,
    label = row_data[,name_col],
    Gene = row_data$Gene
  )

  df <- df_tmp %>%
    data.frame() %>%
    filter(!is.na(signif)) %>%
    arrange(signif)
  name1 <- gsub("_vs_.*", "", contrast)
  name2 <- gsub(".*_vs_", "", contrast)
  # return(df)
  # Plot volcano with or without labels
  if (show_gene) {
    df$ID_new <- paste0(df$Gene, gsub(".*_", "_", df$ID))
    if (!show_other_peptides) {
      p <- ggplot(df, aes(diff, p_values)) +
        geom_vline(xintercept = 0) +
        geom_point(aes(col = signif)) +
        geom_point(data = subset(df, ID %in% c(peptides)), color = "maroon", size= 3) +
        geom_text(data = data.frame(), aes(
          x = c(Inf, -Inf),
          y = c(-Inf, -Inf),
          hjust = c(1, 0),
          vjust = c(-1, -1),
          label = c(name1, name2),
          size = 5,
          fontface = "bold"
        )) +
        geom_text_repel(data = subset(df, ID %in% peptides),
                        color = "maroon",
                        aes(label = ID_new)) +
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
          aes(label = label),
          size = label_size,
          box.padding = unit(0.1, "lines"),
          point.padding = unit(0.1, "lines"),
          segment.size = 0.5
        )
      }
    } else {
      p <- ggplot(df, aes(diff, p_values)) +
        geom_vline(xintercept = 0) +
        geom_point(aes(col = signif)) +
        geom_point(data = subset(df, gsub("_.*", "", ID) %in% c(gsub("_.*", "", peptides))),
                   color = "blue", size= 3) +
        geom_point(data = subset(df, ID %in% peptides),
                   color = "maroon", size= 3) +
        geom_text_repel(data = subset(df, ID %in% peptides),
                        color = "maroon",
                        aes(label = ID_new)) +
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
          aes(label = ID_new),
          size = label_size,
          box.padding = unit(0.1, "lines"),
          point.padding = unit(0.1, "lines"),
          segment.size = 0.5
        )
      }
    }
  } else {
    if (!show_other_peptides) {
      p <- ggplot(df, aes(diff, p_values)) +
        geom_vline(xintercept = 0) +
        geom_point(aes(col = signif)) +
        geom_point(data = subset(df, ID %in% c(peptides)), color = "maroon", size= 3) +
        geom_text_repel(data = subset(df, ID %in% peptides),
                        color = "maroon",
                        aes(label = ID)) +
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
          aes(label = label),
          size = label_size,
          box.padding = unit(0.1, "lines"),
          point.padding = unit(0.1, "lines"),
          segment.size = 0.5
        )
      }
    } else {
      p <- ggplot(df, aes(diff, p_values)) +
        geom_vline(xintercept = 0) +
        geom_point(aes(col = signif)) +
        geom_point(data = subset(df, gsub("_.*", "", ID) %in% c(gsub("_.*", "", peptides))),
                   color = "blue", size= 3) +
        geom_point(data = subset(df, ID %in% peptides),
                   color = "maroon", size= 3) +
        geom_text_repel(data = subset(df, ID %in% peptides),
                        color = "maroon",
                        aes(label = ID)) +
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
          aes(label = label),
          size = label_size,
          box.padding = unit(0.1, "lines"),
          point.padding = unit(0.1, "lines"),
          segment.size = 0.5
        )
      }
    }
  }

  if (adjusted) {
    p <- p + labs(y = expression(-log[10] ~ "Adjusted p-value"))
  } else {
    p <- p + labs(y = expression(-log[10] ~ "P-value"))
  }

  return(p)
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
