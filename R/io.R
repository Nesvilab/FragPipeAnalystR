# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#' Read SummarizedExperiment from RData file
#'
#' Loads a SummarizedExperiment object from an RData file exported from
#' FragPipe-Analyst or saved from a previous R session.
#'
#' @param file Character string specifying the path to the RData file.
#'
#' @return A \code{SummarizedExperiment} object containing the proteomics data
#'   with assay data, rowData, colData, and metadata.
#'
#' @details
#' This function provides a convenient way to load analysis results that were
#' previously saved or exported from the FragPipe-Analyst Shiny application.
#' The first object in the RData file is returned.
#'
#' @examples
#' \dontrun{
#' # Load a previously saved analysis
#' se <- readResultRData("my_analysis.RData")
#'
#' # Check the loaded data
#' se
#' }
#'
#' @seealso \code{\link{make_se_from_files}}
#'
#' @export
readResultRData <- function(file) {
  env <- new.env()
  nm <- load(file, env)[1]
  return(env[[nm]])
}

# make make.unique generates new names start from 1 rather than nothing
make.unique.2 <- function(x, sep = ".") {
  ave(x, x, FUN = function(a) {
    if (length(a) > 1) {
      paste(a, 1:length(a), sep = sep)
    } else {
      a
    }
  })
}

# internal function to read quantification table
readQuantTable <- function(quant_table_path, type = "TMT", level=NULL, log2transform = F, exp_type=NULL, additional_cols=NULL) {
  temp_data <- read.table(quant_table_path,
    header = TRUE,
    fill = TRUE, # to fill any missing data
    sep = "\t",
    quote = "",
    comment.char = "",
    blank.lines.skip = F,
    check.names = F
  )
  colnames(temp_data) <- make.unique.2(colnames(temp_data), "_")
  # validate(maxquant_input_test(temp_data))
  if (type == "TMT") {
    # validate(tmt_input_test(temp_data))
    # convert columns into numeric
    if (is.null(additional_cols)) {
      mut.cols <- colnames(temp_data)[!colnames(temp_data) %in% c("Index", "Gene", "Peptide", "NumberPSM", "ProteinID", "MaxPepProb", "SequenceWindow", "ReferenceIntensity")]
    } else {
      mut.cols <- colnames(temp_data)[!colnames(temp_data) %in% c("Index", "Gene", "Peptide", "NumberPSM", "ProteinID", "MaxPepProb", "SequenceWindow", "ReferenceIntensity", additional_cols)]
    }
    temp_data[mut.cols] <- sapply(temp_data[mut.cols], as.numeric)
  } else if (type == "LFQ") {
    if (level == "peptide") {
      # handle - (dash) in experiment column
      colnames(temp_data) <- gsub("-", ".", colnames(temp_data))
      colnames(temp_data)[colnames(temp_data) == "Protein Description"] <- "Description"
      # validate(fragpipe_input_test(temp_data))
      # remove contam
      temp_data <- temp_data[!grepl("contam", temp_data$Protein),]
      if (!"Modified Sequence" %in% colnames(temp_data)) {
        temp_data$Index <- paste0(temp_data$`Protein ID`, "_", temp_data$`Peptide Sequence`)
      } else { # internally support combined_modified_peptide.tsv
        temp_data$Index <- paste0(temp_data$`Protein ID`, "_", temp_data$`Modified Sequence`)
      }
    } else {
      # handle - (dash) in experiment column
      colnames(temp_data) <- gsub("-", ".", colnames(temp_data))
      # validate(fragpipe_input_test(temp_data))
      # remove contam
      temp_data <- temp_data[!grepl("contam", temp_data$Protein), ]
    }
  } else { # DIA
    if (level == "peptide") {
      if ("SequenceWindow" %in% colnames(temp_data)) {
        print("Error: wrong format. Are you using single-site report? Switch to use level=site to correctly read the file")
        return(NULL)
      } else {
        # validate(fragpipe_DIA_input_test(temp_data))
        # temp_data <- temp_data[!grepl("contam", temp_data$Protein),]
        temp <- melt.data.table(setDT(temp_data[,!colnames(temp_data) %in% c("Proteotypic", "Precursor.Charge", "Precursor.Id",
                                                                             "Modified.Sequence", "First.Protein.Description",
                                                                             "All Mapped Proteins", "All Mapped Genes")]),
                                id.vars = c("Protein.Group", "Protein.Names", "Protein.Ids", "Genes", "Stripped.Sequence"),
                                variable.name = "File.Name")
        temp_data <- as.data.frame(
          dcast.data.table(temp, Protein.Group+Protein.Names+Protein.Ids+Genes+Stripped.Sequence ~ File.Name,
                           value.var = "value", fun.aggregate = function(x) max(x, na.rm=TRUE)))
        temp_data[sapply(temp_data, is.infinite)] <- NA
        temp_data$Index <- paste0(temp_data$Protein.Ids, "_", temp_data$Stripped.Sequence)
        temp_data <- temp_data %>% select(Index, everything())
      }
    } else if (level == "site") {
      if (!"SequenceWindow" %in% colnames(temp_data)) {
        print("Error: no SequenceWindow column found in single-site report.")
        return(NULL)
      }
      if (is.null(additional_cols)) {
        mut.cols <- colnames(temp_data)[!colnames(temp_data) %in% c("Index", "ProteinID", "Gene", "Peptide", "SequenceWindow")]
      } else {
        mut.cols <- colnames(temp_data)[!colnames(temp_data) %in% c("Index", "ProteinID", "Gene", "Peptide", "SequenceWindow", additional_cols)]
      }
      temp_data[mut.cols] <- sapply(temp_data[mut.cols], as.numeric)
    }
  }
  return(temp_data)
}

# internal function to read experiment annotation file
readExpDesign <- function(exp_anno_path, type = "TMT", lfq_type="Intensity", lowercase=F) {
  temp_df <- read.table(exp_anno_path, header = T, sep = "\t", stringsAsFactors = F)
  if (type == "TMT") {
    if (ncol(temp_df) == 1) {
      # submitting annotation.txt (not experiment_annotation.tsv) will crash here
      tryCatch(
        {
          temp_df <- read.table(exp_anno_path,
            header = T,
            sep = " ",
            stringsAsFactors = FALSE
          )
        },
        error = function(e) {
          # validate(need(F,
          #               "Error: coudn't read the experiment_annotation.tsv. Note that experiment_annotation.tsv is not annotation.txt used to denote channel assignment in each plex set."))
        }
      )
    }

    # change it to lower case
    if (lowercase) {
      colnames(temp_df) <- tolower(colnames(temp_df))
    }
    if (!"sample_name" %in% colnames(temp_df)) {
      stop("'sample_name' column is not present in the experiment annotation file.")
    }
    if (!"sample" %in% colnames(temp_df)) {
      stop("'sample' column is not present in the experiment annotation file.")
    }
    if (!is.character(temp_df$sample_name)) {
      temp_df$sample_name <- as.character(temp_df$sample_name)
    }
    if (anyDuplicated(temp_df$sample_name)) {
      stop("Duplicated 'sample_name' detected in the experiment annotation file.")
    }
    if (anyDuplicated(temp_df$sample)) {
      stop("Duplicated 'sample' detected in the experiment annotation file.")
    }
    # to support - (dash) or name starts with number in condition column
    temp_df$condition <- make.names(temp_df$condition)
    # validate(need(try(test_TMT_annotation(temp_df)),
    #               paste0("The input annotation file should have following columns: ",
    #                      "plex, channel, sample, condition, replicate, condition\n",
    #                      "your current input annotation file is with following columns: ", paste(colnames(temp_df), collapse=", "))))
    temp_df$label <- temp_df$sample
    # if duplicate label exists
    if (anyDuplicated(temp_df$label)) {
      # add _number for repeat labels, but need to remove _1
      temp_df$label <- paste(temp_df$label, temp_df$replicate, sep = "_")
      temp_df$label <- gsub("_1$", "", temp_df$label)
      samples_with_replicate <- temp_df$label[grepl("_", temp_df$label)]
      samples_with_replicate <- unique(gsub("_\\d+$", "", samples_with_replicate))
      temp_df[temp_df$label %in% samples_with_replicate, "label"] <- paste0(temp_df[temp_df$label %in% samples_with_replicate, "label"], "_1")
    }
  } else if (type == "LFQ") {
    if (!"sample_name" %in% colnames(temp_df)) {
      stop("'sample_name' column is not present in the experiment annotation file.")
    }
    if (!"sample" %in% colnames(temp_df)) {
      stop("'sample' column is not present in the experiment annotation file.")
    }
    if (!is.character(temp_df$sample_name)) {
      temp_df$sample_name <- as.character(temp_df$sample_name)
    }
    if (anyDuplicated(temp_df$sample_name)) {
      stop("Duplicated 'sample_name' detected in the experiment annotation file.")
    }
    if (anyDuplicated(temp_df$sample)) {
      stop("Duplicated 'sample' detected in the experiment annotation file.")
    }
    # to support - (dash) or name starts with number in condition column
    temp_df$condition <- make.names(temp_df$condition)

    # make sure replicate column is not empty
    if (!all(is.na(temp_df$replicate))) {
      # handle - (dash) in sample (experiment) column
      temp_df$sample <- gsub("-", ".", temp_df$sample)
      temp_df$label <- temp_df$sample
      if (lfq_type == "Intensity") {
        temp_df$label <- paste(temp_df$label, "Intensity", sep = " ")
      } else if (lfq_type == "MaxLFQ") {
        temp_df$label <- paste(temp_df$label, "MaxLFQ.Intensity", sep = " ")
      } else if (lfq_type == "Spectral Count") {
        temp_df$label <- paste(temp_df$label, "Spectral.Count", sep = " ")
      }
    }
  } else if (type == "DIA") {
    temp_df <- read.table(exp_anno_path,
      header = T,
      sep = "\t",
      stringsAsFactors = FALSE
    )
    # change it to lower case
    if (lowercase) {
      colnames(temp_df) <- tolower(colnames(temp_df))
    }
    if (!"file" %in% colnames(temp_df)) {
      stop("'file' column is not present in the experiment annotation file.")
    }
    if (!"sample_name" %in% colnames(temp_df)) {
      stop("'sample_name' column is not present in the experiment annotation file.")
    }
    if (!"sample" %in% colnames(temp_df)) {
      stop("'sample' column is not present in the experiment annotation file.")
    }
    if (!is.character(temp_df$sample_name)) {
      temp_df$sample_name <- as.character(temp_df$sample_name)
    }
    if (anyDuplicated(temp_df$sample_name)) {
      stop("Duplicated 'sample_name' detected in the experiment annotation file.")
    }
    if (anyDuplicated(temp_df$sample)) {
      stop("Duplicated 'sample' detected in the experiment annotation file.")
    }
    # to support - (dash) or name starts with number in condition column
    temp_df$condition <- make.names(temp_df$condition)
    # make sure replicate column is not empty
    if (!all(is.na(temp_df$replicate))) {
      temp_df$label <- temp_df$file
    }
  }
  return(temp_df)
}

#' Create SummarizedExperiment from FragPipe output files
#'
#' Reads FragPipe quantification output and experiment annotation files to
#' create a SummarizedExperiment object for downstream analysis.
#'
#' @param quant_table_path Character string specifying the path to the
#'   quantification table (e.g., gene/protein/peptide-level abundance report).
#' @param exp_anno_path Character string specifying the path to the experiment
#'   annotation file (experiment_annotation.tsv).
#' @param type Character string specifying the quantification method.
#'   Options are:
#'   \itemize{
#'     \item \code{"TMT"}: TMT-based quantification (default)
#'     \item \code{"LFQ"}: Label-free quantification
#'     \item \code{"DIA"}: Data-independent acquisition
#'   }
#' @param level Character string specifying the quantification level.
#'   Options are:
#'   \itemize{
#'     \item \code{"gene"}: Gene-level (default for TMT)
#'     \item \code{"protein"}: Protein-level (default for LFQ/DIA)
#'     \item \code{"peptide"}: Peptide-level
#'     \item \code{"site"}: PTM site-level
#'     \item \code{"glycan"}: Glycan-level
#'   }
#' @param exp_type Character string specifying the experiment type for
#'   specialized analyses. Options include "global", "phospho", "glyco",
#'   "acetyl", "ubiquit". Default is \code{NULL}.
#' @param log2transform Logical indicating whether to log2-transform intensity
#'   values. Default is \code{NULL} (auto-determined based on type).
#' @param lfq_type Character string specifying which LFQ intensity to use.
#'   Options are "Intensity" (default), "MaxLFQ", or "Spectral Count".
#' @param gencode Logical indicating whether to filter for GENCODE/Ensembl
#'   identifiers (DIA only). Default is \code{FALSE}.
#' @param additional_cols Character vector of additional column names to
#'   preserve in rowData. Default is \code{NULL}.
#'
#' @return A \code{SummarizedExperiment} object with:
#'   \itemize{
#'     \item \code{assay}: Log2-transformed intensity matrix
#'     \item \code{rowData}: Feature annotations (gene, protein, peptide info)
#'     \item \code{colData}: Sample metadata (condition, replicate, etc.)
#'     \item \code{metadata}: Experiment type, level, and other settings
#'   }
#'
#' @details
#' This is the primary function for importing FragPipe output into R.
#' The function automatically handles different file formats based on the
#' specified type and level. Zero values are converted to NA, and
#' log2 transformation is applied where appropriate.
#'
#' @examples
#' \dontrun{
#' # Import TMT data
#' se <- make_se_from_files("gene_abundance.tsv", "experiment_annotation.tsv",
#'                          type = "TMT", level = "gene")
#'
#' # Import LFQ data with MaxLFQ intensities
#' se <- make_se_from_files("combined_protein.tsv", "experiment_annotation.tsv",
#'                          type = "LFQ", lfq_type = "MaxLFQ")
#'
#' # Import DIA peptide data
#' se <- make_se_from_files("peptide_report.tsv", "experiment_annotation.tsv",
#'                          type = "DIA", level = "peptide")
#'
#' # Import phosphoproteomics data
#' se <- make_se_from_files("site_abundance.tsv", "experiment_annotation.tsv",
#'                          type = "TMT", level = "site", exp_type = "phospho")
#' }
#'
#' @seealso \code{\link{make_unique}}, \code{\link{make_se_customized}},
#'   \code{\link{readResultRData}}
#'
#' @importFrom SummarizedExperiment SummarizedExperiment colData colData<-
#'   rowData rowData<- metadata metadata<-
#'
#' @export
make_se_from_files <- function(quant_table_path, exp_anno_path, type = "TMT", level = NULL, exp_type=NULL,
                               log2transform = NULL, lfq_type = "Intensity", gencode = F, additional_cols=NULL) {
  if (type == "TMT" & is.null(level)) {
    level <- "gene"
  } else if (is.null(level)) {
    level <- "protein"
  }

  if (type == "DIA" & is.null(log2transform)) {
    log2transform <- T
  } else if (is.null(level)) {
    log2transform <- F
  }

  if (!level %in% c("gene", "protein", "peptide", "site", "glycan")) {
    cat(paste0("The specified level: ", level, " is not a valid level. Available levels are gene, protein, site, and peptide.\n"))
    return(NULL)
  }

  quant_table <- readQuantTable(quant_table_path, type = type, level=level, exp_type=exp_type, additional_cols=additional_cols)
  exp_design <- readExpDesign(exp_anno_path, type = type, lfq_type = lfq_type)
  if (type == "LFQ") {
    quant_table <- quant_table[!grepl("contam", quant_table$Protein),]
    if (level != "peptide") {
      data_unique <- make_unique(quant_table, "Gene", "Protein ID")
      if (lfq_type == "Intensity") {
        lfq_columns <- setdiff(
          grep("Intensity", colnames(data_unique)),
          grep("MaxLFQ", colnames(data_unique))
        )
        lfq_columns <- setdiff(lfq_columns, grep("Total Intensity", colnames(data_unique)))
        lfq_columns <- setdiff(lfq_columns, grep("Unique Intensity", colnames(data_unique)))
      } else if (lfq_type == "MaxLFQ") {
        lfq_columns <- grep("MaxLFQ", colnames(data_unique))
        if (length(lfq_columns) == 0) {
          stop(safeError("No MaxLFQ column available. Please make sure your file has MaxLFQ intensity columns."))
        }
      } else if (lfq_type == "Spectral Count") {
        lfq_columns <- grep("Spectral", colnames(data_unique))
        lfq_columns <- setdiff(lfq_columns, grep("Total Spectral Count", colnames(data_unique)))
        lfq_columns <- setdiff(lfq_columns, grep("Unique Spectral Count", colnames(data_unique)))
      }

      ## Check for matching columns in expression report and experiment manifest file
      # test_match_lfq_column_manifest(data_unique, lfq_columns, exp_design)
      if (lfq_type == "Spectral Count") {
        data_se <- make_se_customized(data_unique, lfq_columns, exp_design, log2transform = F,
                                      exp="LFQ", lfq_type = lfq_type)
      } else {
        data_se <- make_se_customized(data_unique, lfq_columns, exp_design, log2transform = T,
                                      exp="LFQ", lfq_type = lfq_type)
      }
    } else {
      data_unique <- make_unique(quant_table, "Protein ID", "Index")
      if (lfq_type == "Intensity") {
        lfq_columns <- setdiff(grep("Intensity", colnames(data_unique)),
                               grep("MaxLFQ", colnames(data_unique)))
        lfq_columns <- setdiff(lfq_columns, grep("Total Intensity", colnames(data_unique)))
        lfq_columns <- setdiff(lfq_columns, grep("Unique Intensity", colnames(data_unique)))
      } else if (lfq_type == "MaxLFQ") {
        lfq_columns <- grep("MaxLFQ", colnames(data_unique))
        if (length(lfq_columns) == 0) {
          stop(safeError("No MaxLFQ column available. Please make sure your files have MaxLFQ intensity columns."))
        }
      } else if (lfq_type == "Spectral Count") {
        lfq_columns <- grep("Spectral", colnames(data_unique))
        lfq_columns <- setdiff(lfq_columns, grep("Total Spectral Count", colnames(data_unique)))
        lfq_columns <- setdiff(lfq_columns, grep("Unique Spectral Count", colnames(data_unique)))
      }
      # TODO: Check for matching columns in expression report and experiment manifest file
      # test_match_lfq_column_manifest(data_unique, lfq_columns, exp_design)
      if (lfq_type == "Spectral Count") {
        data_se <- make_se_customized(data_unique, lfq_columns, exp_design, log2transform=F, exp="LFQ", lfq_type="Spectral Count", level="peptide")
      } else {
        data_se <- make_se_customized(data_unique, lfq_columns, exp_design, log2transform=T, exp="LFQ", lfq_type=lfq_type, level="peptide")
      }
    }
  } else if (type == "DIA") {
    if (level == "protein") {
      if (gencode) {
        quant_table <- quant_table[grepl("^ENS", quant_table$Protein.Group),]
      }
      data_unique <- make_unique(quant_table, "Genes", "Protein.Group")
      cols <- colnames(data_unique)
      if (is.null(additional_cols)) {
        selected_cols <- which(!(cols %in% c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description", "ID", "name")))
      } else {
        selected_cols <- which(!(cols %in% c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description", "ID", "name",
                                             additional_cols)))
      }
      # TODO: use DIA function
      # test_match_DIA_column_design(data_unique, selected_cols, exp_design)
      data_se <- make_se_customized(data_unique, selected_cols, exp_design,
                                    log2transform = log2transform, exp="DIA", level="protein")
      dimnames(data_se) <- list(dimnames(data_se)[[1]], colData(data_se)$sample_name)
      colData(data_se)$label <- colData(data_se)$sample_name
    } else if (level == "gene") {
      if (gencode) {
        quant_table <- quant_table[grepl("^ENS", quant_table$Genes),]
      }
      quant_table$Index <- quant_table$Genes
      data_unique <- make_unique(quant_table, "Genes", "Index")
      cols <- colnames(data_unique)
      if (is.null(additional_cols)) {
        selected_cols <- which(!(cols %in% c("Genes", "Index", "ID", "name")))
      } else {
        selected_cols <- which(!(cols %in% c("Genes", "Index", "ID", "name", additional_cols)))
      }
      # TODO: use DIA function
      # test_match_DIA_column_design(data_unique, selected_cols, exp_design)
      data_se <- make_se_customized(data_unique, selected_cols, exp_design,
                                    log2transform = log2transform, exp="DIA", level="gene")
      dimnames(data_se) <- list(dimnames(data_se)[[1]], colData(data_se)$sample_name)
      colData(data_se)$label <- colData(data_se)$sample_name
    } else {
      if (level == "site") { # single-site
        data_unique <- make_unique(quant_table, "ProteinID", "Index")
      } else { # level == "peptide"
        data_unique <- make_unique(quant_table, "Protein.Group", "Index")
      }
      cols <- colnames(data_unique)
      if (is.null(additional_cols)) {
        selected_cols <- which(!(cols %in% c("Index", "Protein.Group", "Protein.Ids", "Stripped.Sequence", "Protein.Names", "Genes", "First.Protein.Description", "ID", "name",
                                             "Gene", "ProteinID", "Peptide", "SequenceWindow", "All Mapped Proteins", "All Mapped Genes")))
      } else {
        selected_cols <- which(!(cols %in% c("Index", "Protein.Group", "Protein.Ids", "Stripped.Sequence", "Protein.Names", "Genes", "First.Protein.Description", "ID", "name",
                                             "Gene", "ProteinID", "Peptide", "SequenceWindow", "All Mapped Proteins", "All Mapped Genes", additional_cols)))
      }
      # test_match_DIA_column_design(data_unique, selected_cols, exp_design)
      data_se <- make_se_customized(data_unique, selected_cols, exp_design, log2transform=log2transform, exp="DIA", level=level)
      dimnames(data_se) <- list(dimnames(data_se)[[1]], colData(data_se)$sample_name)
      colData(data_se)$label <- colData(data_se)$sample_name
    }
    if (level == "peptide") {
      rowData(data_se)$Gene <- rowData(data_se)$Genes
    }
  } else { # TMT
    temp_exp_design <- exp_design
    # sample without specified condition will be removed
    temp_exp_design <- temp_exp_design[!is.na(temp_exp_design$condition), ]
    temp_exp_design <- temp_exp_design[!temp_exp_design$condition == "", ]
    # temp_exp_design[is.na(temp_exp_design), "replicate"] <- 1
    # need to handle duplicate columns first
    # filtered_data <- avearrays(filtered_data)
    # filtered_data <- as.data.frame(filtered_data)
    # print(apply(filtered_data, 2, is.numeric))

    # ProteinID is included as Index for protein-level result
    if (level == "protein") {
      quant_table$ProteinID <- quant_table$Index
    }
    data_unique <- make_unique(quant_table, "ProteinID", "Index")
    # handle unmatched columns
    overlapped_samples <- intersect(colnames(data_unique), temp_exp_design$label)
    if (level == "gene") {
      interest_cols <- c("Index", "NumberPSM", "ProteinID", "MaxPepProb", "ReferenceIntensity", "name", "ID")
      if (!is.null(additional_cols)) {
        interest_cols <- c(interest_cols, additional_cols)
      }
      data_unique <- data_unique[, colnames(data_unique) %in% c(interest_cols, overlapped_samples)]
      temp_exp_design <- temp_exp_design[temp_exp_design$label %in% overlapped_samples, ]
      cols <- colnames(data_unique)
      selected_cols <- which(!(cols %in% interest_cols))
    } else if (level == "protein") {
      interest_cols <- c("Index", "NumberPSM", "Gene", "ProteinID", "MaxPepProb", "ReferenceIntensity", "name", "ID")
      if (!is.null(additional_cols)) {
        interest_cols <- c(interest_cols, additional_cols)
      }
      data_unique <- data_unique[, colnames(data_unique) %in% c(interest_cols, overlapped_samples)]
      temp_exp_design <- temp_exp_design[temp_exp_design$label %in% overlapped_samples, ]
      cols <- colnames(data_unique)
      selected_cols <- which(!(cols %in% interest_cols))
    } else if (level == "peptide" | level == "site") {
      interest_cols <- c("Index", "Gene", "Peptide", "NumberPSM", "ProteinID", "SequenceWindow", "MaxPepProb", "ReferenceIntensity", "name", "ID")
      if (!is.null(additional_cols)) {
        interest_cols <- c(interest_cols, additional_cols)
      }
      data_unique <- data_unique[, colnames(data_unique) %in% c(interest_cols, overlapped_samples)]
      temp_exp_design <- temp_exp_design[temp_exp_design$label %in% overlapped_samples, ]
      cols <- colnames(data_unique)
      selected_cols <- which(!(cols %in% interest_cols))
    } else {
      interest_cols <- c("Index", "Gene", "ProteinID", "Peptide", "SequenceWindow", "Start", "End", "MaxPepProb", "ReferenceIntensity", "name", "ID", "Spectrum Number")
      if (!is.null(additional_cols)) {
        interest_cols <- c(interest_cols, additional_cols)
      }
      data_unique <- data_unique[, colnames(data_unique) %in% c(interest_cols, overlapped_samples)]
      temp_exp_design <- temp_exp_design[temp_exp_design$label %in% overlapped_samples, ]
      cols <- colnames(data_unique)
      selected_cols <- which(!(cols %in% interest_cols))
    }
    data_unique[selected_cols] <- apply(data_unique[selected_cols], 2, as.numeric)

    # test_match_tmt_column_design(data_unique, selected_cols, temp_exp_design)
    # TMT-I report is already log2 transformed
    data_se <- make_se_customized(data_unique, selected_cols, temp_exp_design, exp="TMT", level=level)
  }
  supported_types <- c("global", "phospho", "glyco", "acetyl", "ubiquit")
  if (!is.null(exp_type)) {
    if (!(exp_type %in% supported_types)) {
      cat(paste0(exp_type), " is not a customized type. Some functionality might not work.")
    }
  }
  metadata(data_se)$exp_type <- exp_type
  metadata(data_se)$level <- level
  return(data_se)
}


#' Generate unique feature identifiers
#'
#' Creates unique identifiers for a proteomics dataset based on specified
#' name and ID columns. This is required before creating a SummarizedExperiment
#' object to ensure each row has a unique identifier.
#'
#' @param proteins A data.frame containing the protein/peptide data table.
#' @param names Character string specifying the column name containing
#'   feature names (e.g., "Gene" or "Protein ID").
#' @param ids Character string specifying the column name containing
#'   feature identifiers (e.g., "Protein.Group" or "Index").
#' @param delim Character string specifying the delimiter used to separate
#'   multiple names within a single entry. Default is ";".
#'
#' @return A data.frame with two additional columns:
#'   \itemize{
#'     \item \code{name}: Unique feature names (uses ID if name is missing)
#'     \item \code{ID}: Feature identifiers (same as ids column)
#'   }
#'   Names are made unique by appending numbers when duplicates exist.
#'
#' @details
#' This function is adapted from the DEP package. When duplicate names are
#' encountered, they are made unique by appending ".1", ".2", etc.
#' If a name is empty or NA, the ID is used instead.
#'
#' @examples
#' \dontrun{
#' # Make unique identifiers for protein data
#' proteins_unique <- make_unique(protein_table, "Gene", "Protein.Group")
#'
#' # Make unique identifiers for peptide data
#' peptides_unique <- make_unique(peptide_table, "Protein ID", "Index")
#' }
#'
#' @seealso \code{\link{make_se_from_files}}, \code{\link{make_se_customized}}
#'
#' @importFrom tibble is_tibble
#' @importFrom dplyr mutate
#'
#' @export
make_unique <- function(proteins, names, ids, delim = ";") {
  # Show error if inputs are not the required classes
  # assertthat::assert_that(is.data.frame(proteins),
  #                         is.character(names),
  #                         length(names) == 1,
  #                         is.character(ids),
  #                         length(ids) == 1,
  #                         is.character(delim),
  #                         length(delim) == 1)

  col_names <- colnames(proteins)
  # Show error if inputs do not contain required columns
  if (!names %in% col_names) {
    stop("'", names, "' is not a column in '",
      deparse(substitute(proteins)), "'",
      call. = FALSE
    )
  }
  if (!ids %in% col_names) {
    stop("'", ids, "' is not a column in '",
      deparse(substitute(proteins)), "'",
      call. = FALSE
    )
  }

  # If input is a tibble, convert to data.frame
  if (tibble::is_tibble(proteins)) {
    proteins <- as.data.frame(proteins)
  }

  # Select the name and id columns, and check for NAs
  double_NAs <- apply(proteins[, c(names, ids)], 1, function(x) all(is.na(x)))
  if (any(double_NAs)) {
    stop("NAs in both the 'names' and 'ids' columns")
  }

  # Take the first identifier per row and make unique names.
  # If there is no name, the ID will be taken.
  proteins_unique <- proteins %>%
    dplyr::mutate(
      name = get(names),
      ID = get(ids),
      name = make.unique(ifelse(name == "" | is.na(name), ID, name))
    )
  return(proteins_unique)
}

#' Create SummarizedExperiment from data.frame
#'
#' Converts a protein/peptide data.frame and experimental design into a
#' SummarizedExperiment object. This is the core function for creating
#' analysis-ready data structures.
#'
#' @param proteins_unique A data.frame with unique feature identifiers,
#'   typically output from \code{\link{make_unique}()}. Must contain
#'   'name' and 'ID' columns.
#' @param columns Integer vector specifying the column indices containing
#'   the intensity/abundance data.
#' @param expdesign A data.frame containing the experimental design with
#'   required columns: 'label', 'condition', 'replicate', and 'sample_name'.
#' @param log2transform Logical indicating whether to log2-transform the
#'   intensity values. Default is \code{FALSE}.
#' @param exp Character string specifying the quantification method.
#'   Options are "LFQ", "TMT", or "DIA". Default is "LFQ".
#' @param lfq_type Character string specifying the LFQ intensity type
#'   (for LFQ experiments). Default is \code{NULL}.
#' @param level Character string specifying the quantification level
#'   (e.g., "protein", "peptide", "site"). Default is \code{NULL}.
#' @param exp_type Character string specifying the experiment type
#'   (e.g., "global", "phospho"). Default is \code{NULL}.
#'
#' @return A \code{SummarizedExperiment} object containing:
#'   \itemize{
#'     \item \code{assay}: Intensity matrix (optionally log2-transformed)
#'     \item \code{rowData}: Feature annotations
#'     \item \code{colData}: Sample metadata from expdesign
#'     \item \code{metadata}: Analysis settings (exp, level, exp_type, etc.)
#'   }
#'   Zero values in the intensity data are converted to NA.
#'
#' @details
#' This function is adapted from the DEP package. It performs the following:
#' \enumerate{
#'   \item Extracts intensity columns based on indices
#'   \item Converts zeros to NA
#'   \item Optionally log2-transforms the data
#'   \item Matches samples between data and experimental design
#'   \item Creates the SummarizedExperiment structure
#' }
#'
#' @examples
#' \dontrun{
#' # Create SE from prepared data
#' se <- make_se_customized(proteins_unique, columns = 5:20,
#'                          expdesign = exp_design, log2transform = TRUE,
#'                          exp = "LFQ", level = "protein")
#' }
#'
#' @seealso \code{\link{make_unique}}, \code{\link{make_se_from_files}}
#'
#' @importFrom SummarizedExperiment SummarizedExperiment rowData rowData<-
#' @importFrom tibble is_tibble
#'
#' @export
make_se_customized <- function(proteins_unique, columns, expdesign, log2transform = F, exp="LFQ", lfq_type=NULL,
                               level=NULL, exp_type=NULL) {
  # Show error if inputs are not the required classes
  # assertthat::assert_that(is.data.frame(proteins_unique),
  #                         is.integer(columns),
  #                         is.data.frame(expdesign))

  # Show error if inputs do not contain required columns
  if (any(!c("name", "ID") %in% colnames(proteins_unique))) {
    stop("'name' and/or 'ID' columns are not present in '",
      deparse(substitute(proteins_unique)),
      "'.\nRun make_unique() to obtain the required columns",
      call. = FALSE
    )
  }
  if (any(!c("label", "condition", "replicate") %in% colnames(expdesign))) {
    stop("'label', 'condition' and/or 'replicate' columns",
      "are not present in the experimental design",
      call. = FALSE
    )
  }
  if (any(!apply(proteins_unique[, columns], 2, is.numeric))) {
    stop("specified 'columns' should be numeric",
      "\nRun make_se_customized() with the appropriate columns as argument",
      call. = FALSE
    )
  }

  # If input is a tibble, convert to data.frame
  if (tibble::is_tibble(proteins_unique)) {
    proteins_unique <- as.data.frame(proteins_unique)
  }
  if (tibble::is_tibble(expdesign)) {
    expdesign <- as.data.frame(expdesign)
  }

  # Select the assay data
  rownames(proteins_unique) <- proteins_unique$ID
  raw <- proteins_unique[, columns]
  raw[raw == 0] <- NA
  if (log2transform) {
    raw <- log2(raw)
  }
  # Generate the colData from the experimental design
  # and match these with the assay data
  # expdesign <- mutate(expdesign, condition = make.names(condition)) %>%
  #   unite(ID, condition, replicate, remove = FALSE)
  # rownames(expdesign) <- expdesign$ID
  rownames(expdesign) <- expdesign$label
  # print(expdesign)
  # print(colnames(raw))
  matched <- match(
    make.names(expdesign$label),
    make.names(colnames(raw))
  )

  if (any(is.na(matched))) {
    print(make.names(expdesign$label))
    print(make.names(colnames(raw)))
    stop(
      "None of the labels in the experimental design match ",
      "with column names in 'proteins_unique'",
      "\nRun make_se() with the correct labels in the experimental design",
      "and/or correct columns specification"
    )
  }

  rownames(expdesign) <- expdesign$sample_name
  colnames(raw)[matched] <- expdesign$sample_name
  raw <- raw[, !is.na(colnames(raw))][rownames(expdesign)]

  # Select the rowData
  row_data <- proteins_unique[, -columns]
  rownames(row_data) <- row_data$ID
  # Generate the SummarizedExperiment object
  se <- SummarizedExperiment(
    assays = as.matrix(raw),
    colData = expdesign,
    rowData = row_data,
    metadata = list("log2transform"=log2transform, "exp"=exp, "lfq_type"=lfq_type,
                    "exp_type"=exp_type, "level"=level)
  )

  if (exp == "DIA" & !is.null(level)) {
    if (level == "protein") {
      rowData(se)$Index <- rowData(se)$Protein.Group
    }
  }

  return(se)
}

export_se <- function(se, output_file, sep="\t", include_cols=NULL) {
  if (is.null(include_cols)) {
    temp <- cbind(rowData(se), assay(se))
  } else {
    temp <- cbind(rowData(se)[,include_cols], assay(se))
  }
  write.table(temp, output_file, sep=sep, quote=F, row.names = F)
}
