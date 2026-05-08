#' Gene Set Enrichment Analysis (GSEA)
#'
#' Performs Gene Set Enrichment Analysis on ranked protein-level data using
#' the clusterProfiler package.
#'
#' @param se A \code{SummarizedExperiment} object containing protein-level data.
#'   Peptide-level data is not supported.
#' @param col Character string specifying the column name in rowData to use
#'   for ranking (e.g., log2 fold change column like "Treatment_vs_Control_diff").
#' @param database Character string specifying the database for enrichment analysis.
#'   Options are:
#'   \itemize{
#'     \item \code{"GO Biological Process"}: Gene Ontology Biological Process (default)
#'     \item \code{"GO Cellular Component"}: Gene Ontology Cellular Component
#'     \item \code{"GO Molecular Function"}: Gene Ontology Molecular Function
#'     \item \code{"Hallmark"}: MSigDB Hallmark gene sets (requires msigdbr package)
#'     \item \code{"KEGG"}: KEGG pathways
#'   }
#' @param file Character string specifying the path to a custom GMT file.
#'   If provided, this overrides the \code{database} parameter.
#' @param minGSSize Numeric value for the minimum gene set size to include.
#'   Default is 0.
#' @param eps Numeric value for the boundary for calculating p-values.
#'   Default is 0. See \code{fgsea} documentation.
#' @param convert Logical indicating whether to convert the result to a data frame
#'   and add a Count column. Default is \code{TRUE}.
#'
#' @return If \code{convert = TRUE}, returns a data frame with GSEA results including:
#'   \itemize{
#'     \item \code{ID}: Gene set identifier
#'     \item \code{Description}: Gene set description
#'     \item \code{NES}: Normalized Enrichment Score
#'     \item \code{pvalue}: P-value
#'     \item \code{p.adjust}: Adjusted p-value (BH method)
#'     \item \code{core_enrichment}: Leading edge genes
#'     \item \code{Count}: Number of genes in core enrichment
#'   }
#'   If \code{convert = FALSE}, returns a clusterProfiler gseaResult object.
#'
#' @examples
#' \dontrun{
#' # GSEA using GO Biological Process
#' gsea_result <- GSEA_test(se_diff, col = "Treatment_vs_Control_diff",
#'                          database = "GO Biological Process")
#'
#' # GSEA using KEGG pathways
#' gsea_result <- GSEA_test(se_diff, col = "Treatment_vs_Control_diff",
#'                          database = "KEGG")
#'
#' # GSEA using custom GMT file
#' gsea_result <- GSEA_test(se_diff, col = "Treatment_vs_Control_diff",
#'                          file = "my_genesets.gmt")
#' }
#'
#' @seealso \code{\link{plot_GSEA}}, \code{\link{or_test}}
#'
#' @importFrom clusterProfiler GSEA gseGO gseKEGG bitr read.gmt setReadable
#' @importFrom stringr str_count
#'
#' @export
GSEA_test <- function(se, col=NULL, database="GO Biological Process", file=NULL, minGSSize=0, eps=0, convert=T) {
  if (metadata(se)$level %in% c("peptide")) {
    print("Error: Currently, GSEA test doesn't support peptide level data.")
    return(NULL)
  } else if (is.null(col)) {
    cat("Error: No column specified. \n")
    return(NULL)
  } else if (!col %in% colnames(rowData(se))) {
    print(paste0("Error: The column speiciifed: ", col, " doesn't not exist. Valid columns in rowData are ", paste0(colnames(rowData(se)), collapse = ", ")))
    return(NULL)
  }
  geneList <- rowData(se)[[col]]
  names(geneList) <- rowData(se)[["ID"]]
  geneList <- geneList[!is.na(geneList)]
  geneList = sort(geneList, decreasing = TRUE)

  if (!is.null(file)) {
    term2gene <- read.gmt(file)
    result <- GSEA(geneList, TERM2GENE=term2gene, pvalueCutoff = 1, minGSSize = minGSSize, eps=eps)
  } else {
    database_mappings <- c(
      "GO Biological Process" = "BP",
      "GO Cellular Component" = "CC",
      "GO Molecular Function" = "MF"
    )
    if (database %in% c("GO Biological Process", "GO Cellular Component", "GO Molecular Function")){
      result <- gseGO(geneList = geneList,
            OrgDb = org.Hs.eg.db,
            ont = database_mappings[database],
            minGSSize = minGSSize,
            maxGSSize = 500,
            pvalueCutoff = 1,
            eps = eps,
            verbose = F)
    } else if (database == "Hallmark") {
      if(!requireNamespace("msigdbr", quietly = T)) {
        cat("msigdbr package is required to preform the task.")
        return(NULL)
      }
      hallmark <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
        dplyr::select(gs_name, gene_symbol)
      result <- GSEA(geneList, TERM2GENE=hallmark, pvalueCutoff = 1, minGSSize = minGSSize, eps=eps)
    } else if (database == "KEGG") {
      mappings <- bitr(names(geneList), fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
      names(geneList) <- mappings$ENTREZID
      result <- gseKEGG(geneList     = geneList,
                        organism     = "hsa",
                        keyType = "ncbi-geneid",
                        minGSSize    = minGSSize,
                        pvalueCutoff = 1,
                        verbose      = F)
      result <- setReadable(result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    # } else if (database == "Reactome") {
    } else {
      cat(paste0("The database specified: ", database, " is not supported."))
      return(NULL)
    }
  }
  if (convert){
    result <- as.data.frame(result)
    result$Count <- str_count(result$core_enrichment, "/")
  }
  return(result)
}

#' Plot GSEA results
#'
#' Creates a dot plot visualization of Gene Set Enrichment Analysis results,
#' showing the most positively and negatively enriched gene sets.
#'
#' @param gsea_result A data frame containing GSEA results from \code{\link{GSEA_test}}.
#'   Must contain columns: ID, NES, p.adjust, core_enrichment.
#' @param categroies Numeric value specifying the number of top and bottom
#'   enriched gene sets to display. Default is 15 (showing 30 total: 15 positive
#'   and 15 negative NES).
#'
#' @return A \code{ggplot} object showing a dot plot where:
#'   \itemize{
#'     \item X-axis: Normalized Enrichment Score (NES)
#'     \item Y-axis: Gene set names
#'     \item Point size: Number of core enrichment genes
#'     \item Point color: Adjusted p-value (red = low, blue = high)
#'   }
#'
#' @examples
#' \dontrun{
#' # Perform GSEA
#' gsea_result <- GSEA_test(se_diff, col = "Treatment_vs_Control_diff")
#'
#' # Plot top 10 positive and negative enrichments
#' plot_GSEA(gsea_result, categroies = 10)
#' }
#'
#' @seealso \code{\link{GSEA_test}}
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_continuous theme_bw
#' @importFrom stringr str_count
#'
#' @export
plot_GSEA <- function(gsea_result, categroies=15) {
  if (!"Count" %in% colnames(gsea_result)) {
    gsea_result$Count <- str_count(gsea_result$core_enrichment, "/")
  }
  temp <- gsea_result[order(gsea_result$NES, decreasing = T),]
  ID_selected <- temp[c(1:categroies, (dim(temp)[1]-categroies + 1):dim(temp)[1]),"ID"]
  temp <- temp[temp$ID %in% ID_selected,]
  temp$ID <- factor(temp$ID, levels = rev(ID_selected))
  plot <- ggplot(temp, aes(x=NES, y=ID, size=Count, color=p.adjust)) +
    geom_point() +
    scale_color_continuous(low="red", high="blue", name = "p.adjust",
                           guide=guide_colorbar(reverse=T,
                                                label.vjust = 0.5)) +
    theme_bw()
  return(plot)
}

#' Over-representation analysis (ORA)
#'
#' Performs over-representation analysis on significantly differential proteins
#' using either Enrichr or clusterProfiler as the backend.
#'
#' @param se A \code{SummarizedExperiment} object containing differential expression
#'   results from \code{\link{test_limma}} or \code{\link{test_diff}}.
#' @param database Character string specifying the database for enrichment analysis.
#'   Options are:
#'   \itemize{
#'     \item \code{"GO Biological Process"}: Gene Ontology Biological Process (default)
#'     \item \code{"GO Cellular Component"}: Gene Ontology Cellular Component
#'     \item \code{"GO Molecular Function"}: Gene Ontology Molecular Function
#'     \item \code{"Hallmark"}: MSigDB Hallmark gene sets
#'     \item \code{"KEGG"}: KEGG pathways (human)
#'     \item \code{"KEGG (Mouse)"}: KEGG pathways for mouse (clusterProfiler backend only)
#'     \item \code{"WikiPathways"}: WikiPathways (human, clusterProfiler backend only)
#'     \item \code{"WikiPathways (Mouse)"}: WikiPathways for mouse (clusterProfiler backend only)
#'     \item \code{"Reactome"}: Reactome pathways (Enrichr backend only)
#'   }
#' @param backend Character string specifying the analysis backend.
#'   Options are:
#'   \itemize{
#'     \item \code{"enrichr"}: Use Enrichr web service (default)
#'     \item \code{"clusterProfiler"}: Use clusterProfiler R package
#'   }
#' @param direction Character string specifying which proteins to analyze.
#'   Options are:
#'   \itemize{
#'     \item \code{"UP"}: Up-regulated proteins (positive log2 fold change)
#'     \item \code{"DOWN"}: Down-regulated proteins (negative log2 fold change)
#'   }
#' @param log2_threshold Numeric value for the log2 fold change cutoff.
#'   Default is 0.7 (approximately 1.6-fold change).
#' @param alpha Numeric value for the significance threshold on adjusted p-values.
#'   Default is 0.05.
#'
#' @return A data frame with enrichment results including:
#'   \itemize{
#'     \item \code{Term}: Enriched term name
#'     \item \code{Overlap}: Number of genes overlapping with the term
#'     \item \code{P.value}: Raw p-value
#'     \item \code{Adjusted.P.value}: Adjusted p-value
#'     \item \code{Odds.Ratio}: Odds ratio for enrichment
#'     \item \code{log_odds}: Log2 odds ratio (background corrected)
#'     \item \code{contrast}: The contrast tested
#'     \item \code{var}: Database used
#'     \item \code{p_hyper}: Hypergeometric test p-value (background corrected)
#'     \item \code{p.adjust_hyper}: Adjusted hypergeometric p-value
#'   }
#'
#' @examples
#' \dontrun
#' # ORA for up-regulated proteins using Enrichr
#' or_result <- or_test(se_diff, database = "GO Biological Process",
#'                      direction = "UP", alpha = 0.05)
#'
#' # ORA using clusterProfiler backend
#' or_result <- or_test(se_diff, database = "KEGG",
#'                      backend = "clusterProfiler")
#' }
#'
#' @seealso \code{\link{plot_or}}, \code{\link{GSEA_test}}
#'
#' @importFrom clusterProfiler enricher enrichGO enrichKEGG enrichWP bitr setReadable
#'
#' @export
or_test <- function(se, database="GO Biological Process", backend="enrichr", direction="UP", log2_threshold=0.7, alpha=0.05) {
  if (backend == "enrichr") {
    database_mappings <- c(
      "GO Biological Process" = "GO_Biological_Process_2021",
      "GO Cellular Component" = "GO_Cellular_Component_2021",
      "GO Molecular Function" = "GO_Molecular_Function_2021",
      "Hallmark"              = "MSigDB_Hallmark_2020",
      "KEGG"                  = "KEGG_2021_Human",
      "KEGG (Mouse)"          = "KEGG_2019_Mouse",
      "WikiPathways"          = "WikiPathway_2023_Human",
      "WikiPathways (Mouse)"  = "WikiPathways_2019_Mouse",
      "Reactome"              = "Reactome_2022"
    )
    if (!database %in% names(database_mappings)) {
      cat(paste0("The database specified: ", database, " is not supported.\n"))
      return(NULL)
    }
    reverse_mappings <- setNames(names(database_mappings), unname(database_mappings))
    result <- test_ora_mod(se, databases = database_mappings[database], backend = "enrichr",
                           contrasts = TRUE, direction = direction,
                           log2_threshold = log2_threshold, alpha = alpha)
    result$var <- reverse_mappings[result$var]
    return(result)
  } else if (backend == "clusterProfiler") {
    database_mappings <- c(
      "GO Biological Process" = "BP",
      "GO Cellular Component" = "CC",
      "GO Molecular Function" = "MF",
      "Hallmark"              = "Hallmark",
      "KEGG"                  = "KEGG",
      "KEGG (Mouse)"          = "KEGG (Mouse)",
      "WikiPathways"          = "WikiPathways",
      "WikiPathways (Mouse)"  = "WikiPathways (Mouse)"
    )
    if (!database %in% names(database_mappings)) {
      cat(paste0("The database specified: ", database, " is not supported by the clusterProfiler backend.\n"))
      return(NULL)
    }
    reverse_mappings <- setNames(names(database_mappings), unname(database_mappings))
    result <- test_ora_mod(se, databases = database_mappings[database], backend = "clusterProfiler",
                           contrasts = TRUE, direction = direction,
                           log2_threshold = log2_threshold, alpha = alpha)
    result$var <- reverse_mappings[result$var]
    return(result)
  } else {
    cat("Currently, we only support enrichr and clusterProfiler backend.\n")
  }
}

#' Plot over-representation analysis results
#'
#' Creates a dot plot visualization of over-representation analysis results,
#' showing enriched terms ranked by log2 odds ratio.
#'
#' @param or_result A data frame containing ORA results from \code{\link{or_test}}.
#'   Must contain columns: Term, var, contrast, Adjusted.P.value.
#' @param number Numeric value specifying the maximum number of top enriched
#'   terms to display per contrast. Default is 10.
#' @param alpha Numeric value for the significance threshold. Only terms with
#'   p-value below this threshold are shown. Default is 0.05.
#' @param contrasts Character vector specifying which contrasts to plot.
#'   Default is \code{NULL} (plot all contrasts).
#' @param databases Character vector specifying which databases to plot.
#'   Default is \code{NULL} (plot all databases).
#' @param adjust Logical indicating whether to use adjusted p-values (\code{TRUE})
#'   or raw p-values (\code{FALSE}). Default is \code{FALSE}.
#' @param use_whole_proteome Logical indicating whether to use whole proteome
#'   statistics instead of background-corrected statistics. Default is \code{FALSE}.
#' @param nrow Numeric value specifying the number of rows for faceting when
#'   multiple contrasts are plotted. Default is 1.
#' @param term_size Numeric value for term text size (currently unused).
#'   Default is 8.
#'
#' @return A \code{ggplot} object showing a dot plot where:
#'   \itemize{
#'     \item X-axis: Log2 odds ratio
#'     \item Y-axis: Enriched terms
#'     \item Point size: Number of genes in overlap
#'     \item Point color: P-value (red = low/significant, blue = high)
#'   }
#'   If no enrichment is found, returns a plot with "No enrichment found" message.
#'
#' @examples
#' \dontrun{
#' # Perform ORA
#' or_result <- or_test(se_diff, database = "GO Biological Process")
#'
#' # Plot top 15 enriched terms
#' plot_or(or_result, number = 15)
#'
#' # Plot specific contrasts with adjusted p-values
#' plot_or(or_result, contrasts = "Treatment_vs_Control", adjust = TRUE)
#' }
#'
#' @seealso \code{\link{or_test}}
#'
#' @importFrom ggplot2 ggplot aes geom_point facet_wrap scale_color_continuous
#'   labs theme_bw theme element_text annotate
#' @importFrom dplyr filter group_by arrange slice
#' @importFrom readr parse_factor
#'
#' @export
plot_or <- function(or_result, number = 10, alpha = 0.05,
                    contrasts = NULL, databases = NULL,  adjust=F, use_whole_proteome=F,
                    nrow = 1, term_size = 8) {
  assertthat::assert_that(is.data.frame(or_result),
                          is.numeric(number),
                          length(number) == 1,
                          is.numeric(alpha),
                          length(alpha) == 1,
                          is.numeric(term_size),
                          length(term_size) == 1,
                          is.numeric(nrow),
                          length(nrow) == 1)

  # Check or_result object
  if(any(!c("Term", "var",
            "contrast","Adjusted.P.value")
         %in% colnames(or_result))) {
    stop("'", deparse(substitute(or_result)),
         "' does not contain the required columns",
         "\nMake sure that HGNC gene symbols are present",
         "\n in your 'Gene Names' column of Results table",
         call. = FALSE)
  }

  no_enrichment_text <- paste("\n   No enrichment found.\n",
                              "       You can still download enrichment result table. \n")

  if(!is.null(contrasts)) {
    assertthat::assert_that(is.character(contrasts))


    valid_contrasts <- unique(or_result$contrast)

    if(!all(contrasts %in% valid_contrasts)) {
      return(ggplot() +
               annotate("text", x = 4, y = 25, size=8, label = no_enrichment_text) +
               theme_void()
      )
    }
    if(!any(contrasts %in% valid_contrasts)) {
      contrasts <- contrasts[contrasts %in% valid_contrasts]
      message("Not all contrasts found",
              "\n Following contrasts are found: '",
              paste0(contrasts, collapse = "', '"), "'")
    }

    or_result <- filter(or_result, contrast %in% contrasts)
  }
  if(!is.null(databases)) {
    assertthat::assert_that(is.character(databases))

    valid_databases <- unique(or_result$var)

    if(all(!databases %in% valid_databases)) {
      valid_cntrsts_msg <- paste0("Valid databases are: '",
                                  paste0(valid_databases, collapse = "', '"),
                                  "'")
      stop("Not a valid database, please run `plot_gsea()`",
           "with valid databases as argument\n",
           valid_cntrsts_msg,
           call. = FALSE)
    }
    if(any(!databases %in% valid_databases)) {
      databases <- databases[databases %in% valid_databases]
      message("Not all databases found",
              "\nPlotting the following databases: '",
              paste0(databases, collapse = "', '"), "'")
    }

    or_result <- filter(or_result, var %in% databases)
  }

  # Get top enriched gene sets
  if (!use_whole_proteome) {
    if (adjust) {
      terms <- or_result %>%
        dplyr::group_by(contrast, var) %>%
        dplyr::filter(p.adjust_hyper <= alpha) %>%
        dplyr::arrange(p.adjust_hyper) %>%
        dplyr::slice(seq_len(number)) %>%
        .$Term
      subset <- or_result %>%
        dplyr::filter(Term %in% terms) %>%
        dplyr::arrange(var, p.adjust_hyper)
    } else {
      terms <- or_result %>%
        dplyr::group_by(contrast, var) %>%
        dplyr::filter(p_hyper <= alpha) %>%
        dplyr::arrange(p_hyper) %>%
        dplyr::slice(seq_len(number)) %>%
        .$Term
      subset <- or_result %>%
        dplyr::filter(Term %in% terms) %>%
        dplyr::arrange(var, p_hyper)
    }
  } else {
    if (adjust) {
      terms <- or_result %>%
        dplyr::group_by(contrast, var) %>%
        dplyr::filter(Adjusted.P.value <= alpha) %>%
        dplyr::arrange(Adjusted.P.value) %>%
        dplyr::slice(seq_len(number)) %>%
        .$Term
      subset <- or_result %>%
        dplyr::filter(Term %in% terms) %>%
        dplyr::arrange(var, Adjusted.P.value)
    } else {
      terms <- or_result %>%
        dplyr::group_by(contrast, var) %>%
        dplyr::filter(P.value <= alpha) %>%
        dplyr::arrange(P.value) %>%
        dplyr::slice(seq_len(number)) %>%
        .$Term
      subset <- or_result %>%
        dplyr::filter(Term %in% terms) %>%
        dplyr::arrange(var, P.value)
    }
  }
  subset$Term <- readr::parse_factor(subset$Term, levels = unique(subset$Term))
  subset$var <- readr::parse_factor(subset$var, levels = unique(subset$var))

  if (nrow(subset) == 0) {
    return(ggplot() +
             annotate("text", x = 4, y = 25, size=8, label = no_enrichment_text) +
             theme_void()
    )
  } else {
    # Plot top enriched gene sets
    subset$Overlap_ratio <- sapply(subset$Overlap, function(x) eval(parse(text=x)))
    # return(ggplot(subset, aes(y = reorder(Term, Overlap_ratio), x=Overlap_ratio, size=IN, color=Adjusted.P.value)) +
    #   geom_point() +
    #   facet_wrap(~contrast, nrow = nrow) +
    #   scale_color_continuous(low="red", high="blue", name = "Adjusted.P.value",
    #                            guide=guide_colorbar(reverse=TRUE)) +
    #   labs(y = "Term") +
    #   theme_bw() +
    #   theme(legend.position = "top", legend.text = element_text(size = 9))
    # )
    if (!use_whole_proteome) {
      if (adjust){
        return(ggplot(subset, aes(y = reorder(Term, log_odds), x=log_odds, size=IN, color=p.adjust_hyper)) +
                 geom_point() +
                 facet_wrap(~contrast, nrow = nrow) +
                 scale_color_continuous(low="red", high="blue", name = "p.adjust",
                                        guide=guide_colorbar(reverse=T,
                                                             label.theme = element_text(angle = 90),
                                                             label.vjust = 0.5)) +
                 labs(y = "Term", x = "log2 Odds ratio", size = "size") +
                 theme_bw() +
                 theme(legend.position = "top", legend.text = element_text(size = 9)))
      } else {
        return(ggplot(subset, aes(y = reorder(Term, log_odds), x=log_odds, size=IN, color=p_hyper)) +
                 geom_point() +
                 facet_wrap(~contrast, nrow = nrow) +
                 scale_color_continuous(low="red", high="blue", name = "p",
                                        guide=guide_colorbar(reverse=T,
                                                             label.theme = element_text(angle = 90),
                                                             label.vjust = 0.5)) +
                 labs(y = "Term", x = "log2 Odds ratio", size = "size") +
                 theme_bw() +
                 theme(legend.position = "top", legend.text = element_text(size = 9)))
      }
    } else {
      if (adjust){
        return(ggplot(subset, aes(y = reorder(Term, Odds.Ratio), x=log2(Odds.Ratio), size=IN, color=Adjusted.P.value)) +
                 geom_point() +
                 facet_wrap(~contrast, nrow = nrow) +
                 scale_color_continuous(low="red", high="blue", name = "Adjusted.P.value",
                                        guide=guide_colorbar(reverse=T,
                                                             label.theme = element_text(angle = 90),
                                                             label.vjust = 0.5)) +
                 labs(y = "Term", x = "log2 Odds ratio", size = "size") +
                 theme_bw() +
                 theme(legend.position = "top", legend.text = element_text(size = 9)))
      } else {
        return(ggplot(subset, aes(y = reorder(Term, Odds.Ratio), x=log2(Odds.Ratio), size=IN, color=P.value)) +
                 geom_point() +
                 facet_wrap(~contrast, nrow = nrow) +
                 scale_color_continuous(low="red", high="blue", name = "P.value",
                                        guide=guide_colorbar(reverse=T,
                                                             label.theme = element_text(angle = 90),
                                                             label.vjust = 0.5)) +
                 labs(y = "Term", x = "log2 Odds ratio", size = "size") +
                 theme_bw() +
                 theme(legend.position = "top", legend.text = element_text(size = 9)))
      }
    }
  }
}

test_ora_mod <- function(dep,
                         databases, backend = "clusterProfiler",
                         contrasts = TRUE, direction = "UP", log2_threshold = 0.7,
                         alpha = 0.05, adjust_alpha = TRUE) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.character(databases),
                          is.logical(contrasts),
                          length(contrasts) == 1)


  row_data <- rowData(dep, use.names = FALSE)
  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(row_data))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(dep)),
         "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if(length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in '",
         deparse(substitute(dep)),
         "'\nRun test_diff() to obtain the required columns",
         call. = FALSE)
  }



  # Run background list
  message("Background")
  if (metadata(dep)$exp == "TMT" & metadata(dep)$level == "protein") {
    background <- unique(row_data$Gene)
  } else if (metadata(dep)$exp == "TMT" & metadata(dep)$level == "gene") {
    background <- unique(row_data$ID)
  } else if (metadata(dep)$level == "protein") {
      background <- unique(gsub("[.].*", "", row_data$name))
  } else if (metadata(dep)$exp == "LFQ" & metadata(dep)$level == "peptide") {
    background <- unique(row_data$Gene)
  } else if (metadata(dep)$exp == "TMT" & metadata(dep)$level %in% c("peptide", "site")) {
    background <- unique(row_data$Gene)
  } else if (metadata(dep)$exp == "DIA" & metadata(dep)$level == "site") {
    background <- unique(row_data$Gene)
  } else if (metadata(dep)$exp == "DIA" & metadata(dep)$level == "peptide") {
    if ("Gene" %in% colnames(row_data)) {
      background <- unique(row_data$Gene)
    } else {
      background <- unique(row_data$Genes)
    }
  }

  if (backend == "clusterProfiler") {

    if (contrasts) {
      df <- row_data %>%
        as.data.frame() %>%
        # select(name, ends_with("_significant")) %>%
        mutate(name = gsub("[.].*", "", name))

      constrast_columns <- df %>% select(ends_with("_significant")) %>% colnames()
      constrasts <- gsub("_significant", "", constrast_columns)

      df_enrich <- data.frame(ID = character(),
                              Description = character(),
                              GeneRatio = character(),
                              BgRatio = character(),
                              pvalue = character(),
                              p.adjust = character(),
                              qvalue = character(),
                              geneID = character(),
                              Count = character(),
                              stringsAsFactors = FALSE)

      # --- Pre-compute background ID mappings & db-specific setup ONCE ---
      organism_db_map <- ifelse(endsWith(databases, "(Mouse)"), "org.Mm.eg.db", "org.Hs.eg.db")
      bg_entrez    <- NULL
      hallmark     <- NULL
      organism_map <- NULL
      if (startsWith(databases, "KEGG") || startsWith(databases, "WikiPathways")) {
        bg_entrez <- bitr(background, fromType = "SYMBOL",
                          toType = c("ENTREZID"), OrgDb = organism_db_map)$ENTREZID
      }
      if (startsWith(databases, "KEGG")) {
        organism_map <- ifelse(endsWith(databases, "(Mouse)"), "mmu", "hsa")
      } else if (startsWith(databases, "WikiPathways")) {
        organism_map <- ifelse(endsWith(databases, "(Mouse)"), "Mus musculus", "Homo sapiens")
      } else if (databases == "Hallmark") {
        hallmark <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H") %>%
          dplyr::select(gs_name, gene_symbol)
      } else if (!databases %in% c("MF", "BP", "CC")) {
        stop("Not a valid database, please choose from 'KEGG', 'WikiPathways', 'Hallmark', 'MF', 'BP', 'CC'",
             call. = FALSE)
      }

      for(contrast in constrast_columns) {
        df[is.na(df[[contrast]]), contrast] <- FALSE
        significant <- df
        if (direction == "UP") {
          significant <- significant[(significant[gsub("_significant", "_diff", contrast)] > log2_threshold), ]
        } else if (direction == "DOWN") {
          significant <- significant[significant[gsub("_significant", "_diff", contrast)] < -log2_threshold, ]
        }
        if (adjust_alpha) {
          significant <- significant[significant[gsub("_significant", "_p.adj", contrast)] < alpha, ]
        } else {
          significant <- significant[significant[gsub("_significant", "_p.val", contrast)] < alpha, ]
        }

        if (metadata(dep)$exp == "TMT" & metadata(dep)$level == "protein") {
          genes <- unique(significant$Gene)
        } else if (metadata(dep)$exp == "TMT" & metadata(dep)$level == "gene") {
          genes <- unique(significant$ID)
        } else if (metadata(dep)$level == "protein") {
          genes <- significant$name
        } else if (metadata(dep)$exp == "LFQ" & metadata(dep)$level == "peptide") {
          genes <- unique(significant$Gene)
        } else if (metadata(dep)$exp == "TMT" & metadata(dep)$level %in% c("peptide", "site")) {
          genes <- unique(significant$Gene)
        } else if (metadata(dep)$exp == "DIA" & metadata(dep)$level == "peptide") {
          if ("Gene" %in% colnames(row_data)) {
            genes <- unique(significant$Gene)
          } else {
            genes <- unique(significant$Genes)
          }
        } else if (metadata(dep)$exp == "DIA" & metadata(dep)$level == "site") {
          genes <- unique(significant$Gene)
        }

        if (startsWith(databases, "KEGG")) {
          mappings <- bitr(genes, fromType = "SYMBOL",
                           toType = c("ENTREZID"), OrgDb = organism_db_map)
          result <- enrichKEGG(gene = mappings$ENTREZID,
                               universe = bg_entrez,
                               keyType = "ncbi-geneid",
                               organism = organism_map,
                               pvalueCutoff = 1,
                               qvalueCutoff = 1)
          result <- setReadable(result, OrgDb = organism_db_map, keyType = "ENTREZID")
        } else if (startsWith(databases, "WikiPathways")) {
          mappings <- bitr(genes, fromType = "SYMBOL",
                           toType = c("ENTREZID"), OrgDb = organism_db_map)
          result <- enrichWP(gene = mappings$ENTREZID,
                             universe = bg_entrez,
                             organism = organism_map,
                             pvalueCutoff = 1,
                             qvalueCutoff = 1)
          result <- setReadable(result, OrgDb = organism_db_map, keyType = "ENTREZID")
        } else if (databases == "Hallmark") {
          result <- enricher(genes, universe = background,
                             TERM2GENE = hallmark, pvalueCutoff = 1, qvalueCutoff = 1)
        } else if (databases %in% c("MF", "BP", "CC")) {
          result <- enrichGO(gene = genes,
                             universe = background,
                             OrgDb = organism_db_map,
                             ont = databases,
                             keyType = "SYMBOL",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 1,
                             qvalueCutoff  = 1)
        }
        temp <- as.data.frame(result)

        temp$contrast <- gsub("_significant", "", contrast)
        temp$OUT <- length(genes) - temp$Count
        df_enrich <- rbind(df_enrich, temp)
      }
      lookup <- c("Term" = "Description", "P.value" = "pvalue", "Adjusted.P.value" = "p.adjust", "IN" = "Count")
      df_enrich <- dplyr::rename(df_enrich, all_of(lookup))
      df_enrich$var <- databases
      df_enrich$Overlap <- paste0(df_enrich$IN, "/", as.numeric(gsub("/.*", "", df_enrich$BgRatio)))
      # Compute log odds from GeneRatio and BgRatio
      bg_IN  <- as.numeric(gsub("/.*", "", df_enrich$BgRatio))
      bg_OUT <- as.numeric(gsub(".*/", "", df_enrich$BgRatio)) - bg_IN
      df_enrich$Odds.Ratio <- (df_enrich$IN * bg_OUT) / (df_enrich$OUT * bg_IN)
      df_enrich$log_odds   <- log2(df_enrich$Odds.Ratio)
      # Use clusterProfiler's p-values directly (already corrected for universe = background)
      df_enrich$p_hyper        <- df_enrich$P.value
      df_enrich$p.adjust_hyper <- df_enrich$Adjusted.P.value

      return(df_enrich)
    } else {
      # contrasts == FALSE: run enrichment on all significant genes at once
      significant <- row_data %>%
        as.data.frame() %>%
        select(name, significant) %>%
        filter(significant) %>%
        mutate(name = gsub("[.].*", "", name))

      if (metadata(dep)$exp == "TMT" & metadata(dep)$level == "protein") {
        genes <- unique(significant$Gene)
      } else if (metadata(dep)$exp == "TMT" & metadata(dep)$level == "gene") {
        genes <- unique(significant$ID)
      } else if (metadata(dep)$level == "protein") {
        genes <- significant$name
      } else if (metadata(dep)$level == "site" |
                 (metadata(dep)$exp == "TMT" & metadata(dep)$level == "peptide")) {
        genes <- unique(significant$Gene)
      } else { # DIA-peptide
        if ("Gene" %in% colnames(row_data)) {
          genes <- unique(significant$Gene)
        } else {
          genes <- unique(significant$Genes)
        }
      }

      organism_db_map <- ifelse(endsWith(databases, "(Mouse)"), "org.Mm.eg.db", "org.Hs.eg.db")
      bg_entrez    <- NULL
      hallmark     <- NULL
      organism_map <- NULL
      if (startsWith(databases, "KEGG") || startsWith(databases, "WikiPathways")) {
        bg_entrez <- bitr(background, fromType = "SYMBOL",
                          toType = c("ENTREZID"), OrgDb = organism_db_map)$ENTREZID
      }

      if (startsWith(databases, "KEGG")) {
        organism_map <- ifelse(endsWith(databases, "(Mouse)"), "mmu", "hsa")
        mappings <- bitr(genes, fromType = "SYMBOL",
                         toType = c("ENTREZID"), OrgDb = organism_db_map)
        result <- enrichKEGG(gene = mappings$ENTREZID,
                             universe = bg_entrez,
                             keyType = "ncbi-geneid",
                             organism = organism_map,
                             pvalueCutoff = 1,
                             qvalueCutoff = 1)
        result <- setReadable(result, OrgDb = organism_db_map, keyType = "ENTREZID")
      } else if (startsWith(databases, "WikiPathways")) {
        organism_map <- ifelse(endsWith(databases, "(Mouse)"), "Mus musculus", "Homo sapiens")
        mappings <- bitr(genes, fromType = "SYMBOL",
                         toType = c("ENTREZID"), OrgDb = organism_db_map)
        result <- enrichWP(gene = mappings$ENTREZID,
                           universe = bg_entrez,
                           organism = organism_map,
                           pvalueCutoff = 1,
                           qvalueCutoff = 1)
        result <- setReadable(result, OrgDb = organism_db_map, keyType = "ENTREZID")
      } else if (databases == "Hallmark") {
        hallmark <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H") %>%
          dplyr::select(gs_name, gene_symbol)
        result <- enricher(genes, universe = background,
                           TERM2GENE = hallmark, pvalueCutoff = 1, qvalueCutoff = 1)
      } else if (databases %in% c("MF", "BP", "CC")) {
        result <- enrichGO(gene = genes,
                           universe = background,
                           OrgDb = organism_db_map,
                           ont = databases,
                           keyType = "SYMBOL",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 1,
                           qvalueCutoff  = 1)
      } else {
        stop("Not a valid database, please choose from 'KEGG', 'WikiPathways', 'Hallmark', 'MF', 'BP', 'CC'",
             call. = FALSE)
      }

      df_enrich <- as.data.frame(result)
      df_enrich$contrast <- "significant"
      df_enrich$OUT <- length(genes) - df_enrich$Count

      lookup <- c("Term" = "Description", "P.value" = "pvalue", "Adjusted.P.value" = "p.adjust", "IN" = "Count")
      df_enrich <- dplyr::rename(df_enrich, all_of(lookup))
      df_enrich$var <- databases
      df_enrich$Overlap <- paste0(df_enrich$IN, "/", as.numeric(gsub("/.*", "", df_enrich$BgRatio)))
      bg_IN  <- as.numeric(gsub("/.*", "", df_enrich$BgRatio))
      bg_OUT <- as.numeric(gsub(".*/", "", df_enrich$BgRatio)) - bg_IN
      df_enrich$Odds.Ratio <- (df_enrich$IN * bg_OUT) / (df_enrich$OUT * bg_IN)
      df_enrich$log_odds   <- log2(df_enrich$Odds.Ratio)
      df_enrich$p_hyper        <- df_enrich$P.value
      df_enrich$p.adjust_hyper <- df_enrich$Adjusted.P.value

      return(df_enrich)
    }


  } else if (backend == "enrichr") {
    background_enriched <- enrichr_mod(background, databases)
    df_background <- NULL
    for(database in databases) {
      temp <- background_enriched[database][[1]] %>%
        mutate(var = database)
      df_background <- rbind(df_background, temp)
    }
    df_background$contrast <- "background"
    df_background$n <- length(background)

    OUT <- df_background %>%
      mutate(bg_IN = as.numeric(gsub("/.*", "", Overlap)),
             bg_OUT = n - bg_IN) %>%
      select(Term, bg_IN, bg_OUT)

    if(contrasts) {
      # Get gene symbols
      df <- row_data %>%
        as.data.frame() %>%
        # select(name, ends_with("_significant")) %>%
        mutate(name = gsub("[.].*", "", name))

      constrast_columns <- df %>% select(ends_with("_significant")) %>% colnames()
      constrasts <- gsub("_significant", "", constrast_columns)

      # Run enrichR for every contrast
      df_enrich <- NULL
      for(contrast in constrast_columns) {
        message(gsub("_significant", "", contrast))
        # contrast column might have NA
        df[is.na(df[[contrast]]), contrast] <- FALSE
        significant <- df
        if (direction == "UP") {
          significant <- significant[(significant[gsub("_significant", "_diff", contrast)] > log2_threshold), ]
        } else if (direction == "DOWN") {
          significant <- significant[significant[gsub("_significant", "_diff", contrast)] < -log2_threshold, ]
        }
        if (adjust_alpha) {
          significant <- significant[significant[gsub("_significant", "_p.adj", contrast)] < alpha, ]
        } else {
          significant <- significant[significant[gsub("_significant", "_p.val", contrast)] < alpha, ]
        }

        if (metadata(dep)$exp == "TMT" & metadata(dep)$level == "protein") {
          genes <- unique(significant$Gene)
        } else if (metadata(dep)$exp == "TMT" & metadata(dep)$level == "gene") {
          genes <- unique(significant$ID)
        } else if (metadata(dep)$level == "protein") {
          genes <- significant$name
        } else if (metadata(dep)$exp == "LFQ" & metadata(dep)$level == "peptide") {
          genes <- unique(significant$Gene)
        } else if (metadata(dep)$exp == "TMT" & metadata(dep)$level %in% c("peptide", "site")) {
          genes <- unique(significant$Gene)
        } else if (metadata(dep)$exp == "DIA" & metadata(dep)$level == "peptide") {
          if ("Gene" %in% colnames(row_data)) {
            genes <- unique(significant$Gene)
          } else {
            genes <- unique(significant$Genes)
          }
        } else if (metadata(dep)$exp == "DIA" & metadata(dep)$level == "site") {
          genes <- unique(significant$Gene)
        }

        message(paste0(length(genes), " genes are submitted"))
        if (length(genes) != 0) {
          enriched <- enrichr_mod(genes, databases)
          # Tidy output
          contrast_enrich <- NULL
          for(database in databases) {
            temp <- enriched[database][[1]] %>%
              mutate(var = database)
            contrast_enrich <- rbind(contrast_enrich, temp)
          }
          if (nrow(contrast_enrich) != 0) { # has enrichment
            contrast_enrich$contrast <- contrast
            contrast_enrich$n <- length(genes)
            # Background correction
            cat("Background correction... ")
            contrast_enrich <- contrast_enrich %>%
              mutate(IN = as.numeric(gsub("/.*", "", Overlap)),
                     OUT = n - IN) %>%
              select(-n) %>%
              left_join(OUT, by = "Term") %>%
              mutate(log_odds = log2((IN * bg_OUT) / (OUT * bg_IN)))
            cat("Done.")
          }
          df_enrich <- rbind(df_enrich, contrast_enrich) %>%
            mutate(contrast = gsub("_significant", "", contrast))
        } else {
          cat("No significant genes for enrichment analysis")
        }
      }
    } else {
      # Get gene symbols
      significant <- row_data %>%
        as.data.frame() %>%
        select(name, significant) %>%
        filter(significant) %>%
        mutate(name = gsub("[.].*", "", name))

      # Run enrichR
      if (metadata(dep)$exp == "TMT" & metadata(dep)$level == "protein") {
        genes <- unique(significant$Gene)
      } else if (metadata(dep)$exp == "TMT" & metadata(dep)$level == "gene") {
        genes <- unique(significant$ID)
      } else if (metadata(dep)$level == "protein") {
        genes <- significant$name
      } else if (metadata(dep)$level == "site" |
                 (metadata(dep)$exp == "TMT" & metadata(dep)$level == "peptide")) {
        genes <- unique(significant$Gene)
      } else { # DIA-peptide
        if ("Gene" %in% colnames(row_data)) {
          genes <- unique(significant$Gene)
        } else {
          genes <- unique(significant$Genes)
        }
      }

      enriched <- enrichr_mod(genes, databases)

      # Tidy output
      df_enrich <- NULL
      for(database in databases) {
        temp <- enriched[database][[1]] %>%
          mutate(var = database)
        df_enrich <- rbind(df_enrich, temp)
      }
      df_enrich$contrast <- "significant"
      df_enrich$n <- length(genes)

      # Background correction
      cat("Background correction... ")
      df_enrich <- df_enrich %>%
        mutate(IN = as.numeric(gsub("/.*", "", Overlap)),
               OUT = n - IN) %>%
        select(-n) %>%
        left_join(OUT, by = "Term") %>%
        mutate(log_odds = log2((IN * bg_OUT) / (OUT * bg_IN)))
      cat("Done.")
    }

    if (nrow(df_enrich) != 0) {
      df_enrich$p_hyper = phyper(q=(df_enrich$IN-1), m = df_enrich$bg_IN, n = df_enrich$bg_OUT, k = (df_enrich$IN+df_enrich$OUT),
                                 lower.tail = F )
      df_enrich$p.adjust_hyper = p.adjust(df_enrich$p_hyper, method = "BH")
    }
    return(df_enrich)

  } else {
    stop("Not a valid backend, please choose 'enrichr' or 'clusterProfiler'",
         call. = FALSE)
  }
}

enrichr_mod <- function(genes, databases = NULL) {
  if (length(genes) != 0) {
    if (all(startsWith(genes, "ENSG"))) {
      genes_map <- ensembldb::select(EnsDb.Hsapiens.v86,
                                     keys = genes, keytype = "GENEID", columns = c("SYMBOL", "GENEID"))
      genes <- genes_map$SYMBOL
    }
    httr::set_config(httr::config(ssl_verifypeer = 0L))
    cat("Uploading data to Enrichr... ")
    if (is.vector(genes) & !all(genes == "") & length(genes) != 0) {
      gene_str <- paste(genes, collapse = "\n")
    } else if (is.data.frame(genes)) {
      gene_str <- paste(paste(genes[,1], genes[,2], sep = ","), collapse = "\n")
    } else {
      warning("genes must be a non-empty vector of gene names or a dataframe with genes and score.")
      return(NULL)
    }
    temp <- httr::POST(url = "https://maayanlab.cloud/Enrichr/addList",
                      body = list(list = gene_str, description = ""),
                      encode = "multipart", httr::timeout(30))
    Sys.sleep(1)
    user_list_id <- httr::content(temp, as = "parsed", type = "application/json")$userListId
    cat("Done.\n")
    dbs <- as.list(databases)
    dfSAF <- options()$stringsAsFactors
    options(stringsAsFactors = FALSE)
    result <- lapply(dbs, function(x) {
      cat("  Querying ", x, "... ", sep = "")
      r <- httr::GET(url = "https://maayanlab.cloud/Enrichr/enrich",
                     query = list(userListId = user_list_id, backgroundType = x), httr::timeout(30))
      data <- httr::content(r, as = "parsed", type = "application/json")[[x]]
      if (is.null(data) || length(data) == 0) {
        cat("Done.\n")
        return(data.frame(Term = character(), Overlap = character(),
                          P.value = double(), Adjusted.P.value = double(),
                          Old.P.value = double(), Old.Adjusted.P.value = double(),
                          Odds.Ratio = double(), Combined.Score = double(),
                          Genes = character()))
      }
      # JSON entry: [rank, term, p_value, Odds ratio, combined_score, [genes], adj_p_value, old_p_value, old_adj_p_value]
      Sys.sleep(1)
      df <- do.call(rbind, lapply(data, function(entry) {
        data.frame(
          Term                 = entry[[2]],
          Overlap              = as.character(length(entry[[6]])),
          P.value              = entry[[3]],
          Adjusted.P.value     = entry[[7]],
          Old.P.value          = if (length(entry) >= 8) entry[[8]] else NA_real_,
          Old.Adjusted.P.value = if (length(entry) >= 9) entry[[9]] else NA_real_,
          Odds.Ratio           = entry[[4]],
          Combined.Score       = entry[[5]],
          Genes                = paste(entry[[6]], collapse = ";"),
          stringsAsFactors     = FALSE
        )
      }))
      cat("Done.\n")
      return(df)
    })
    options(stringsAsFactors = dfSAF)
    cat("Parsing results... ")
    names(result) <- dbs
    cat("Done.\n")
  } else {
    result <- data.frame(Term = character(),
                         Overlap = character(),
                         P.value = double(),
                         Adjusted.P.value = double(),
                         Old.P.value = double(),
                         Old.Adjusted.P.value = double(),
                         Odds.Ratio = double(),
                         Combined.Score = double(),
                         Genes = character())
  }
  return(result)
}

#' Prepare data for PTM-SEA analysis
#'
#' Exports phosphoproteomics data in GCT format for PTM Signature Enrichment
#' Analysis (PTM-SEA) using ssGSEA2.
#'
#' @param se A \code{SummarizedExperiment} object containing phosphoproteomics data.
#' @param score_col Character string specifying the column name in rowData
#'   containing the score values to export (e.g. a log2 fold change column like
#'   \code{"Treatment_vs_Control_diff"}).
#' @param id_col Character string specifying the column name in rowData containing
#'   peptide/site identifiers. Default is \code{"SequenceWindow"}.
#' @param outfile Character string specifying the output file path for the GCT file.
#'
#' @return Writes a GCT format file to the specified path. Returns \code{NULL}
#'   invisibly.
#'
#' @details
#' Peptide identifiers are formed by converting \code{id_col} values to uppercase
#' and appending \code{"-p"}.  Rows where \code{score_col} is \code{NA} are
#' dropped.  Duplicate identifiers are made unique with \code{make.unique()}.
#'
#' @examples
#' \dontrun{
#' prepare_PTMSEA(se_phospho, score_col = "Treatment_vs_Control_diff",
#'                outfile = "ptmsea_input.gct")
#' }
#'
#' @seealso \code{\link{visualize_PTMSEA}}
#'
#' @importFrom cmapR GCT
#'
#' @export
prepare_PTMSEA <- function(se, score_col, id_col = "SequenceWindow", outfile) {
  temp <- data.frame(as.data.frame(rowData(se))[id_col],
                     as.data.frame(rowData(se))[score_col])
  colnames(temp) <- c("peptide", score_col)
  temp <- temp[!is.na(temp[[score_col]]), ]
  gct <- new("GCT", mat = as.matrix(temp[, score_col, drop = FALSE]),
             rid = paste0(make.unique(toupper(temp$peptide)), "-p"))
  write_gct(gct, outfile, appenddim = FALSE)
  invisible(NULL)
}

write_gct <- function(ds, ofile, precision=4, appenddim=TRUE, ver=3) {
  if (!methods::is(ds, "GCT")) {
    stop("ds must be a GCT object")
  }
  # make sure it's valid
  methods::validObject(ds)

  # extract the components
  m <- mat(ds)
  rdesc <- meta(ds)
  cdesc <- meta(ds, dimension="column")
  rid <- ids(ds)
  cid <- ids(ds, dimension="column")

  # append the dimensions of the data set, if desired
  if (appenddim) ofile <- append_dim(ofile, m, extension="gct")

  precision <- floor(precision)
  cat("Saving file to ", ofile, "\n")
  nr <- nrow(m)
  nc <- ncol(m)
  cat(sprintf("Dimensions of matrix: [%dx%d]\n", nr, nc))
  cat(sprintf("Setting precision to %d\n", precision))
  # open file and write
  if (ver == 3) {
    # remove the 'id' columns
    cdesc$id <- NULL
    rdesc$id <- NULL
    # get the counts of meta data fields
    nrdesc <- ncol(rdesc)
    ncdesc <- ncol(cdesc)
    colkeys <- names(cdesc)
    # append header
    cat(sprintf("#1.%d\n%d\t%d\t%d\t%d", ver, nr, nc, nrdesc, ncdesc),
        file=ofile,sep='\n')
    # line 3: sample row desc keys and sample names
    cat(paste(c("id", names(rdesc), cid), collapse="\t"),
        file=ofile, sep="\n", append=TRUE)
    # line 4 + ncdesc: sample desc
    filler <- 'na'
    if (ncdesc > 0) {
      for (ii in seq_len(ncdesc)) {
        if (is.numeric(cdesc[, ii])) {
          cat(paste(c(colkeys[ii], rep(filler, nrdesc),
                      round(cdesc[, ii], precision)),
                    collapse="\t"),
              file=ofile, sep="\n", append=TRUE)
        } else {
          cat(paste(c(colkeys[ii], rep(filler, nrdesc),
                      cdesc[, ii]),
                    collapse="\t"),
              file=ofile, sep="\n", append=TRUE)
        }
      }
    }
    file_o <-file(ofile, "a")
    temp <- cbind(rid, rdesc[1:nr, ], round(m, precision))
    temp <- temp %>% unite("merge", 1:ncol(temp), sep="\t")
    writeLines(temp$merge, file_o)
    close(file_o)
  } else {
    # assume ver 1.2 and below, ignore descriptors
    # append header
    cat(sprintf("#1.%d\n%d\t%d", ver, nr, nc),
        file=ofile, sep="\n")
    # line 3: sample row desc keys and sample names
    cat(paste(c("id", "Description", cid), collapse="\t"),
        file=ofile, sep="\n", append=TRUE)
    for (ii in seq_len(nr)) {
      # print rows
      cat(paste(c(rid[ii],
                  rdesc[ii, 2],
                  round(m[ii, ], precision)), collapse="\t"),
          sep="\n", file=ofile, append=TRUE)
    }
  }
  cat("Saved.\n")
}

#' Visualize PTM-SEA results
#'
#' Creates a dot plot visualization of PTM Signature Enrichment Analysis results
#' from a GCT output file.
#'
#' @param gct_file Character string specifying the path to the PTM-SEA output
#'   GCT file (typically generated by ssGSEA2).
#' @param score_col Character string specifying the score column to visualize.
#' @param selected_collections Character vector of PTMsigDB collection prefixes
#'   to include. Default is \code{c('PERT', 'PATH', 'KINASE', 'DISEASE')}.
#' @param selected_concepts Character vector of specific concept/signature IDs
#'   to display. Default is \code{NULL} (show top/bottom enriched concepts).
#' @param num_concepts Numeric value specifying the number of top and bottom
#'   enriched concepts to display when \code{selected_concepts} is \code{NULL}.
#'   Default is 5 (showing up to 10 total).
#' @param direction Character string to filter by enrichment direction.
#'   Options: \code{"Both"} (default), \code{"Up"} (positive scores only),
#'   \code{"Down"} (negative scores only).
#' @param fdr_pvalue_cutoff Numeric value for the FDR-adjusted p-value threshold.
#'   Default is 0.05.
#' @param score_cutoff Numeric value for the minimum absolute enrichment score.
#'   Default is 1.
#'
#' @return A \code{ggplot} object showing a dot plot where:
#'   \itemize{
#'     \item X-axis: Enrichment score
#'     \item Y-axis: PTM signature names
#'     \item Point size: Signature set overlap percentage
#'     \item Point color: FDR-adjusted p-value
#'   }
#'   Returns \code{NULL} if the specified column is not found or no signatures
#'   pass the filters.
#'
#' @examples
#' \dontrun{
#' # Visualize all significant PTM-SEA results
#' p <- visualize_PTMSEA("ptmsea_output.gct", score_col = "Treatment_vs_Control")
#'
#' # Show only kinase signatures
#' p <- visualize_PTMSEA("ptmsea_output.gct", score_col = "Treatment_vs_Control",
#'                       selected_collections = "KINASE")
#'
#' # Show specific signatures
#' p <- visualize_PTMSEA("ptmsea_output.gct", score_col = "Treatment_vs_Control",
#'                       selected_concepts = c("KINASE-PSP_AKT1", "KINASE-PSP_AKT2"))
#' }
#'
#' @seealso \code{\link{prepare_PTMSEA}}
#'
#' @importFrom cmapR parse_gctx mat meta
#' @importFrom ggplot2 ggplot aes geom_point scale_size_continuous
#'   scale_color_continuous theme_bw
#'
#' @export
visualize_PTMSEA <- function(gct_file, score_col,
                             selected_collections = c("PERT", "PATH", "KINASE", "DISEASE"),
                             selected_concepts = NULL,
                             num_concepts = 5,
                             direction = "Both",
                             fdr_pvalue_cutoff = 0.05,
                             score_cutoff = 1) {
  gct <- parse_gctx(gct_file)
  score_cols <- colnames(mat(gct))
  if (!score_col %in% score_cols) {
    return(NULL)
  }
  data <- cbind(meta(gct, dimension = "row"), mat(gct))
  data <- data[grepl(paste0("^(", paste(selected_collections, collapse = "|"), ")"), data$id), ]
  if (is.null(selected_concepts)) {
    data <- data[order(data[[score_col]]), ]
    data <- data[c(seq_len(num_concepts), (nrow(data) - num_concepts + 1):nrow(data)), ]
  } else {
    data <- data[data$id %in% selected_concepts, ]
    data <- data[order(data[[score_col]]), ]
  }
  if (direction == "Up") {
    data <- data[data[[score_col]] > 0, ]
  } else if (direction == "Down") {
    data <- data[data[[score_col]] < 0, ]
  }
  fdr_col <- paste0("fdr.pvalue.", score_col)
  data <- data[abs(data[[score_col]]) >= score_cutoff & data[[fdr_col]] <= fdr_pvalue_cutoff, ]
  if (nrow(data) == 0) {
    return(NULL)
  }
  data$id <- factor(data$id, levels = data$id)
  overlap_col <- paste0("Signature.set.overlap.percent.", score_col)
  p <- ggplot(data, aes(x = .data[[score_col]], y = .data[["id"]],
                        color = .data[[fdr_col]],
                        size = .data[[overlap_col]])) +
    geom_point() +
    scale_size_continuous(range = c(0.5, 11), name = "Overlap percentage") +
    scale_color_continuous(low = "red", high = "blue", name = "FDR p-value",
                           guide = guide_colorbar(reverse = TRUE, label.vjust = 0.5)) +
    theme_bw()
  return(p)
}

#' Prepare data for kinome enrichment analysis
#'
#' Exports phosphoproteomics data in a format compatible with the PhosphoSitePlus
#' Kinase Library enrichment analysis tool.
#'
#' @param se A \code{SummarizedExperiment} object containing phosphoproteomics data.
#'   Must have \code{exp_type = "phospho"} in metadata and a "SequenceWindow"
#'   column in rowData.
#' @param col Character string specifying the column name in rowData containing
#'   log2 fold change values.
#' @param outfile Character string specifying the output file path.
#' @param format Character string specifying the phosphosite format in the output.
#'   Options are:
#'   \itemize{
#'     \item \code{"asterisk"}: Mark phosphosites with asterisks (e.g., "s*", "t*", "y*")
#'     \item \code{"central"}: Convert to uppercase without markers
#'   }
#' @param p_col Character string specifying the column name containing p-values.
#'   Default is \code{NULL} (p-values not included in output).
#'
#' @return Writes a tab-separated file to the specified path containing:
#'   \itemize{
#'     \item Column 1: Site sequence
#'     \item Column 2: Log2 fold change
#'     \item Column 3: P-value (if \code{p_col} specified)
#'   }
#'
#' @details
#' The output file can be uploaded to the PhosphoSitePlus Kinase Library
#' enrichment analysis tool at https://kinase-library.phosphosite.org/ea
#'
#' @examples
#' \dontrun{
#' # Prepare kinome input with asterisk format
#' prepare_kinome(se_phospho, col = "Treatment_vs_Control_diff",
#'                outfile = "kinome_input.tsv", format = "asterisk")
#'
#' # Include p-values
#' prepare_kinome(se_phospho, col = "Treatment_vs_Control_diff",
#'                p_col = "Treatment_vs_Control_p.val",
#'                outfile = "kinome_input.tsv")
#' }
#'
#' @seealso \code{\link{visualize_kinome}}
#'
#' @export
prepare_kinome <- function(se, col, outfile, format="asterisk", p_col=NULL) {
  # Generate input for kinome analysis: https://kinase-library.phosphosite.org/ea
  # Column 1: site sequence
  # Column 2: log fold change
  # Column 3: p-value (optional)
  if (metadata(se)$exp_type == "phospho") {
    if (is.null(p_col)) {
      temp <- data.frame(as.data.frame(rowData(se))["SequenceWindow"],
                         as.data.frame(rowData(se))[col])
      colnames(temp) <- c("peptide", "log2fc")
    } else {
      temp <- data.frame(as.data.frame(rowData(se))["SequenceWindow"],
                         as.data.frame(rowData(se))[col],
                         as.data.frame(rowData(se))[p_col])
      colnames(temp) <- c("peptide", "log2fc", "p.value")
    }
    temp <- temp[!is.na(temp$log2fc),]
    if (format == "asterisk") {
      temp$peptide <- gsub("s", "s*", temp$peptide)
      temp$peptide <- gsub("t", "t*", temp$peptide)
      temp$peptide <- gsub("y", "y*", temp$peptide)
    } else if (format == "central"){
      temp$peptide <- toupper(temp$peptide)
    } else {
      return(NULL)
    }
    write.table(temp, outfile, quote = F, row.names = F, col.names = T, sep="\t")
  } else {
    cat("Kinome analysis requires phosphorylation dataset. \n")
  }
}

#' Visualize kinome enrichment analysis results
#'
#' Creates a volcano-style plot of kinase enrichment results from the
#' PhosphoSitePlus Kinase Library analysis.
#'
#' @param tsv_file Character string specifying the path to the kinase enrichment
#'   results file (downloaded from PhosphoSitePlus Kinase Library).
#' @param labels Character vector of kinase names to highlight in red.
#'   Default is \code{NULL} (highlight based on thresholds).
#' @param log2fc Numeric value for the log2 enrichment score threshold for
#'   labeling kinases. Default is 1.
#' @param pval Numeric value for the adjusted p-value threshold for labeling
#'   kinases. Default is 0.05.
#' @param legacy Logical indicating whether to use legacy file format parsing.
#'   Default is \code{FALSE}.
#'
#' @return A \code{ggplot} object showing a scatter plot where:
#'   \itemize{
#'     \item X-axis: Log2 enrichment score
#'     \item Y-axis: -log10 adjusted p-value
#'     \item Labels: Kinase names (either based on thresholds or specified labels)
#'   }
#'
#' @examples
#' \dontrun{
#' # Visualize kinome results with default thresholds
#' p <- visualize_kinome("kinase_enrichment_results.tsv")
#'
#' # Highlight specific kinases
#' p <- visualize_kinome("kinase_enrichment_results.tsv",
#'                       labels = c("AKT1", "MAPK1", "CDK1"))
#'
#' # Use stricter thresholds
#' p <- visualize_kinome("kinase_enrichment_results.tsv",
#'                       log2fc = 2, pval = 0.01)
#' }
#'
#' @seealso \code{\link{prepare_kinome}}
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_text_repel theme_bw theme
#'   element_blank element_line labs
#' @importFrom data.table fread
#'
#' @export
visualize_kinome <- function(tsv_file, labels=NULL, log2fc=1, pval=0.05, legacy=F) {
  if (legacy) {
    data <- read.csv(tsv_file, sep = "\t", stringsAsFactors = F)
  } else {
    data <- fread(tsv_file, sep="\t", stringsAsFactors = F, data.table = F)
    colnames(data)[colnames(data) %in% c("Gene Name")] <- "kinase"
  }
  if (is.null(labels)){
    p <- ggplot(data, aes(x=dominant_enrichment_value_log2, y=dominant_adjusted_p_value_log10_abs, label=kinase)) +
      geom_point() +
      geom_text_repel(data=subset(data, abs(dominant_enrichment_value_log2) > log2fc & (dominant_adjusted_p_value_log10_abs) > -log(pval))) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      labs(x="Enrichment Score", y="-log10(adjusted p-value)")
  } else {
    p <- ggplot(data, aes(x=dominant_enrichment_value_log2, y=dominant_adjusted_p_value_log10_abs, label=kinase)) +
      geom_point(data=subset(data, !kinase %in% labels), color="black") +
      geom_point(data=subset(data, kinase %in% labels), color="red") +
      geom_text_repel(data=subset(data, kinase %in% labels), color="red") +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      labs(x="Enrichment Score", y="-log10(adjusted p-value)")
  }
  return(p)
}

#' Kinase activity Z-score inference (RoKAI)
#'
#' Matches DE phosphosite data to a curated kinase-substrate library via
#' flanking sequences, then computes a Z-score for each kinase across all
#' contrasts (RoKAI formulation: \code{Z = mean(FC_substrates) * sqrt(m) / sd(all FC)}).
#'
#' @param dep A \code{SummarizedExperiment} after \code{test_limma()}.
#'   Must have a \code{SequenceWindow} column in \code{rowData}.
#' @param ks_library data.frame with columns: \code{source} (kinase name),
#'   \code{target} (substrate), \code{sequence} (15-char flanking window),
#'   \code{mor} (mode of regulation, +1/-1).
#' @param min_targets Minimum number of matched substrates required per kinase.
#' @return data.frame with columns: \code{kinase}, \code{n_substrates},
#'   \code{score} (Z), \code{p_value}, \code{adj_p_value}, \code{contrast}.
#' @seealso \code{\link{build_ks_data}}, \code{\link{prepare_PTMSEA}}
#' @importFrom SummarizedExperiment rowData
#' @export
run_kinase_zscore <- function(dep, ks_library, min_targets = 3L) {
  rd <- as.data.frame(SummarizedExperiment::rowData(dep))
  if (!"SequenceWindow" %in% colnames(rd))
    stop("SequenceWindow column missing — kinase activity requires site-level data.")

  all_contrasts <- gsub("_significant$", "",
                        grep("_significant$", colnames(rd), value = TRUE))
  if (length(all_contrasts) == 0)
    stop("No contrast columns found. Run test_limma() first.")

  sw_upper  <- toupper(trimws(rd$SequenceWindow))
  stat_list <- lapply(all_contrasts, function(ct) {
    dc <- paste0(ct, "_diff")
    if (!dc %in% colnames(rd)) return(NULL)
    rd[[dc]]
  })
  names(stat_list) <- all_contrasts
  stat_list <- Filter(Negate(is.null), stat_list)
  mat <- do.call(cbind, stat_list)
  rownames(mat) <- sw_upper

  if (anyDuplicated(rownames(mat))) {
    row_mean_abs <- rowMeans(abs(mat), na.rm = TRUE)
    mat <- mat[order(-row_mean_abs), , drop = FALSE]
    mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
  }
  valid <- !is.na(rownames(mat)) & nchar(rownames(mat)) > 0 &
    rowSums(!is.na(mat)) > 0
  mat <- mat[valid, , drop = FALSE]
  mat[is.na(mat)] <- 0

  lib_seq <- toupper(trimws(ks_library$sequence))
  matched <- lib_seq %in% rownames(mat) & !is.na(lib_seq)
  if (sum(matched) == 0)
    stop("No kinase-substrate library entries match the input data.")

  network <- data.frame(
    source = ks_library$source[matched],
    target = lib_seq[matched],
    mor    = as.numeric(ks_library$mor[matched]),
    stringsAsFactors = FALSE
  )
  network <- network[!duplicated(paste0(network$source, "|", network$target)), ]
  message(sprintf("[KA] %d KS pairs (%d kinases) matched to %d sites",
                  nrow(network), length(unique(network$source)), nrow(mat)))

  min_targets    <- as.integer(min_targets)
  kinase_targets <- split(network$target, network$source)
  kinase_targets <- kinase_targets[lengths(kinase_targets) >= min_targets]
  if (length(kinase_targets) == 0)
    stop("No kinases with >= ", min_targets, " matched substrates.")

  contrast_names <- colnames(mat)
  res_list <- list()
  for (j in seq_along(contrast_names)) {
    fc_vec <- mat[, j]
    delta  <- stats::sd(fc_vec, na.rm = TRUE)
    if (is.na(delta) || delta == 0) delta <- 1
    for (kinase in names(kinase_targets)) {
      sub_fc <- fc_vec[kinase_targets[[kinase]]]
      sub_fc <- sub_fc[!is.na(sub_fc)]
      m      <- length(sub_fc)
      if (m < min_targets) next
      z    <- mean(sub_fc) * sqrt(m) / delta
      pval <- stats::pnorm(-abs(z))
      res_list[[length(res_list) + 1]] <- data.frame(
        kinase       = kinase,
        n_substrates = m,
        score        = z,
        p_value      = pval,
        contrast     = contrast_names[j],
        stringsAsFactors = FALSE
      )
    }
  }
  res <- do.call(rbind, res_list)
  if (is.null(res) || nrow(res) == 0)
    stop("Kinase activity inference produced no results.")

  res <- do.call(rbind, lapply(split(res, res$contrast), function(df) {
    df$adj_p_value <- stats::p.adjust(df$p_value, method = "BH")
    df
  }))
  rownames(res) <- NULL
  res
}
