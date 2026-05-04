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
#'     \item \code{"KEGG"}: KEGG pathways
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
#' @export
or_test <- function(se, database="GO Biological Process", backend="enrichr", direction="UP", log2_threshold=0.7, alpha=0.05) {
  if (backend == "enrichr") {
    database_mappings <- c(
      "GO Biological Process" = "GO_Biological_Process_2021",
      "GO Cellular Component" = "GO_Cellular_Component_2021",
      "GO Molecular Function" = "GO_Molecular_Function_2021",
      "Hallmark" = "MSigDB_Hallmark_2020",
      "KEGG" = "KEGG_2021_Human",
      "Reactome" = "Reactome_2022"
    )
    reverse_database_mappings <- c(
      "GO_Biological_Process_2021" = "GO Biological Process",
      "GO_Cellular_Component_2021" = "GO Cellular Component",
      "GO_Molecular_Function_2021" = "GO Molecular Function",
      "MSigDB_Hallmark_2020" = "Hallmark",
      "KEGG_2021_Human" = "KEGG",
      "Reactome_2022" = "Reactome"
    )
    result <- test_ora_mod(se, databases = database_mappings[database], contrasts = T,
                           direction = direction, log2_threshold = log2_threshold, alpha = alpha)
    result$var <- reverse_database_mappings[result$var]
    return(result)
  } else if (backend == "clusterProfiler") {
    row_data <- rowData(se, use.names = FALSE)
    df <- row_data %>%
      as.data.frame() %>%
      # select(name, ends_with("_significant")) %>%
      mutate(name = gsub("[.].*", "", name))

    constrast_columns <- df %>% select(ends_with("_significant")) %>% colnames()
    constrasts <- gsub("_significant", "", constrast_columns)
    if (database %in% c("KEGG", "GO Biological Process", "GO Cellular Component", "GO Molecular Function", "Hallmark")) {
      combined_df <- data.frame(ID=character(),
                       Description=character(),
                       GeneRatio=character(),
                       BgRatio=character(),
                       pvalue=character(),
                       p.adjust=character(),
                       qvalue=character(),
                       geneID=character(),
                       Count=character(),
                       stringsAsFactors=FALSE)
      for(contrast in constrast_columns) {
        df[is.na(df[[contrast]]),contrast] <- F
        significant <- df
        significant <- significant[!is.na(significant[gsub("_significant", "_diff", contrast)]),]
        if (direction == "UP"){
          significant <- significant[(significant[gsub("_significant", "_diff", contrast)] > log2_threshold),]
        } else if (direction == "DOWN") {
          significant <- significant[significant[gsub("_significant", "_diff", contrast)] < -log2_threshold,]
        }

        if (metadata(dep)$level == "protein" & metadata(dep)$exp == "TMT") {
          genes <- unique(significant$Gene)
        } else if (metadata(dep)$level == "protein" | metadata(dep)$level == "gene") {
          genes <- significant$name
        } else if (metadata(dep)$level %in% c("peptide", "site", "glycan")) {
          genes <- unique(significant$Gene)
        }

        background <- df$name
        significant <- significant[significant[gsub("_significant", "_p.adj", contrast)] < alpha,]
        if (database == "KEGG") {
          mappings <- bitr(gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
          result <- enrichKEGG(gene = mappings$ENTREZID,
                               keyType = "ncbi-geneid",
                               organism = "hsa",
                               pvalueCutoff = 1,
                               qvalueCutoff = 1)
          bg_mappings <- bitr(background, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
          bg_result <- enrichKEGG(gene = bg_mappings$ENTREZID,
                                 keyType = "ncbi-geneid",
                                 organism = "hsa",
                                 pvalueCutoff = 1,
                                 qvalueCutoff = 1)
          temp <- as.data.frame(setReadable(result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"))
          bg_temp <- as.data.frame(bg_result)
          bg_temp <- bg_temp %>%
            mutate(bg_IN = Count,
                   bg_OUT = length(background) - Count) %>%
            select(ID, bg_IN, bg_OUT)
          temp <- temp %>% left_join(bg_temp, by = "ID")
        } else {
          if (database == "Hallmark") {
            hallmark <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
              dplyr::select(gs_name, gene_symbol)
            result <- enricher(gene, TERM2GENE = hallmark, pvalueCutoff = 1, qvalueCutoff = 1)
            bg_result <- enricher(background, TERM2GENE = hallmark, pvalueCutoff = 1, qvalueCutoff = 1)
          } else{
            database_mappings <- c(
              "GO Biological Process" = "BP",
              "GO Cellular Component" = "CC",
              "GO Molecular Function" = "MF"
            )
            result <- enrichGO(gene = gene,
                               OrgDb = org.Hs.eg.db,
                               ont = database_mappings[database],
                               keyType = "SYMBOL",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 1,
                               qvalueCutoff  = 1)
            bg_result <- enrichGO(gene = background,
                                  OrgDb = org.Hs.eg.db,
                                  ont = database_mappings[database],
                                  keyType = "SYMBOL",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff  = 1,
                                  qvalueCutoff  = 1)
          }
          temp <- as.data.frame(result)
          bg_temp <- as.data.frame(bg_result)
          bg_temp <- bg_temp %>%
            mutate(bg_IN = Count,
                   bg_OUT = length(background) - bg_IN) %>%
            select(ID, bg_IN, bg_OUT)
          temp <- temp %>% left_join(bg_temp, by = "ID")
        }
        temp$contrast <- gsub("_significant", "", contrast)
        temp$OUT <- dim(significant)[1] - temp$Count
        combined_df <- rbind(combined_df, temp)
      }
      lookup <- c("Term"="Description", "P.value"="pvalue", "Adjusted.P.value"="p.adjust", "IN"="Count")
      combined_df <- rename(combined_df, all_of(lookup))
      combined_df$var <- database
      combined_df$Overlap <- paste0(combined_df$IN, "/", as.numeric(gsub("/.*", "", combined_df$BgRatio)))
      # combined_df$Odds.Ratio <- (combined_df$IN * (as.numeric(gsub(".*/", "", combined_df$BgRatio)) - as.numeric(gsub("/.*", "", combined_df$BgRatio)))) / (combined_df$OUT * as.numeric(gsub("/.*", "", combined_df$BgRatio)))
      # inspired from enrichr https://github.com/MaayanLab/enrichr_issues/issues/3#issuecomment-780078054
      combined_df$Odds.Ratio <- (combined_df$IN * (20000 - as.numeric(gsub("/.*", "", combined_df$BgRatio)))) / (combined_df$OUT * as.numeric(gsub("/.*", "", combined_df$BgRatio)))
      combined_df$log_odds <- log2((combined_df$IN * combined_df$bg_OUT) / (combined_df$OUT * combined_df$bg_IN))
      combined_df$p_hyper = phyper(q=(combined_df$IN-1), m = combined_df$bg_IN, n = combined_df$bg_OUT, k = (combined_df$IN+combined_df$OUT),
                                 lower.tail = F)
      combined_df$p.adjust_hyper = p.adjust(combined_df$p_hyper, method = "BH")
      return(combined_df)
    } else {
      return(NULL)
    }
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
                         databases,
                         contrasts = T, direction="UP", log2_threshold=0.7, alpha=0.05) {
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
  } else if (metadata(dep)$level %in% c("peptide", "site", "glycan")) {
    background <- unique(row_data$Gene)
  }

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
      df[is.na(df[[contrast]]),contrast] <- F
      significant <- df
      significant <- significant[!is.na(significant[gsub("_significant", "_diff", contrast)]),]
      if (direction == "UP"){
        significant <- significant[(significant[gsub("_significant", "_diff", contrast)] > log2_threshold),]
      } else if (direction == "DOWN") {
        significant <- significant[significant[gsub("_significant", "_diff", contrast)] < -log2_threshold,]
      }
      significant <- significant[significant[gsub("_significant", "_p.adj", contrast)] < alpha,]

      if (metadata(dep)$exp == "TMT" & metadata(dep)$level == "protein") {
        genes <- unique(significant$Gene)
      } else if (metadata(dep)$exp == "TMT" & metadata(dep)$level == "gene") {
        genes <- unique(significant$ID)
      } else if (metadata(dep)$level == "protein") {
        genes <- significant$name
      } else if (metadata(dep)$level %in% c("peptide", "site", "glycan")) {
        genes <- unique(significant$Gene)
      }

      message(paste0(length(genes), " genes are submitted"))
      if (length(genes) != 0){
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
    } else if (metadata$level %in% c("peptide", "site", "glycan")) {
      genes <- unique(significant$Gene)
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
  df_enrich$p_hyper = phyper(q=(df_enrich$IN-1), m = df_enrich$bg_IN, n = df_enrich$bg_OUT, k = (df_enrich$IN+df_enrich$OUT),
                             lower.tail = F )
  df_enrich$p.adjust_hyper = p.adjust(df_enrich$p_hyper, method = "BH")
  return(df_enrich)
}

enrichr_mod <- function(genes, databases = NULL) {
  # check gene type
  if (length(genes) != 0){
    if (all(startsWith(genes, "ENSG"))) {
      genes <- gsub("\\..*", "", genes)
      genes_map <- ensembldb::select(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
                                     keys= genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
      genes <- genes_map$SYMBOL
    }
    httr::set_config(httr::config(ssl_verifypeer = 0L))
    cat("Uploading data to Enrichr... ")
    if (is.vector(genes) & ! all(genes == "") & length(genes) != 0) {
      temp <- httr::POST(url="http://maayanlab.cloud/Enrichr/enrich",
                   body=list(list=paste(genes, collapse="\n")))
    } else if (is.data.frame(genes)) {
      temp <- httr::POST(url="http://maayanlab.cloud/Enrichr/enrich",
                   body=list(list=paste(paste(genes[,1], genes[,2], sep=","),
                                        collapse="\n")))
    } else {
      warning("genes must be a non-empty vector of gene names or a dataframe with genes and score.")
    }
    GET(url="http://maayanlab.cloud/Enrichr/share")
    cat("Done.\n")
    dbs <- as.list(databases)
    dfSAF <- options()$stringsAsFactors
    options(stringsAsFactors = FALSE)
    result <- lapply(dbs, function(x) {
      cat("  Querying ", x, "... ", sep="")
      r <- httr::GET(url="http://maayanlab.cloud/Enrichr/export",
               query=list(file="API", backgroundType=x))
      r <- gsub("&#39;", "'", intToUtf8(r$content))
      tc <- textConnection(r)
      r <- read.table(tc, sep = "\t", header = TRUE, quote = "", comment.char="")
      close(tc)
      cat("Done.\n")
      return(r)
    })
    options(stringsAsFactors = dfSAF)
    cat("Parsing results... ")
    names(result) <- dbs
    cat("Done.\n")
  }
  else { # no genes provided
    result <- data.frame(Term=character(),
                         Overlap=character(),
                         P.value=double(),
                         Adjusted.P.value=double(),
                         Old.P.value=double(),
                         Old.Adjusted.P.value=double(),
                         Odds.Ratio=double(),
                         Combined.Score=double(),
                         Genes=character())
  }
  return(result)
}

#' Prepare data for PTM-SEA analysis
#'
#' Exports phosphoproteomics data in GCT format for PTM Signature Enrichment
#' Analysis (PTM-SEA) using ssGSEA.
#'
#' @param se A \code{SummarizedExperiment} object containing phosphoproteomics data.
#'   Must have \code{exp_type = "phospho"} in metadata and a "SequenceWindow"
#'   column in rowData.
#' @param col Character string specifying the column name in rowData containing
#'   the values to export (e.g., log2 fold change column).
#' @param outfile Character string specifying the output file path for the GCT file.
#'
#' @return Writes a GCT format file to the specified path. Returns \code{NULL} invisibly.
#'
#' @details
#' PTM-SEA is a method for identifying dysregulated PTM signatures based on
#' phosphoproteomics data. The output GCT file can be used with the ssGSEA
#' module in GenePattern or other compatible tools.
#'
#' The function creates peptide identifiers by converting sequence windows to
#' uppercase and appending "-p" to indicate phosphorylation sites.
#'
#' @examples
#' \dontrun{
#' # Prepare data for PTM-SEA analysis
#' prepare_PTMSEA(se_phospho, col = "Treatment_vs_Control_diff",
#'                outfile = "ptmsea_input.gct")
#' }
#'
#' @seealso \code{\link{visualize_PTMSEA}}
#'
#' @importFrom cmapR GCT
#'
#' @export
prepare_PTMSEA <- function(se, col, outfile) {
  if (metadata(se)$exp_type == "phospho") {
    temp <- data.frame(as.data.frame(rowData(se))["SequenceWindow"], as.data.frame(rowData(se))[col])
    colnames(temp) <- c("peptide", col)
    temp <- temp[!is.na(temp[[col]]),]
    gct <- new("GCT", mat=as.matrix(temp[,col, drop=F]), rid=paste0(make.unique(toupper(temp$peptide)), "-p"))
    write_gct(gct, outfile, appenddim = F)
  }
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
#'   GCT file (typically generated by ssGSEA).
#' @param col Character string specifying the score column to visualize.
#' @param selected_concepts Character vector of specific concept/signature IDs
#'   to display. Default is \code{NULL} (show top enriched concepts).
#' @param num_concepts Numeric value specifying the number of top and bottom
#'   enriched concepts to display when \code{selected_concepts} is \code{NULL}.
#'   Default is 5 (showing 10 total).
#'
#' @return A \code{ggplot} object showing a dot plot where:
#'   \itemize{
#'     \item X-axis: Enrichment score
#'     \item Y-axis: PTM signature names
#'     \item Point size: Signature set overlap percentage
#'     \item Point color: FDR-adjusted p-value
#'   }
#'   Returns \code{NULL} if the specified column is not found.
#'
#' @examples
#' \dontrun{
#' # Visualize PTM-SEA results
#' p <- visualize_PTMSEA("ptmsea_output.gct", col = "Treatment_vs_Control")
#'
#' # Show specific kinase signatures
#' p <- visualize_PTMSEA("ptmsea_output.gct", col = "Treatment_vs_Control",
#'                       selected_concepts = c("KINASE_A", "KINASE_B"))
#' }
#'
#' @seealso \code{\link{prepare_PTMSEA}}
#'
#' @importFrom cmapR parse_gctx mat meta
#' @importFrom ggplot2 ggplot aes_string geom_point scale_size_continuous
#'   scale_color_continuous theme_bw
#'
#' @export
visualize_PTMSEA <- function(gct_file, col, selected_concepts=NULL, num_concepts=5) {
  gct <- parse_gctx(gct_file)
  score_cols <- colnames(mat(gct))
  if (!col %in% score_cols) {
    return(NULL)
  }
  data <- cbind(meta(gct, dimension = "row"), mat(gct))
  if (is.null(selected_concepts)) {
    data <- data[order(data[[col]], decreasing = T),]
    data <- data[c(c(1:num_concepts), (dim(data)[1] - num_concepts + 1):dim(data)[1]),]
  } else {
    data <- data[data$id %in% selected_concepts,]
    data <- data[order(data[[col]], decreasing = F),]
  }
  data$id <- factor(data$id, levels=data$id)
  p <- ggplot(data, aes_string(x = col, y = "id", color=paste0("fdr.pvalue.", col), size=paste0("Signature.set.overlap.percent.", col))) +
    geom_point() +
    scale_size_continuous(range = c(0.5, 11), name="Overlap percentage") +
    scale_color_continuous(low="red", high="blue", name = "p.adjust",
                           guide=guide_colorbar(reverse=T, label.vjust = 0.5)) +
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
#' @importFrom ggplot2 ggplot aes geom_point theme_bw theme
#'   element_blank element_line labs
#' @importFrom ggrepel geom_text_repel
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
