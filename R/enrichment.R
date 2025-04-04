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

#' @export
plot_GSEA <- function(gsea_result, categroies=15) {
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
    write.table(temp, outfile, quote = F, row.names = F, col.names = F, sep="\t")
  } else {
    cat("Kinome analysis requires phosphorylation dataset. \n")
  }
}

#' @export
visualize_kinome <- function(tsv_file, labels=NULL) {
  data <- read.csv(tsv_file, sep = "\t", stringsAsFactors = F)
  if (is.null(labels)){
    p <- ggplot(data, aes(x=dominant_enrichment_value_log2, y=dominant_adjusted_p_value_log10_abs, label=kinase)) +
      geom_point() +
      geom_text_repel(data=subset(data, abs(dominant_enrichment_value_log2) > 1 & (dominant_adjusted_p_value_log10_abs) > 1.3)) +
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
