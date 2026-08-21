# These symbols are column names referenced via non-standard evaluation
# (dplyr/data.table verbs, ggplot2 aes() mappings, magrittr pipe placeholder)
# rather than actual global variables. Declaring them here silences the
# `R CMD check` NOTEs "no visible binding for global variable" without
# affecting runtime behavior.
utils::globalVariables(c(
  ".", ".SD", "rowname",
  "contrast", "var", "p.adjust_hyper", "Term", "p_hyper",
  "Adjusted.P.value", "P.value", "log_odds", "IN", "Odds.Ratio",
  "ID", "ID_new", "label", "name",
  "correlation",
  "UMAP1", "UMAP2",
  "NA_count", "missPercent", "FeatureCount",
  "plex", "total_psm",
  "ProteinID", "Index", "caseID"
))
