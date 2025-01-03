---
title: "FragPipeAnalystR TMT Global Proteome"
output:
  html_document:
    keep_md: yes
---

# Introduction

FragPipeAnalystR is a R package intended for FragPipe downstream analysis. We also make it compatible with the result obtained from FragPipe-Analyst. Users are able to reproduce and customize the plots generated in FragPipe-Analyst.

## Quick Start Example

``` r
library(FragPipeAnalystR)
ccrcc <- make_se_from_files("/Users/hsiaoyi/Documents/workspace/FragPipeR_manuscript/data/TMT_4plex/abundance_protein_MD.tsv",
                         "/Users/hsiaoyi/Documents/workspace/FragPipeR_manuscript/data/TMT_4plex/experiment_annotation_clean.tsv",
                         type = "TMT", level = "protein")
```


``` r
plot_pca(ccrcc)
```

![](global_TMT_prot_tutorial_files/figure-html/unnamed-chunk-2-1.png)<!-- -->


``` r
plot_correlation_heatmap(ccrcc)
```

![](global_TMT_prot_tutorial_files/figure-html/unnamed-chunk-3-1.png)<!-- -->


``` r
plot_missval_heatmap(ccrcc)
```

```
## `use_raster` is automatically set to TRUE for a matrix with more than
## 2000 rows. You can control `use_raster` argument by explicitly setting
## TRUE/FALSE to it.
## 
## Set `ht_opt$message = FALSE` to turn off this message.
```

```
## 'magick' package is suggested to install to give better rasterization.
## 
## Set `ht_opt$message = FALSE` to turn off this message.
```

![](global_TMT_prot_tutorial_files/figure-html/unnamed-chunk-4-1.png)<!-- -->


``` r
plot_feature_numbers(ccrcc)
```

![](global_TMT_prot_tutorial_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

You may want to check some of known markers through box plots:


``` r
plot_feature(ccrcc,  c("Q16790", # CA9
                      "Q8IVF2", # AHNAK2
                      "P19404", # NDUFV2
                      "P01833" # PIGR
                      ))
```

![](global_TMT_prot_tutorial_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

This could be done via `Gene` column as well:


``` r
plot_feature(ccrcc, c("CA9", "AHNAK2", "NDUFV2", "PIGR"), index="Gene")
```

![](global_TMT_prot_tutorial_files/figure-html/unnamed-chunk-7-1.png)<!-- -->


``` r
de_result <- test_limma(ccrcc, type = "all")
```

```
## Tested contrasts: Tumor_vs_NAT
```

``` r
de_result_updated <- add_rejections(de_result)
```

Volcano plot is designed for visualizing differential expression analysis result:


``` r
plot_volcano(de_result_updated, "Tumor_vs_NAT")
```

![](global_TMT_prot_tutorial_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

It could be labelled by different column available in the `rowData(de_result_updated)` such as `Gene`:


``` r
plot_volcano(de_result_updated, "Tumor_vs_NAT", name_col="Gene")
```

![](global_TMT_prot_tutorial_files/figure-html/unnamed-chunk-10-1.png)<!-- -->


``` r
or_result <- or_test(de_result_updated, database = "Hallmark", direction = "UP")
```

```
## Background
```

```
## Uploading data to Enrichr... Done.
##   Querying MSigDB_Hallmark_2020... Done.
## Parsing results... Done.
```

```
## Tumor_vs_NAT
```

```
## 774 genes are submitted
```

```
## Uploading data to Enrichr... Done.
##   Querying MSigDB_Hallmark_2020... Done.
## Parsing results... Done.
## Background correction... Done.
```

``` r
plot_or(or_result)
```

![](global_TMT_prot_tutorial_files/figure-html/unnamed-chunk-11-1.png)<!-- -->


``` r
or_result <- or_test(de_result_updated, database = "Hallmark", direction = "DOWN")
```

```
## Background
```

```
## Uploading data to Enrichr... Done.
##   Querying MSigDB_Hallmark_2020... Done.
## Parsing results... Done.
```

```
## Tumor_vs_NAT
```

```
## 1432 genes are submitted
```

```
## Uploading data to Enrichr... Done.
##   Querying MSigDB_Hallmark_2020... Done.
## Parsing results... Done.
## Background correction... Done.
```

``` r
plot_or(or_result)
```

![](global_TMT_prot_tutorial_files/figure-html/unnamed-chunk-12-1.png)<!-- -->


``` r
sessionInfo()
```

```
## R version 4.3.1 Patched (2023-10-12 r85331)
## Platform: aarch64-apple-darwin20 (64-bit)
## Running under: macOS Ventura 13.4
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRblas.0.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## time zone: America/Detroit
## tzcode source: internal
## 
## attached base packages:
## [1] stats     graphics  grDevices datasets  utils     methods   base     
## 
## other attached packages:
## [1] FragPipeAnalystR_0.1.5
## 
## loaded via a namespace (and not attached):
##   [1] bitops_1.0-7                fdrtool_1.2.17             
##   [3] rlang_1.1.3                 magrittr_2.0.3             
##   [5] clue_0.3-65                 GetoptLong_1.0.5           
##   [7] matrixStats_1.3.0           compiler_4.3.1             
##   [9] png_0.1-8                   vctrs_0.6.5                
##  [11] stringr_1.5.1               ProtGenerics_1.34.0        
##  [13] pkgconfig_2.0.3             shape_1.4.6.1              
##  [15] crayon_1.5.2                fastmap_1.2.0              
##  [17] XVector_0.42.0              labeling_0.4.3             
##  [19] utf8_1.2.4                  rmarkdown_2.27             
##  [21] tzdb_0.4.0                  preprocessCore_1.64.0      
##  [23] purrr_1.0.2                 xfun_0.44                  
##  [25] zlibbioc_1.48.2             cachem_1.1.0               
##  [27] SNFtool_2.3.1               GenomeInfoDb_1.38.8        
##  [29] jsonlite_1.8.8              ExPosition_2.8.23          
##  [31] highr_0.10                  DelayedArray_0.28.0        
##  [33] BiocParallel_1.36.0         parallel_4.3.1             
##  [35] cluster_2.1.4               R6_2.5.1                   
##  [37] stringi_1.8.4               bslib_0.7.0                
##  [39] RColorBrewer_1.1-3          limma_3.58.1               
##  [41] GenomicRanges_1.54.1        jquerylib_0.1.4            
##  [43] assertthat_0.2.1            Rcpp_1.0.12                
##  [45] SummarizedExperiment_1.32.0 iterators_1.0.14           
##  [47] knitr_1.46                  readr_2.1.5                
##  [49] flowCore_2.14.2             IRanges_2.36.0             
##  [51] Matrix_1.6-1.1              tidyselect_1.2.1           
##  [53] rstudioapi_0.16.0           abind_1.4-5                
##  [55] yaml_2.3.8                  doParallel_1.0.17          
##  [57] codetools_0.2-19            affy_1.80.0                
##  [59] curl_5.2.1                  lattice_0.21-9             
##  [61] tibble_3.2.1                plyr_1.8.9                 
##  [63] withr_3.0.0                 Biobase_2.62.0             
##  [65] evaluate_0.23               ConsensusClusterPlus_1.66.0
##  [67] circlize_0.4.16             pillar_1.9.0               
##  [69] affyio_1.72.0               BiocManager_1.30.23        
##  [71] MatrixGenerics_1.14.0       renv_0.17.0                
##  [73] foreach_1.5.2               stats4_4.3.1               
##  [75] plotly_4.10.4               MSnbase_2.28.1             
##  [77] MALDIquant_1.22.2           ncdf4_1.22                 
##  [79] generics_0.1.3              RCurl_1.98-1.14            
##  [81] hms_1.1.3                   S4Vectors_0.40.2           
##  [83] ggplot2_3.5.1               munsell_0.5.1              
##  [85] scales_1.3.0                glue_1.7.0                 
##  [87] lazyeval_0.2.2              tools_4.3.1                
##  [89] data.table_1.15.4           mzID_1.40.0                
##  [91] vsn_3.70.0                  mzR_2.36.0                 
##  [93] XML_3.99-0.16.1             grid_4.3.1                 
##  [95] impute_1.76.0               tidyr_1.3.1                
##  [97] RProtoBufLib_2.14.1         prettyGraphs_2.1.6         
##  [99] MsCoreUtils_1.14.1          colorspace_2.1-0           
## [101] GenomeInfoDbData_1.2.11     cmapR_1.14.0               
## [103] cli_3.6.2                   fansi_1.0.6                
## [105] viridisLite_0.4.2           cytolib_2.14.1             
## [107] S4Arrays_1.2.1              ComplexHeatmap_2.18.0      
## [109] dplyr_1.1.4                 pcaMethods_1.94.0          
## [111] gtable_0.3.5                sass_0.4.9                 
## [113] digest_0.6.35               BiocGenerics_0.48.1        
## [115] ggrepel_0.9.5               SparseArray_1.2.4          
## [117] farver_2.1.2                htmlwidgets_1.6.4          
## [119] rjson_0.2.21                htmltools_0.5.8.1          
## [121] lifecycle_1.0.4             httr_1.4.7                 
## [123] alluvial_0.1-2              GlobalOptions_0.1.2        
## [125] statmod_1.5.0               MASS_7.3-60
```
