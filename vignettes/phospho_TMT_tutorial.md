---
title: "Phosphopeptide Analysis using FragPipeAnalystR"
output:
  html_document:
    keep_md: yes
---


``` r
library(FragPipeAnalystR)
```

## Introduction
One of the reasons we created `FragPipeAnalystR` is to support peptide-level analysis of post-translational modifications (PTM). Here, we used three plex sets of TMT phosphoproteomics experiemnt from the previously published [clear cell renal cell carcinoma study (ccRCC)](https://www.sciencedirect.com/science/article/pii/S0092867419311237?via%3Dihub). After the phosphoproteomics data is processed by FragPipe, you should be able to get the reports produced by TMT-Integrator. Or you can download the report used in the tutorial [here](https://drive.google.com/drive/u/2/folders/1x8DCxGKdsZQmYGsv3vSWfj1IYRPzykmT).

## Read the data

The first step is always the same to read the data. Note that we specify the `level` and `exp_type` here.


``` r
se <- make_se_from_files("/Users/hsiaoyi/Documents/workspace/FragPipeR_manuscript/data/TMT_4plex/ratio_protein_MD.tsv", "/Users/hsiaoyi/Documents/workspace/FragPipeR_manuscript/data/TMT_4plex/experiment_annotation_clean.tsv", level="protein", type = "TMT")
```


``` r
se_phospho <- make_se_from_files("/Users/hsiaoyi/Documents/workspace/FragPipeR_manuscript/data/TMT_phospho_4plex/ratio_single-site_MD.tsv", "/Users/hsiaoyi/Documents/workspace/FragPipeR_manuscript/data/TMT_4plex/experiment_annotation_clean.tsv", level="peptide", type = "TMT", exp_type="phospho")
```

## QC plots

Then, you should able to observe the data through various plots supported. For example, the PCA plot 

``` r
plot_pca(se_phospho)
```

![](phospho_TMT_tutorial_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

or heatmap

``` r
plot_correlation_heatmap(se_phospho)
```

![](phospho_TMT_tutorial_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

Both of these plots show that protein phosphorylation status of tumor and normal adjacent tumors (NATs) is quite different.

## Normalization

``` r
normalized_se <- PTM_normalization(se_phospho, se)
```

```
## correlation with prot after normalization:  -1.926138e-16 
## get sites_windows 
## get dm 
## get subpsite
```

``` r
pca_plot_normalized <- plot_pca(normalized_se, n=0, ID_col="label", exp="TMT", interactive = F)
```

### Boxplot before normalization

``` r
plot_feature(se_phospho, c("P14618_Y148", # PKM_Y146
                           "P28482_Y187",
                           "Q13541_S65"))
```

![](phospho_TMT_tutorial_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

### Boxplot after normalization

``` r
plot_feature(normalized_se, c("P14618_Y148", # PKM_Y146
                                          "P28482_Y187",
                                          "Q13541_S65"))
```

![](phospho_TMT_tutorial_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

### Boxplot of parent proteins

``` r
plot_feature(se, c("P14618",
                   "P28482",
                   "Q13541"))
```

![](phospho_TMT_tutorial_files/figure-html/unnamed-chunk-8-1.png)<!-- -->


## Differential expression (DE) analysis

One of the frequent analysis we do is differential expression analysis to understand the difference between tumor and NAT. It can be performed through following commands:


``` r
de_result <-test_limma(se_phospho, type = "all")
```

```
## Tested contrasts: Tumor_vs_NAT
```

``` r
de_result_updated <- add_rejections(de_result)
plot_volcano(de_result_updated, "Tumor_vs_NAT")
```

![](phospho_TMT_tutorial_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

As you can see, there are many phosphosites upregulated and downregulated in this comparison. This gives a lot of research opportunities. For example, CAK2 (P51636) has been associated with [maintaining kidney cancer malignant](https://pubmed.ncbi.nlm.nih.gov/30288056/). Further investigating this phoshphosite (S36) might help us understand its mechanism. One thing that needs to be noted here is that the abundance of phosphopeptides is usually correlated with its protein abundance, so it might be just because the protein abundance of CAK2 is upregulated in ccRCC. It might be worth to check with the global proteome available as well. 


``` r
normalized_de_result <- test_limma(normalized_se, type = "all")
```

```
## Tested contrasts: Tumor_vs_NAT
```

``` r
normalized_de_result_updated <- add_rejections(normalized_de_result)
plot_volcano(de_result_updated, "Tumor_vs_NAT")
```

![](phospho_TMT_tutorial_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

## Enrichment analysis
One of the key difference between peptide-level analysis of PTMs we demonstrated here and protein level analysis is that PTMs usually act in a site-specific manner. Here, we provide the way to help users perform enrichment analysis site-specifically through creating the input file you needed for [PTM-SEA](https://doi.org/10.1074/mcp.tir118.000943).

``` r
prepare_PTMSEA(normalized_de_result_updated, "Tumor_vs_NAT_diff", "/Users/hsiaoyi/Documents/workspace/FragPipeR_manuscript/result/PTMSEA_new/result.gct")
```

```
## Saving file to  /Users/hsiaoyi/Documents/workspace/FragPipeR_manuscript/result/PTMSEA_new/result.gct 
## Dimensions of matrix: [34051x1]
## Setting precision to 4
## Saved.
```

Then, you can run PTM-SEA though [ssGSEA2 R package](https://github.com/nicolerg/ssGSEA2) like this:

``` r
library(ssGSEA2)
res <- run_ssGSEA2("/Users/hsiaoyi/Documents/workspace/FragPipeR_manuscript/result/PTMSEA_new/result.gct",
                   output.prefix = "ccRCC",
                   gene.set.databases = "/Users/hsiaoyi/Documents/workspace/FragPipeR_manuscript/db_ptm.sig.db.all.v2.0.0/ptm.sig.db.all.flanking.human.v2.0.0.gmt",
                   output.directory = "/Users/hsiaoyi/Documents/workspace/FragPipeR_manuscript/result/PTMSEA_new/result",
                   sample.norm.type = "none", 
                   weight = 0.75, 
                   correl.type = "rank", 
                   statistic = "area.under.RES",
                   output.score.type = "NES", 
                   nperm = 1000, 
                   min.overlap = 5, 
                   extended.output = T,
                   global.fdr = FALSE,
                   par=T,
                   spare.cores=4,
                   log.file ="/Users/hsiaoyi/Documents/workspace/FragPipeR_manuscript/result/PTMSEA_new/ccRCC_PTMSEA.log")
```

and visualize the result through:


``` r
visualize_PTMSEA("/Users/hsiaoyi/Documents/workspace/FragPipeR_manuscript/result/PTMSEA_new/result/ccRCC-combined.gct", "Tumor_vs_NAT_diff")
```

```
## parsing as GCT v1.3
```

```
## /Users/hsiaoyi/Documents/workspace/FragPipeR_manuscript/result/PTMSEA_new/result/ccRCC-combined.gct 603 rows, 1 cols, 7 row descriptors, 0 col descriptors
```

![](phospho_TMT_tutorial_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

You may also select particular gene set of interests. For example, following code snippets demonstrate the result on PKC and AKT protein kinases, respectively.


``` r
visualize_PTMSEA("/Users/hsiaoyi/Documents/workspace/FragPipeR_manuscript/result/PTMSEA_new/result/ccRCC-combined.gct",
                                  "Tumor_vs_NAT_diff",
                                  selected_concepts=c("KINASE-PSP_PKCA/PRKCA", "KINASE-PSP_PKCB/PRKCB", "KINASE-PSP_PKCG/PRKCG", "KINASE-PSP_PKCB_iso2/PRKCB", "KINASE-PSP_PKCD/PRKCD", "KINASE-PSP_PKCI/PRKCI", "KINASE-PSP_PKCT/PRKCQ",
                                                      "KINASE-PSP_PKCH/PRKCH", "KINASE-PSP_PKCE/PRKCE", "KINASE-PSP_PKCZ/PRKCZ"))
```

```
## parsing as GCT v1.3
```

```
## /Users/hsiaoyi/Documents/workspace/FragPipeR_manuscript/result/PTMSEA_new/result/ccRCC-combined.gct 603 rows, 1 cols, 7 row descriptors, 0 col descriptors
```

![](phospho_TMT_tutorial_files/figure-html/unnamed-chunk-14-1.png)<!-- -->


``` r
visualize_PTMSEA("/Users/hsiaoyi/Documents/workspace/FragPipeR_manuscript/result/PTMSEA_new/result/ccRCC-combined.gct",
                 "Tumor_vs_NAT_diff",
                 selected_concepts=c("KINASE-PSP_Akt1/AKT1", "KINASE-PSP_Akt3/AKT3", "KINASE-PSP_Akt2/AKT2",
                                     "KINASE-iKiP_AKT3", "KINASE-iKiP_AKT2", "PATH-WP_PI3K-Akt_signaling_pathway", "KINASE-iKiP_AKT1"))
```

```
## parsing as GCT v1.3
```

```
## /Users/hsiaoyi/Documents/workspace/FragPipeR_manuscript/result/PTMSEA_new/result/ccRCC-combined.gct 603 rows, 1 cols, 7 row descriptors, 0 col descriptors
```

![](phospho_TMT_tutorial_files/figure-html/unnamed-chunk-15-1.png)<!-- -->
                                                      
Alternatively, FragPipeAnalystR also provides ways to generate input and visualization for the [Kinase Library](https://kinase-library.phosphosite.org/). Here is the example of generating the input file:


``` r
prepare_kinome(de_result_updated, "Tumor_vs_NAT_diff", "/Users/hsiaoyi/Documents/workspace/FragPipeR_manuscript/result/kinome_new/kinome_input_asterisk.tsv")
```


And here is the result visualization:

``` r
visualize_kinome("/Users/hsiaoyi/Documents/workspace/FragPipeR_manuscript/result/kinome_result/enrichment-analysis-result-table.txt")
```

![](phospho_TMT_tutorial_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

Similar to `visualize_PTMSEA`, users could specify the concepts of interests.


``` r
AKT_labels <- c("AKT1", "AKT2", "ATK3")
visualize_kinome("/Users/hsiaoyi/Documents/workspace/FragPipeR_manuscript/result/kinome_result/enrichment-analysis-result-table.txt",
                 labels=AKT_labels)
```

![](phospho_TMT_tutorial_files/figure-html/unnamed-chunk-18-1.png)<!-- -->


``` r
PKC_labels <- c("PKCI", "PKCH", "PKCT", "PKCZ", "PKCE", "PKCD", "PKCG", "PKCA", "PKCB")
visualize_kinome("/Users/hsiaoyi/Documents/workspace/FragPipeR_manuscript/result/kinome_result/enrichment-analysis-result-table.txt",
                 labels=PKC_labels)
```

![](phospho_TMT_tutorial_files/figure-html/unnamed-chunk-19-1.png)<!-- -->




