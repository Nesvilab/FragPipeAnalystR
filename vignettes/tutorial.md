---
title: "FragPipeAnalystR Tutorial"
output:
  html_document:
    keep_md: yes
---

# Introduction

FragPipeAnalystR is a R package intended for FragPipe downstream analysis. We also make it compatible with the result obtained from FragPipe-Analyst. Users are able to reproduce and customize the plots generated in FragPipe-Analyst.

## Quick Start Example

```r
library(FragPipeAnalystR)
data("ccrcc", package = "FragPipeAnalystR")
```


```r
plot_pca(ccrcc, n=500, ID_col="label", exp="TMT")
```

![](tutorial_files/figure-html/unnamed-chunk-2-1.png)<!-- -->


```r
plot_correlation_heatmap(ccrcc, indicate="condition", exp="TMT")
```

![](tutorial_files/figure-html/unnamed-chunk-3-1.png)<!-- -->


```r
plot_missval_heatmap(ccrcc)
```

![](tutorial_files/figure-html/unnamed-chunk-4-1.png)<!-- -->
