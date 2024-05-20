# FragPipeAnalystR

R package for [FragPipe](https://fragpipe.nesvilab.org/) downstream analysis. FragPipeAnalystR is the core of [FragPipe-Analyst Shiny App](https://fragpipe-analyst-doc.nesvilab.org/).

## Installation

We recommend to use [renv](https://rstudio.github.io/renv/index.html) for managing your R environment dependencies

``` r
renv::install("bioc::ComplexHeatmap")
renv::install("bioc::limma")
renv::install("bioc::MSnbase")
renv::install("bioc::SummarizedExperiment")
renv::install("bioc::cmapR")
renv::install("bioc::ConsensusClusterPlus")
renv::install("Nesvilab/FragPipeAnalystR")

# optional
renv::install("nicolerg/ssGSEA2")
```

## Tutorial

- [LFQ AP-MS protein level analysis](global_LFQ_prot_tutorial.html)
- [DIA ccRCC protein level analysis](global_DIA_prot_tutorial.html)

## Reference
Analysis and visualization of quantitative proteomics data using FragPipe-Analyst
Yi Hsiao, Haijian Zhang, Ginny Xiaohe Li, Yamei Deng, Fengchao Yu, Hossein Valipour Kahrood, Joel R. Steele, Ralf B. Schittenhelm, Alexey I. Nesvizhskii
bioRxiv 2024.03.05.583643; doi: https://doi.org/10.1101/2024.03.05.583643
