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
# latest version
renv::install("Nesvilab/FragPipeAnalystR@1.0.2")
# or dev version
renv::install("Nesvilab/FragPipeAnalystR")

# optional
renv::install("nicolerg/ssGSEA2")
```

## Get started
To use FragPipeAnalystR, follow the same instruction for [preparing the input file for FragPipe-Analyst](https://fragpipe-analyst-doc.nesvilab.org/Formatting.html).

## Tutorial

- [LFQ AP-MS protein level analysis](global_LFQ_prot_tutorial.html) 
- [TMT ccRCC protein level analysis](global_TMT_prot_tutorial.html)
- [DIA ccRCC protein level analysis](global_DIA_prot_tutorial.html)

- [TMT ccRCC phosphosite analysis](phospho_TMT_tutorial.html)

## Reference
[Yi Hsiao, Haijian Zhang, Ginny Xiaohe Li, Yamei Deng, Fengchao Yu, Hossein Valipour Kahrood, Joel R. Steele, Ralf B. Schittenhelm, and Alexey I. Nesvizhskii
Journal of Proteome Research (2024), DOI: 10.1021/acs.jproteome.4c00294](https://pubs.acs.org/doi/10.1021/acs.jproteome.4c00294)
