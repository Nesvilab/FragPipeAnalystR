# FragPipeAnalystR

R package for FragPipe downstream analysis

## Example

```
library(FragPipeR)
data("ccrcc", package = "FragPipeR")
```

```
plot_pca(ccrcc, n=500, ID_col="label", exp="TMT")
```

![PCA](vignettes/tutorial_files/figure-html/unnamed-chunk-2-1.png)

```
plot_missval_heatmap(ccrcc)
```

![missing value heatmap](vignettes/tutorial_files/figure-html/unnamed-chunk-4-1.png)

```
plot_correlation_heatmap(ccrcc, indicate="condition", exp="TMT")
```

![correlation heatmap](vignettes/tutorial_files/figure-html/unnamed-chunk-3-1.png)


## Developer guide
- `devtools::build()` and `devtools::test()`
- `devtools::document()` when you add a function
- `usethis::use_import_from("package name", "function name")` to include dependencies
- Use `lintr` and `styler` to improve code quality
```
lintr::lint_package()
styler::style_pkg()
```
