# 🌱 Tumor Microenviroment Analysis

## Availiable TME Analysis Method in tigeR

10 open source tumor microenvironment deconvolution methods are indlcued in tigeR.

::: {style="width:780px; height:200px; overflow-y: scroll; overflow-x: hidden;"}
|                                 Algorithm                                  |                                                    license                                                    |                         PMID                          |  allow custom signature matrix  |
|:--------------------:|:------------------------------:|:----------------:|:--------------:|
|                    [TIMER](http://cistrome.org/TIMER/)                     |                                                free (GPL 2.0)                                                 | [27549193](https://pubmed.ncbi.nlm.nih.gov/27549193/) | √ |
|                [CIBERSORT](https://cibersort.stanford.edu/)                |                                       free for non-commercial use only                                        | [25822800](https://pubmed.ncbi.nlm.nih.gov/25822800/) | √ |
|             [MCPCounter](https://github.com/ebecht/MCPcounter)             |               free ([GPL 3.0](https://github.com/ebecht/MCPcounter/blob/master/Source/License))               | [27765066](https://pubmed.ncbi.nlm.nih.gov/27765066/) | √ |
|                      [xCell](http://xcell.ucsf.edu/)                       |                  free ([GPL 3.0](https://github.com/dviraran/xCell/blob/master/DESCRIPTION))                  | [29141660](https://pubmed.ncbi.nlm.nih.gov/29141660/) | × |
|             [IPS](https://github.com/icbi-lab/Immunophenogram)             |                                                  free (BSD)                                                   | [28052254](https://pubmed.ncbi.nlm.nih.gov/28052254/) | × |
|             [EPIC](https://gfellerlab.shinyapps.io/EPIC_1-1/)              | free for non-commercial use only ([Academic License](https://github.com/GfellerLab/EPIC/blob/master/LICENSE)) | [29130882](https://pubmed.ncbi.nlm.nih.gov/29130882/) | √ |
|           [ESTIMATE](https://gfellerlab.shinyapps.io/EPIC_1-1/)            |               free ([GPL 2.0](https://bioinformatics.mdanderson.org/public-software/estimate/))               | [24113773](https://pubmed.ncbi.nlm.nih.gov/24113773/) | × |
|              [ABIS](https://giannimonaco.shinyapps.io/ABIS/)               |                            free ([GPL 2.0](https://github.com/giannimonaco/ABIS))                             | [30726743](https://pubmed.ncbi.nlm.nih.gov/30726743/) | √ |
| [ConsensusTME](https://olliecast.shinyapps.io/Deconvolution_Benchmarking/) |              free ([GPL 3.0](https://github.com/cansysbio/ConsensusTME/blob/master/LICENSE.md))               | [31641033](https://pubmed.ncbi.nlm.nih.gov/31641033/) | × |
|       [quanTIseq](http://icbi.at/software/quantiseq/doc/index.html)        |                                                  free (BSD)                                                   | [31126321](https://pubmed.ncbi.nlm.nih.gov/31126321/) | √ |
:::

## Derive the Proportions of Different TME Cell Types

 You can use the function **deconv_TME()** to derive the proportions of different tumor microenvironment cell types from gene expression data with these tools.


```
devtools::install_github('dviraran/xCell')
devtools::install_github("GfellerLab/EPIC")
devtools::install_github("cansysbio/ConsensusTME")
devtools::install_github("federicomarini/quantiseqr")
```
::: {style="width:780px; height:200px; overflow-y: scroll; overflow-x: hidden;"}

```r
## TIMER
frac1 <- deconv_TME(MEL_GSE78220,method="TIMER")
```

```
## Found 130 genes with uniform expression within a single batch (all zeros); these will not be adjusted for batch.
```

```r
## CIBERSORT
frac2 <- deconv_TME(MEL_GSE78220,method="CIBERSORT")
```

```
## Called from: CIBERSORT(SE = SE, ...)
## debug: isList <- is.list(SE)
## debug: exp_mtr <- bind_mtr(SE, isList)
## debug: if (missing(sig_matrix)) {
##     LM22 <- NULL
##     data(LM22, package = "tigeR", envir = current_env())
##     sig_matrix <- LM22
##     result <- Ciber(sig_matrix, exp_mtr, perm, QN)
##     TME_data <- as.data.frame(result[, 1:22])
## } else {
##     result <- Ciber(sig_matrix, exp_mtr, perm, QN)
##     TME_data <- as.data.frame(result)
## }
## debug: LM22 <- NULL
## debug: data(LM22, package = "tigeR", envir = current_env())
## debug: sig_matrix <- LM22
## debug: result <- Ciber(sig_matrix, exp_mtr, perm, QN)
## debug: TME_data <- as.data.frame(result[, 1:22])
## debug: TME_data$group <- bind_meta(SE, isList)$response_NR
## debug: TME_data$sample <- rownames(TME_data)
## debug: TME_New <- reshape2::melt(TME_data)
## debug: colnames(TME_New) <- c("Group", "Sample", "Celltype", "Composition")
## debug: plot_order <- TME_New[TME_New$Group == "R", ] %>% dplyr::group_by(.data$Celltype) %>% 
##     dplyr::summarise(m = stats::median(.data$Composition)) %>% 
##     dplyr::arrange(dplyr::desc(.data$m)) %>% dplyr::pull(.data$Celltype)
## debug: TME_New$Celltype <- factor(TME_New$Celltype, levels = plot_order)
## debug: TME_New <- TME_New[!TME_New$Group == "UNK", ]
## debug: rs <- c()
## debug: warning_status <- 0
## debug: for (i in levels(TME_New$Celltype)) {
##     m <- TME_New[TME_New$Celltype == i, ]
##     R_S <- m[m$Group == "R", 4]
##     N_S <- m[m$Group == "N", 4]
##     if (any(R_S %in% N_S)) 
##         warning_status <- 1
##     rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## }
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: if (warning_status) message("There are identical relative abundance values in groups R and N for the '", 
##     i, "'. The p value of the Wilcoxon signed-rank test may not be precise due to ties in the data.")
## debug: message("There are identical relative abundance values in groups R and N for the '", 
##     i, "'. The p value of the Wilcoxon signed-rank test may not be precise due to ties in the data.")
## debug: selected_cells <- levels(TME_New$Celltype)
## debug: if (length(selected_cells) < 5) {
##     names(rs) <- seq_along(rs)
##     selected_cells <- levels(TME_New$Celltype)[as.numeric(names(sort(rs)[1:5]))]
## }
## debug: ciber_theme <- ggplot2::theme(plot.title = element_text(size = 12, 
##     color = "black", hjust = 0.5), axis.title = element_text(size = 10, 
##     color = "black"), axis.text = element_text(size = 10, color = "black"), 
##     axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top", 
##     legend.text = element_text(size = 12), legend.title = element_text(size = 12))
## debug: if (style == "raw") {
##     box_TME <- ggplot2::ggplot(TME_New[TME_New$Celltype %in% 
##         selected_cells, ], aes(x = .data$Celltype, y = .data$Composition)) + 
##         ggplot2::labs(y = "Cell composition", x = NULL, title = "TME Cell composition") + 
##         ggplot2::geom_boxplot(aes(fill = .data$Group), position = position_dodge(0.5), 
##             width = 0.5, outlier.alpha = 0) + ggplot2::scale_fill_manual(values = group_color) + 
##         ggplot2::theme_classic() + ciber_theme
##     y_max <- max(ggplot_build(box_TME)$data[[1]]$ymax)
##     box_TME <- box_TME + ggpubr::stat_compare_means(ggplot2::aes(group = .data$Group), 
##         label = "p.signif", method = "wilcox.test", hide.ns = TRUE, 
##         label.y = y_max * 1.1) + coord_cartesian(ylim = c(0, 
##         y_max * 1.1))
## }
## debug: if (style == "elegant") {
##     box_TME <- ggplot(TME_New, aes(x = .data$Celltype, y = .data$Composition, 
##         fill = .data$Group)) + stat_boxplot(data = TME_New, geom = "errorbar", 
##         width = 1, color = "black", linetype = "solid", position = position_dodge(0.8), 
##         linewidth = 0.7) + stat_boxplot(geom = "boxplot", color = "black", 
##         linetype = "solid", position = position_dodge(0.8), linewidth = 0.7, 
##         width = 0.8, outlier.shape = 19) + ggpubr::theme_classic2() + 
##         ciber_theme + scale_fill_manual(values = group_color)
##     y_max <- max(ggplot_build(box_TME)$data[[1]]$ymax)
##     box_TME <- box_TME + ggpubr::stat_compare_means(ggplot2::aes(group = .data$Group), 
##         label = "p.signif", method = "wilcox.test", hide.ns = FALSE, 
##         label.y = y_max * 1.1, size = 3) + coord_cartesian(ylim = c(0, 
##         y_max * 1.1))
## }
## debug: box_TME <- ggplot(TME_New, aes(x = .data$Celltype, y = .data$Composition, 
##     fill = .data$Group)) + stat_boxplot(data = TME_New, geom = "errorbar", 
##     width = 1, color = "black", linetype = "solid", position = position_dodge(0.8), 
##     linewidth = 0.7) + stat_boxplot(geom = "boxplot", color = "black", 
##     linetype = "solid", position = position_dodge(0.8), linewidth = 0.7, 
##     width = 0.8, outlier.shape = 19) + ggpubr::theme_classic2() + 
##     ciber_theme + scale_fill_manual(values = group_color)
## debug: y_max <- max(ggplot_build(box_TME)$data[[1]]$ymax)
## debug: box_TME <- box_TME + ggpubr::stat_compare_means(ggplot2::aes(group = .data$Group), 
##     label = "p.signif", method = "wilcox.test", hide.ns = FALSE, 
##     label.y = y_max * 1.1, size = 3) + coord_cartesian(ylim = c(0, 
##     y_max * 1.1))
## debug: list(t(result), box_TME)
```

```r
## MCPCounter
frac3 <- deconv_TME(MEL_GSE78220,method="MCPCounter")

## xCell
frac4 <- deconv_TME(MEL_GSE78220,method="xCell")
```

```
## [1] "Num. of genes: 10808"
## Setting parallel calculations through a MulticoreParam back-end
## with workers=4 and tasks=100.
## Estimating ssGSEA scores for 489 gene sets.
## [1] "Calculating ranks..."
## [1] "Calculating absolute values from ranks..."
##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   4%  |                                                                              |=====                                                                 |   7%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  18%  |                                                                              |===============                                                       |  21%  |                                                                              |==================                                                    |  25%  |                                                                              |====================                                                  |  29%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  39%  |                                                                              |==============================                                        |  43%  |                                                                              |================================                                      |  46%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  54%  |                                                                              |========================================                              |  57%  |                                                                              |==========================================                            |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  75%  |                                                                              |=======================================================               |  79%  |                                                                              |==========================================================            |  82%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  96%  |                                                                              |======================================================================| 100%
```

```r
## IPS
frac5 <- deconv_TME(MEL_GSE78220,method="IPS")

## EPIC
frac6 <- deconv_TME(MEL_GSE78220,method="epic")

## ESTIMATE
frac7 <- deconv_TME(MEL_GSE78220,method="ESTIMATE")
```

```
## [1] "Merged dataset includes 10130 genes (282 mismatched)."
## [1] "1 gene set: StromalSignature  overlap= 137"
## [1] "2 gene set: ImmuneSignature  overlap= 141"
```

```r
## ABIS
frac8 <- deconv_TME(MEL_GSE78220,method="ABIS")

## ConsensusTME
frac9 <- deconv_TME(MEL_GSE78220,method="ConsensusTME")
```

```
## Producing ConsensusTME Estimates Using The Following Parameters:
##  Statistical Framework: "ssgsea"
##  Gene Sets For Cancer Type: "SKCM"
##  Sample Size: 28
## Estimating ssGSEA scores for 19 gene sets.
## [1] "Calculating ranks..."
## [1] "Calculating absolute values from ranks..."
##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   4%  |                                                                              |=====                                                                 |   7%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  18%  |                                                                              |===============                                                       |  21%  |                                                                              |==================                                                    |  25%  |                                                                              |====================                                                  |  29%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  39%  |                                                                              |==============================                                        |  43%  |                                                                              |================================                                      |  46%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  54%  |                                                                              |========================================                              |  57%  |                                                                              |==========================================                            |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  75%  |                                                                              |=======================================================               |  79%  |                                                                              |==========================================================            |  82%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  96%  |                                                                              |======================================================================| 100%
## 
## [1] "Normalizing..."
```

```r
## quanTIseq
frac10 <- deconv_TME(MEL_GSE78220,method="quanTIseq")
```
:::

## Visualization and Comparing the Cell Proportions


```r
cell1 <- c("T cells CD4","Neutrophil", "Macrophage","mDCs","B cells", "T cells CD8")
pie1 <- fraction_pie(cell_name_filter(frac1),feature=factor(cell1, levels = cell1))

cell2 <- c("DCs resting", "T cells CD8", "T cells CD4 naive", "Macrophages M2", "Yd T cells", "Monocytes","Mast cells resting", "Neutrophils", "Tregs","B cells naive")
pie2 <- fraction_pie(cell_name_filter(frac2[[1]][1:22,]),feature=factor(cell2, levels = cell2))
pie1
```

<img src="03-TME-Analysis_files/figure-html/unnamed-chunk-3-1.png" width="2496" />

```r
pie2
```

<img src="03-TME-Analysis_files/figure-html/unnamed-chunk-3-2.png" width="2496" />

## Searching for Key Cell Types Associated with Immunotherapy Response


```r
## TIMER
TM <- deconv_TME(MEL_GSE91061,method = "TIMER")
```

```
## Found 125 genes with uniform expression within a single batch (all zeros); these will not be adjusted for batch.
```

```r
TM_SE <- SummarizedExperiment(assays=SimpleList(TM),
                               colData=colData(MEL_GSE91061))
browse_biomk(SE=TM_SE)
```

<img src="03-TME-Analysis_files/figure-html/unnamed-chunk-4-1.png" width="672" />

## Refining a Custom Reference Matrix for TME Deconvolution from Single-cell Sequencing Data

```r
library(Seurat)
```

```
## Loading required package: SeuratObject
```

```
## Loading required package: sp
```

```
## 
## Attaching package: 'sp'
```

```
## The following object is masked from 'package:IRanges':
## 
##     %over%
```

```
## 
## Attaching package: 'SeuratObject'
```

```
## The following object is masked from 'package:SummarizedExperiment':
## 
##     Assays
```

```
## The following object is masked from 'package:GenomicRanges':
## 
##     intersect
```

```
## The following object is masked from 'package:GenomeInfoDb':
## 
##     intersect
```

```
## The following object is masked from 'package:IRanges':
## 
##     intersect
```

```
## The following object is masked from 'package:S4Vectors':
## 
##     intersect
```

```
## The following object is masked from 'package:BiocGenerics':
## 
##     intersect
```

```
## The following object is masked from 'package:base':
## 
##     intersect
```

```
## 
## Attaching package: 'Seurat'
```

```
## The following object is masked from 'package:SummarizedExperiment':
## 
##     Assays
```

```r
library(magrittr)
```

```
## 
## Attaching package: 'magrittr'
```

```
## The following object is masked from 'package:GenomicRanges':
## 
##     subtract
```

```r
library(GenomicFeatures)
```

```
## Warning: package 'GenomicFeatures' was built under R version 4.3.3
```

```
## Loading required package: AnnotationDbi
```

```r
pbmc <- readRDS(system.file("extdata","pbmc.rds",
                            package = "tigeR",
                            mustWork = TRUE))
assay_Seu <- GetAssayData(pbmc, layer = "count")
count <- as.matrix(assay_Seu)
tpm <- count2tpm(count)
tpm_assay <- CreateAssayObject(counts = tpm)
pbmc[["TPM"]] <- tpm_assay
DefaultAssay(pbmc) <- "TPM"
pbmc[["RNA"]] <- NULL


pbmc <- FindVariableFeatures(pbmc, 
                             selection.method = "vst", 
                             nfeatures = 2000, verbose = FALSE)
pbmc <- NormalizeData(pbmc, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000, verbose = FALSE)
pbmc <- ScaleData(pbmc, features =  rownames(pbmc), verbose = FALSE)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE)
pbmc <- FindNeighbors(pbmc, dims = 1:10, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution = 0.5, verbose = FALSE)
new.cluster.ids <- c("Naive_CD4T", "CD14_Mono", "Memory_CD4T", "B_cell", "CD8T", "FCGR3A_Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc$seurat_clusters)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
PCAPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
```

<img src="03-TME-Analysis_files/figure-html/unnamed-chunk-5-1.png" width="672" />

```r
pbmc$celltype <- Idents(pbmc)
ref_sig_mtr <- refine_Reference(pbmc)
```

```
## Calculating cluster Naive_CD4T
```

```
## For a (much!) faster implementation of the Wilcoxon Rank Sum Test,
## (default method for FindMarkers) please install the presto package
## --------------------------------------------
## install.packages('devtools')
## devtools::install_github('immunogenomics/presto')
## --------------------------------------------
## After installation of presto, Seurat will automatically use the more 
## efficient implementation (no further action necessary).
## This message will be shown once per session
```

```
## Calculating cluster CD14_Mono
```

```
## Calculating cluster Memory_CD4T
```

```
## Calculating cluster B_cell
```

```
## Calculating cluster CD8T
```

```
## Calculating cluster FCGR3A_Mono
```

```
## Calculating cluster NK
```

```
## Calculating cluster DC
```

```
## Calculating cluster Platelet
```

```
## The following functions and any applicable methods accept the dots: CreateSeuratObject
```

```
## As of Seurat v5, we recommend using AggregateExpression to perform pseudo-bulk analysis.
## Names of identity class contain underscores ('_'), replacing with dashes ('-')
## This message is displayed once per session.
```

```r
cell_fraction <- deconv_TME(MEL_GSE100797,method = "CIBERSORT", sig_matrix=ref_sig_mtr)
```

```
## Called from: CIBERSORT(SE = SE, ...)
## debug: isList <- is.list(SE)
## debug: exp_mtr <- bind_mtr(SE, isList)
## debug: if (missing(sig_matrix)) {
##     LM22 <- NULL
##     data(LM22, package = "tigeR", envir = current_env())
##     sig_matrix <- LM22
##     result <- Ciber(sig_matrix, exp_mtr, perm, QN)
##     TME_data <- as.data.frame(result[, 1:22])
## } else {
##     result <- Ciber(sig_matrix, exp_mtr, perm, QN)
##     TME_data <- as.data.frame(result)
## }
## debug: result <- Ciber(sig_matrix, exp_mtr, perm, QN)
## debug: TME_data <- as.data.frame(result)
## debug: TME_data$group <- bind_meta(SE, isList)$response_NR
## debug: TME_data$sample <- rownames(TME_data)
## debug: TME_New <- reshape2::melt(TME_data)
```

```
## Using group, sample as id variables
```

```
## debug: colnames(TME_New) <- c("Group", "Sample", "Celltype", "Composition")
## debug: plot_order <- TME_New[TME_New$Group == "R", ] %>% dplyr::group_by(.data$Celltype) %>% 
##     dplyr::summarise(m = stats::median(.data$Composition)) %>% 
##     dplyr::arrange(dplyr::desc(.data$m)) %>% dplyr::pull(.data$Celltype)
## debug: TME_New$Celltype <- factor(TME_New$Celltype, levels = plot_order)
## debug: TME_New <- TME_New[!TME_New$Group == "UNK", ]
## debug: rs <- c()
## debug: warning_status <- 0
## debug: for (i in levels(TME_New$Celltype)) {
##     m <- TME_New[TME_New$Celltype == i, ]
##     R_S <- m[m$Group == "R", 4]
##     N_S <- m[m$Group == "N", 4]
##     if (any(R_S %in% N_S)) 
##         warning_status <- 1
##     rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## }
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: m <- TME_New[TME_New$Celltype == i, ]
## debug: R_S <- m[m$Group == "R", 4]
## debug: N_S <- m[m$Group == "N", 4]
## debug: if (any(R_S %in% N_S)) warning_status <- 1
## debug: warning_status <- 1
## debug: rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
## debug: if (warning_status) message("There are identical relative abundance values in groups R and N for the '", 
##     i, "'. The p value of the Wilcoxon signed-rank test may not be precise due to ties in the data.")
## debug: message("There are identical relative abundance values in groups R and N for the '", 
##     i, "'. The p value of the Wilcoxon signed-rank test may not be precise due to ties in the data.")
```

```
## There are identical relative abundance values in groups R and N for the 'Platelet'. The p value of the Wilcoxon signed-rank test may not be precise due to ties in the data.
```

```
## debug: selected_cells <- levels(TME_New$Celltype)
## debug: if (length(selected_cells) < 5) {
##     names(rs) <- seq_along(rs)
##     selected_cells <- levels(TME_New$Celltype)[as.numeric(names(sort(rs)[1:5]))]
## }
## debug: ciber_theme <- ggplot2::theme(plot.title = element_text(size = 12, 
##     color = "black", hjust = 0.5), axis.title = element_text(size = 10, 
##     color = "black"), axis.text = element_text(size = 10, color = "black"), 
##     axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top", 
##     legend.text = element_text(size = 12), legend.title = element_text(size = 12))
## debug: if (style == "raw") {
##     box_TME <- ggplot2::ggplot(TME_New[TME_New$Celltype %in% 
##         selected_cells, ], aes(x = .data$Celltype, y = .data$Composition)) + 
##         ggplot2::labs(y = "Cell composition", x = NULL, title = "TME Cell composition") + 
##         ggplot2::geom_boxplot(aes(fill = .data$Group), position = position_dodge(0.5), 
##             width = 0.5, outlier.alpha = 0) + ggplot2::scale_fill_manual(values = group_color) + 
##         ggplot2::theme_classic() + ciber_theme
##     y_max <- max(ggplot_build(box_TME)$data[[1]]$ymax)
##     box_TME <- box_TME + ggpubr::stat_compare_means(ggplot2::aes(group = .data$Group), 
##         label = "p.signif", method = "wilcox.test", hide.ns = TRUE, 
##         label.y = y_max * 1.1) + coord_cartesian(ylim = c(0, 
##         y_max * 1.1))
## }
## debug: if (style == "elegant") {
##     box_TME <- ggplot(TME_New, aes(x = .data$Celltype, y = .data$Composition, 
##         fill = .data$Group)) + stat_boxplot(data = TME_New, geom = "errorbar", 
##         width = 1, color = "black", linetype = "solid", position = position_dodge(0.8), 
##         linewidth = 0.7) + stat_boxplot(geom = "boxplot", color = "black", 
##         linetype = "solid", position = position_dodge(0.8), linewidth = 0.7, 
##         width = 0.8, outlier.shape = 19) + ggpubr::theme_classic2() + 
##         ciber_theme + scale_fill_manual(values = group_color)
##     y_max <- max(ggplot_build(box_TME)$data[[1]]$ymax)
##     box_TME <- box_TME + ggpubr::stat_compare_means(ggplot2::aes(group = .data$Group), 
##         label = "p.signif", method = "wilcox.test", hide.ns = FALSE, 
##         label.y = y_max * 1.1, size = 3) + coord_cartesian(ylim = c(0, 
##         y_max * 1.1))
## }
## debug: box_TME <- ggplot(TME_New, aes(x = .data$Celltype, y = .data$Composition, 
##     fill = .data$Group)) + stat_boxplot(data = TME_New, geom = "errorbar", 
##     width = 1, color = "black", linetype = "solid", position = position_dodge(0.8), 
##     linewidth = 0.7) + stat_boxplot(geom = "boxplot", color = "black", 
##     linetype = "solid", position = position_dodge(0.8), linewidth = 0.7, 
##     width = 0.8, outlier.shape = 19) + ggpubr::theme_classic2() + 
##     ciber_theme + scale_fill_manual(values = group_color)
## debug: y_max <- max(ggplot_build(box_TME)$data[[1]]$ymax)
## debug: box_TME <- box_TME + ggpubr::stat_compare_means(ggplot2::aes(group = .data$Group), 
##     label = "p.signif", method = "wilcox.test", hide.ns = FALSE, 
##     label.y = y_max * 1.1, size = 3) + coord_cartesian(ylim = c(0, 
##     y_max * 1.1))
## debug: list(t(result), box_TME)
```

## 📝 More Details About TME Analysis {.unnumbered}

 **TIMER** is a comprehensive tool for systematical analysis of immune infiltrates across diverse cancer types.

 **CIBERSORT** is an analytical tool from the Alizadeh Lab and Newman Lab to impute gene expression profiles and provide an estimation of the abundances of member cell types in a mixed cell population, using gene expression data.

 **MCPCounter** is called the Microenvironment Cell Populations-counter. It allows the robust quantification of the absolute abundance of eight immune and two stromal cell populations in heterogeneous tissues from transcriptomic data

 **xCell** is a gene signature-based method learned from thousands of pure cell types from various sources. xCell applies a novel technique for reducing associations between closely related cell types.

 **IPS** uses an analytical strategy to provide comprehensive view of 28 TIL subpopulations including effector and memory T cells and immunosuppressive cells (Tregs, MDSCs).

 **EPIC** is called Estimating the Proportion of Immune and Cancer cells. It compares the level of expression of genes in a tumor with a library of the gene expression profiles from specific cell types that can be found in tumors and uses this information to predict how many of each type of cell are present.

 **ESTIMATE** is described as ‘Estimation of STromal and Immune cells in MAlignant Tumours using Expression data'. It is a method that utilizes gene expression signatures to infer the proportion of stromal and immune cells within tumor samples.

 **ABIS** is called ABsolute Immune Signal (ABIS) deconvolution. ABIS performs absolute deconvolution on RNA-Seq and microarray data.
 
 **Consensus**<sup>TME</sup> uses a consensus approch to generate cancer-specific signatres for TME cell types, and utilizes ssGSEA framework to estimate the relative abundance of these cell types.

 **quanTIseq** is a computational pipeline for the quantification of the Tumor Immune contexture from human RNA-seq data.
