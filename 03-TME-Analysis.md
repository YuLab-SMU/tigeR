# 🌱 Tumor Microenviroment Analysis

## Availiable TME Analysis Method in tigeR

10 open source tumor microenvironment deconvolution methods are indlcued in tigeR.

::: {style="width:780px; height:200px; overflow-y: scroll; overflow-x: hidden;"}
|                                 Algorithm                                  |                                                    license                                                    |                         PMID                          |
|:--------------------:|:------------------------------:|:----------------:|
|                    [TIMER](http://cistrome.org/TIMER/)                     |                                                free (GPL 2.0)                                                 | [27549193](https://pubmed.ncbi.nlm.nih.gov/27549193/) |
|                [CIBERSORT](https://cibersort.stanford.edu/)                |                                       free for non-commercial use only                                        | [25822800](https://pubmed.ncbi.nlm.nih.gov/25822800/) |
|             [MCPCounter](https://github.com/ebecht/MCPcounter)             |               free ([GPL 3.0](https://github.com/ebecht/MCPcounter/blob/master/Source/License))               | [27765066](https://pubmed.ncbi.nlm.nih.gov/27765066/) |
|                      [xCell](http://xcell.ucsf.edu/)                       |                  free ([GPL 3.0](https://github.com/dviraran/xCell/blob/master/DESCRIPTION))                  | [29141660](https://pubmed.ncbi.nlm.nih.gov/29141660/) |
|             [IPS](https://github.com/icbi-lab/Immunophenogram)             |                                                  free (BSD)                                                   | [28052254](https://pubmed.ncbi.nlm.nih.gov/28052254/) |
|             [EPIC](https://gfellerlab.shinyapps.io/EPIC_1-1/)              | free for non-commercial use only ([Academic License](https://github.com/GfellerLab/EPIC/blob/master/LICENSE)) | [29130882](https://pubmed.ncbi.nlm.nih.gov/29130882/) |
|           [ESTIMATE](https://gfellerlab.shinyapps.io/EPIC_1-1/)            |               free ([GPL 2.0](https://bioinformatics.mdanderson.org/public-software/estimate/))               | [24113773](https://pubmed.ncbi.nlm.nih.gov/24113773/) |
|              [ABIS](https://giannimonaco.shinyapps.io/ABIS/)               |                            free ([GPL 2.0](https://github.com/giannimonaco/ABIS))                             | [30726743](https://pubmed.ncbi.nlm.nih.gov/30726743/) |
| [ConsensusTME](https://olliecast.shinyapps.io/Deconvolution_Benchmarking/) |              free ([GPL 3.0](https://github.com/cansysbio/ConsensusTME/blob/master/LICENSE.md))               | [31641033](https://pubmed.ncbi.nlm.nih.gov/31641033/) |
|       [quanTIseq](http://icbi.at/software/quantiseq/doc/index.html)        |                                                  free (BSD)                                                   | [31126321](https://pubmed.ncbi.nlm.nih.gov/31126321/) |
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
```

```
## Called from: fraction_pie(cell_name_filter(frac1), feature = factor(cell1, 
##     levels = cell1))
## debug: ggplot() + scatterpie::geom_scatterpie(data = mt, aes(x, y, r = radius), 
##     pie_scale = pie_scale, cols = colnames(mt)[1:nrow(mtr)], 
##     color = color, label_radius = label_radius, label_show_ratio = label_show_ratio, 
##     label_threshold = label_threshold, fontsize = fontsize, ...) + 
##     scale_fill_manual(values = colors) + theme_void() + theme(legend.position = c(0.5, 
##     1.2), legend.title = element_blank(), legend.text = element_text(size = 18), 
##     plot.background = element_rect(fill = "transparent", color = "transparent"), 
##     panel.background = element_rect(fill = "transparent", color = "transparent")) + 
##     coord_map() + guides(fill = guide_legend(nrow = 1))
```

```r
cell2 <- c("DCs resting", "T cells CD8", "T cells CD4 naive", "Macrophages M2", "Yd T cells", "Monocytes","Mast cells resting", "Neutrophils", "Tregs","B cells naive")
pie2 <- fraction_pie(cell_name_filter(frac2[[1]][1:22,]),feature=factor(cell2, levels = cell2))
```

```
## Called from: fraction_pie(cell_name_filter(frac2[[1]][1:22, ]), feature = factor(cell2, 
##     levels = cell2))
## debug: ggplot() + scatterpie::geom_scatterpie(data = mt, aes(x, y, r = radius), 
##     pie_scale = pie_scale, cols = colnames(mt)[1:nrow(mtr)], 
##     color = color, label_radius = label_radius, label_show_ratio = label_show_ratio, 
##     label_threshold = label_threshold, fontsize = fontsize, ...) + 
##     scale_fill_manual(values = colors) + theme_void() + theme(legend.position = c(0.5, 
##     1.2), legend.title = element_blank(), legend.text = element_text(size = 18), 
##     plot.background = element_rect(fill = "transparent", color = "transparent"), 
##     panel.background = element_rect(fill = "transparent", color = "transparent")) + 
##     coord_map() + guides(fill = guide_legend(nrow = 1))
```

```r
pie1
```

```
## Warning: Removed 3 rows containing missing values or values outside the scale range
## (`geom_text()`).
```

<img src="03-TME-Analysis_files/figure-html/unnamed-chunk-3-1.png" width="2496" />

```r
pie2
```

```
## Warning: Removed 146 rows containing missing values or values outside the scale range
## (`geom_text()`).
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
