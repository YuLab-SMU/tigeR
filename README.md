# tigeR
`tigeR` is an R package designed for the analysis of gene expression in tumor immunotherapy.

<p align="center">
    <img src="https://raw.githubusercontent.com/Chengxugorilla/tigeR.extra/main/logo.png">
</p>

<p align="center">
    <img src="https://raw.githubusercontent.com/Chengxugorilla/tigeR.extra/main/logo.pdf">
</p>

## Requirements
`install.packages(c("devtools", "ggplot2", "pROC"))`

## Install
devtools::install_github("Chengxugorilla/tigeR")

## Quick Start
### 1. Load packages and demo data  
```
library(tigeR)
Dataloader(pick=c(4,5,13,14,18))
```
`pick` a number(1-20) or a numeric vector specifying the corresponding dataset(s) you wish to load. Alternatively, you can use Dataloader() with pick=NULL to get an overview of all available datasets.

When the user enters a number between 1 and 20, this function will load the corresponding dataset into the current_env(). If pick is NULL (is.null(pick) == TRUE), Dataloader() will return a data.frame containing an overview of all the datasets.

### 2. Biomarkder
```
integrate_analysis(SE=MEL_GSE91061, geneSet="CD274")

```
`SE` the dataset you wish to use to perform Differential Expression Analysis and survival analysis. A SummarizedExperiment (SE) object, which can be either a single SE object or a list of SE objects. Note that for each SE object, the colData must contain treatment information under the column name Treatment.

`geneSet` the gene you are interested in.

`Immunotherapy_Response()` will return a list, with the following elements: the first element represents the result of the differential expression analysis between Responder and Non-Responder, the second element represents the result of the differential expression analysis between Pre-Treatment and Post-Treatment, and the third element represents the result of the survival analysis.

Firstly, you need to library ggplot2.
`library(ggplot2)`

You can use diff_biomk to visualize your analysis.

```
diff_biomk(SE=MEL_PRJEB23709,gene='CXCL13',type='Response')  +
  ggtitle("MEL-GSE91061") +
  theme(legend.position = "none")
```
<p align="center">
    <img src="https://raw.githubusercontent.com/Chengxugorilla/tigeR.extra/main/Response.svg" alt="Screenshot">
</p>

```
idx_CTLA <- MEL_GSE115821$Therapy=="anti-PD-1"
diff_biomk(MEL_GSE115821[,MEL_GSE115821$Therapy=="anti-PD-1"],
             gene = "CXCL13",type = "Treatment",p.round=3,
             log_sc = TRUE,p.pos = c(0.05,0.60),textcol="black") +
  ylim(0,6.2) +
  ggtitle("MEL-GSE115821") +
  theme(legend.position = "none")
```
<p align="center">
    <img src="https://raw.githubusercontent.com/Chengxugorilla/tigeR.extra/main/Treatment.svg" alt="Screenshot">
</p>

`SE` the data set or data sets.

`gene` the gene you interest in.

`type` the type of analysis you want to perform, which could be either ‘Treatment’ or ‘Response’. This determines whether you want to compare Responder vs Non-Responder or Pre-Treatment vs Post-Treatment.”


You can also visualization survival analysis using plt_surv() function.
```
surv_biomk(MEL_PRJEB23709,gene = "CXCL13",lg.pos=c(0.8,0.92),
             val.pos = c(0,0.2),lg.text = "specific")$plot +
  theme(plot.margin = unit(c(3, 1, 1, 1), "lines"),
        legend.key.height = unit(0,"cm"),
        legend.key.spacing.y = unit(0,"cm"),
        legend.key.size = unit(0,"cm")) +
  ggtitle("MEL-PRJEB23709")
```

<p align="center">
    <img src="https://raw.githubusercontent.com/Chengxugorilla/tigeR.extra/main/Survival.svg" alt="Screenshot">
</p>

### 3. Calculate signature scores of existing immunotherapy  and  Assess Signature using existing data

```
Sig_scores <- Signature_calculation(SE=MEL_GSE78220)
```
`SE` a SummarizedExperiment object for which you want to calculate the Signature Score.

&emsp;By employing the Signature_calculation() function, you can obtain a comprehensive signature score matrix for the 23 signatures in TigeR. In this matrix, the columns represent the signature scores, and the rows denote the sample names.

```
result <- Signature_assessment(MEL_PRJEB23709,
                               Weighted_mean_Sigs$Tcell_inflamed_GEP,
                               rmBE=TRUE,
                               response_NR=TRUE)
result[[1]]
result[[2]]
```
`SE` the dataset you wish to use to test your Signature. A SummarizedExperiment (SE) object, which can be either a single SE object or a list of SE objects. Note that for each SE object, the colData must contain treatment information under the column name Treatment.

`Signature` the gene set which you want to assess.

`rmBE` whether remove batch effect between different data set using internal Combat method.

`response_NR`	a logical variable.If TRUE, the function will automatically convert the patient's drug response (such as PR, NR, SD, etc. to binary value NR (non-responder) or R (Responder)).


By employing the `Signature_assessment()` function, you can assess the performance of Signature(including user-built Signature) for response prediction in different datasets. The function will return a "roc" object, a list of class "roc".

<p align="center">
    <img src="https://raw.githubusercontent.com/Chengxugorilla/tigeR.extra/main/Sig_ROC.png" alt="Screenshot">
</p>


### 4. Tumor Microenvironment Deconvolution
tigeR integrates 10 open-source TME deconvolution method, namely CIBERSORT, TIMER, ESTIMATE, IPS, xCell, EPIC, ConsensusTME, ABIS, quanTIseq and MCPCounter.

|Algorithm |license |citation |
|-----------------------------------|--------------|---------------------------------|
|[TIMER](http://cistrome.org/TIMER/)|free (GPL 2.0)|Li, B., Severson, E., Pignon, J.-C., Zhao, H., Li, T., Novak, J., … Liu, X. S. (2016). Comprehensive analyses of tumor immunity: implications for cancer immunotherapy. Genome Biology, 17(1), 174. https://doi.org/10.1186/s13059-016-1028-7|
|[CIBERSORT](https://cibersort.stanford.edu/)|free for non-commerical use only|Newman, A. M., Liu, C. L., Green, M. R., Gentles, A. J., Feng, W., Xu, Y., … Alizadeh, A. A. (2015). Robust enumeration of cell subsets from tissue expression profiles. Nature Methods, 12(5), 453–457. https://doi.org/10.1038/nmeth.3337|
|[MCPCounter](https://github.com/ebecht/MCPcounter) |free ([GPL 3.0](https://github.com/ebecht/MCPcounter/blob/master/Source/License))|Becht, E., Giraldo, N. A., Lacroix, L., Buttard, B., Elarouci, N., Petitprez, F., … de Reyniès, A. (2016). Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression. Genome Biology, 17(1), 218. https://doi.org/10.1186/s13059-016-1070-5|
|[xCell](http://xcell.ucsf.edu/)|free ([GPL 3.0](https://github.com/dviraran/xCell/blob/master/DESCRIPTION))|Aran, D., Hu, Z., & Butte, A. J. (2017). xCell: digitally portraying the tissue cellular heterogeneity landscape. Genome Biology, 18(1), 220. https://doi.org/10.1186/s13059-017-1349-1|
|[IPS](https://github.com/icbi-lab/Immunophenogram)| free (BSD)|P. Charoentong et al., Pan-cancer Immunogenomic Analyses Reveal Genotype-Immunophenotype Relationships and Predictors of Response to Checkpoint Blockade. Cell Reports 18, 248-262 (2017). https://doi.org/10.1016/j.celrep.2016.12.019|
|[EPIC](https://gfellerlab.shinyapps.io/EPIC_1-1/)|free for non-commercial use only ([Academic License](https://github.com/GfellerLab/EPIC/blob/master/LICENSE))|Racle, J., de Jonge, K., Baumgaertner, P., Speiser, D. E., & Gfeller, D. (2017). Simultaneous enumeration of cancer and immune cell types from bulk tumor gene expression data. ELife, 6, e26476. https://doi.org/10.7554/eLife.26476|
|[ESTIMATE](https://gfellerlab.shinyapps.io/EPIC_1-1/) | free ([GPL 2.0](https://bioinformatics.mdanderson.org/public-software/estimate/))|Yoshihara, K., Shahmoradgoli, M., Martínez, E., Vegesna, R., Kim, H., Torres-Garcia, W., Treviño, V., Shen, H., Laird, P. W., Levine, D. A., Carter, S. L., Getz, G., Stemke-Hale, K., Mills, G. B., & Verhaak, R. G. (2013). Inferring tumour purity and stromal and immune cell admixture from expression data. Nature communications, 4, 2612. https://doi.org/10.1038/ncomms3612|
|[ABIS](https://giannimonaco.shinyapps.io/ABIS/) |free ([GPL 2.0](https://github.com/giannimonaco/ABIS)) |Monaco, G., Lee, B., Xu, W., Mustafah, S., Hwang, Y. Y., ..., Larbi, A. (2019). RNA-Seq Signatures Normalized by mRNA Abundance Allow Absolute Deconvolution of Human Immune Cell Types. Cell reports, 26(6), 1627–1640.e7. https://doi.org/10.1016/j.celrep.2019.01.041|
|[ConsensusTME](https://olliecast.shinyapps.io/Deconvolution_Benchmarking/)|free ([GPL 3.0](https://github.com/cansysbio/ConsensusTME/blob/master/LICENSE.md))|Jiménez-Sánchez, A., Cast, O., & Miller, M. L. (2019). Comprehensive Benchmarking and Integration of Tumor Microenvironment Cell Estimation Methods. Cancer research, 79(24), 6238–6246. https://doi.org/10.1158/0008-5472.CAN-18-3560|
|[quanTIseq](http://icbi.at/software/quantiseq/doc/index.html)|free (BSD)|Finotello, F., Mayer, C., Plattner, C., Laschober, G., Rieder, D., Hackl, H., ..., Sopper, S. (2019). Molecular and pharmacological modulators of the tumor immune contexture revealed by deconvolution of RNA-seq data. Genome medicine, 11(1), 34. https://doi.org/10.1186/s13073-019-0638-6|


```
result <- CIBERSORT(sig_matrix=LM22,SE=MEL_GSE91061,perm=10, QN=T)
result[[2]]
```
`sig_matrix` gene expression matrix from isolated cells.

`SE` the bulk RNA-seq dataset that you want to use for deconvolution and obtaining its cell fraction.

`perm` the number of permutations.

`PT` whether perform quantile normalization or not (TRUE/FALSE).

`CIBERSORT()` function will return a list with the following elements: the first element is a matrix representing the cell fraction of each sample, and the second element is a ggplot object that visualizes the difference in cell fraction between Responders and Non-Responders.”

<p align="center">
    <img src="https://raw.githubusercontent.com/Chengxugorilla/tigeR.extra/main/CIBERSORT.png" alt="Screenshot">
</p>

### 5. Build machine learning model for immunotherapy prognosis prediction
#### build model
```
train_set <- list(MEL_GSE91061, MEL_phs000452, RCC_Braun_2020)
mymodel <- build_Model(Model='NB', SE=train_set, feature_genes=Stem.Sig, response_NR = TRUE)
```
`Model` represents the type of model you want to build. You have several options to choose from: "NB" for Naive Bayes, "SVM" for Support Vector Machine, "RF" for Random Forest, "CC" for Cancerclass, "ADB" for Adaboost, "LGB" for Logitboost, and "LGT" for Logistics.

`SE` the dataset you wish to use to build your model. A SummarizedExperiment (SE) object, which can be either a single SE object or a list of SE objects. Note that for each SE object, the colData must contain treatment information under the column name Treatment.

`feature_gene` refers to the specific set of genes you wish to use for model construction.

`rmBE` refers to the option of whether to remove batch effects between different datasets using the internal Combat method.

`response_NR` a logical variable.If TRUE, the function will automatically convert the patient's drug response (such as PR, NR, SD, etc. to binary value NR (non-responder) or R (Responder)).


In this case, `build_Model()` will return a "naiveBayes" object.
#### test model
```
test_set <- list(MEL_GSE78220, MEL_PRJEB23709)
Result <- test_Model(Model=mymodel, SE=test_set)

## ROC curve and AUC
Result[[2]]

## Predictive response
Result[[3]]
```
`Model` the model that you want to test, which is generated by the build_Model() function.

`SE` the dataset you wish to use to test your model. A SummarizedExperiment (SE) object, which can be either a single SE object or a list of SE objects. Note that for each SE object, the colData must contain treatment information under the column name Treatment.

`test_Model()` will return an "roc" object. You can use the plot() function to plot the ROC curve and the auc() function to calculate the Area Under the Curve (AUC) of the ROC.

<p align="center">
    <img src="https://raw.githubusercontent.com/Chengxugorilla/tigeR.extra/main/ROC.png" alt="Screenshot">
</p>


## TIGER web server
http://tiger.canceromics.org/#/
