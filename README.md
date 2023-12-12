# tigeR
`tigeR` is an R package designed for the analysis of gene expression in tumor immunotherapy.
![Screenshot](https://github.com/Chengxugorilla/tigeR/raw/main/man/figures/logo.png)

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
`pick` a number(1-20) or a numeric vector specify the corresponding dataset(s) you wish to load. Alternatively, you can use Dataloader() with pick=NULL to get an overview of all available datasets.

When the user enters a number between 1 and 20, this function will load the corresponding dataset into the current_env(). If pick is NULL (is.null(pick) == TRUE), Dataloader() will return a data.frame containing an overview of all the datasets.

### 2. Calculate signature scores of existing immunotherapy 
```
Sig_scores <- calculate_Signature_Score(exp_mtr=assay(MEL_GSE78220))
```
`exp_mtr` an expression matrix for which you want to calculate the Signature Score.

&emsp;By employing the calculate_Signature_Score() function, you can obtain a comprehensive signature score matrix for the 23 signatures in TigeR. In this matrix, the columns represent the signature scores, and the rows denote the sample names.
### 3. Assess Signature using existing data
```
result <- Signature_assessment(MEL_PRJEB23709,
                               Weighted_mean_Sigs$Tcell_inflamed_GEP,
                               rmBE=TRUE,
                               response_NR=TRUE)
plot(result)
```
`SE` the dataset you wish to use to test your Signature. A SummarizedExperiment (SE) object, which can be either a single SE object or a list of SE objects. Note that for each SE object, the colData must contain treatment information under the column name Treatment.

`Signature` the gene set which you want to assess.

`rmBE` whether remove batch effect between different data set using internal Combat method.

`response_NR`	if TRUE, only use R or NR to represent Immunotherapy response of patients.

By employing the `Signature_assessment()` function, you can assess the performance of Signature(including user-built Signature) in different datasets. The function will return a "roc" object, a list of class "roc".
![Screenshot](https://github.com/Chengxugorilla/tigeR/raw/main/man/figures/Sig_ROC.png)

### 4. Build machine learning model for immunotherapy prognosis prediction
#### build model
```
train_set <- list(MEL_GSE91061, MEL_phs000452, RCC_Braun_2020)
mymodel <- build_Model(Model='NB', SE=train_set, feature_genes=Stem.Sig, response_NR = TRUE)
```
`Model` represents the type of model you want to build. You have several options to choose from: "NB" for Naive Bayes, "SVM" for Support Vector Machine, "RF" for Random Forest, "CC" for Cancerclass, "ADB" for Adaboost, "LGB" for Logitboost, and "LGT" for Logistics.

`SE` the dataset you wish to use to build your model. A SummarizedExperiment (SE) object, which can be either a single SE object or a list of SE objects. Note that for each SE object, the colData must contain treatment information under the column name Treatment.

`feature_gene` refers to the specific set of genes you wish to use for model construction.

`rmBE` refers to the option of whether to remove batch effects between different datasets using the internal Combat method.

`response_NR` if TRUE, only use R or NR to represent Immunotherapy response of patients.

In this case, `build_Model()` will return a "naiveBayes" object.
#### test model
```
test_set <- list(MEL_GSE78220, MEL_PRJEB23709)
ROC <- test_Model(Model=mymodel, SE=test_set)

library(pROC)
plot(ROC)
auc(ROC)
## Area under the curve: 0.5983
```
`Model` the model that you want to test, which is generated by the build_Model() function.

`SE` the dataset you wish to use to test your model. A SummarizedExperiment (SE) object, which can be either a single SE object or a list of SE objects. Note that for each SE object, the colData must contain treatment information under the column name Treatment.

`test_Model()` will return an "roc" object. You can use the plot() function to plot the ROC curve and the auc() function to calculate the Area Under the Curve (AUC) of the ROC.
![Screenshot](https://github.com/Chengxugorilla/tigeR/raw/main/man/figures/ROC.png)
### 5. Immunotherapy Response
```
Immunotherapy_Response(SE=MEL_GSE91061, gene="CD274")
```
`SE` the dataset you wish to use to perform DEA and survival analysis. A SummarizedExperiment (SE) object, which can be either a single SE object or a list of SE objects. Note that for each SE object, the colData must contain treatment information under the column name Treatment.

`gene` the gene you are interested in.

`Immunotherapy_Response()` will return a list, with the following elements: the first element represents the result of the differential expression analysis between Responder and Non-Responder, the second element represents the result of the differential expression analysis between Pre-Treatment and Post-Treatment, and the third element represents the result of the survival analysis.

### 6. Visualization
Firstly, you need to library ggplot2.
`library(ggplot2)`

You can use plt_diff and plt_surv to visualize your analysis.

```
plt_diff(SE=MEL_GSE91061,gene='CD274',type='Treatment') 
```
`SE` the data set or data sets.

`gene` the gene you interest in.

`type` the type of analysis you want to perform, which could be either ‘Treatment’ or ‘Response’. This determines whether you want to compare Responder vs Non-Responder or Pre-Treatment vs Post-Treatment.”
![Screenshot](https://github.com/Chengxugorilla/tigeR/raw/main/man/figures/Treatment.png)
```
plt_diff(SE=MEL_GSE91061,gene='CD274',type='Response') 
```
![Screenshot](https://github.com/Chengxugorilla/tigeR/raw/main/man/figures/Response.png)
You can also visualization survival analysis using plt_surv() function.
```
plt_surv(SE=MEL_GSE91061,gene=c('CD274','5S_rRNA'),method = "GSVA)
```
![Screenshot](https://github.com/Chengxugorilla/tigeR/raw/main/man/figures/Survival.png)
### 7. Cibersort
Cibersort algorithm is embeded in tigeR package. Users can use `CIBERSORT()` function for cell fraction deconvolution.
```
data("LM22",package = "tigeR")

result <- CIBERSORT(sig_matrix=LM22,SE=MEL_GSE78220,perm=10, QN=T)
result[[2]]
```
`sig_matrix` gene expression matrix from isolated cells.

`SE` the bulk RNA-seq dataset that you want to use for deconvolution and obtaining its cell fraction.

`perm` the number of permutations.

`PT` whether perform quantile normalization or not (TRUE/FALSE).

`CIBERSORT()` function will return a list with the following elements: the first element is a matrix representing the cell fraction of each sample, and the second element is a ggplot object that visualizes the difference in cell fraction between Responders and Non-Responders.”
![Screenshot](https://github.com/Chengxugorilla/tigeR/raw/main/man/figures/CIBERSORT.png)

## TIGER web server
http://tiger.canceromics.org/#/
