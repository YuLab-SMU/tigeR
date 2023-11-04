# tigeR
`tigeR` is a R package for Tumor Immunotherapy Gene Expression Analysis
![Screenshot](https://github.com/Chengxugorilla/tigeR/raw/main/man/figures/logo.png)

## Requirements
`install.packages(c("devtools", "ggplot2", "pROC"))`

## Install
devtools::install_github("Chengxugorilla/tigeR")

## Quick Start
`tigeR` generally supports the quantification and visualization of tumor immunotherapy data. 

### 1. Load packages and demo data
The demo data include 4 Melanoma RNA-seq datasets and 1 Renal Cell Carcinoma dataset. Use `Dataloader()` to get an overview of all datasets.
```
library(tigeR)
Dataloader(c(4,5,13,14,18))
```
### 2. Calculate signature scores of existing immunotherapy prognosis signatures
```
Sig_scores <- calculate_Signature_Score(assay(MEL_GSE78220))
View(Sig_scores)
```
### 3. Build machine learning model for immunotherapy prognosis prediction
  tigeR allows users to build machine learning prediction models for immunotherapy prognosis. There 7 model including Naive-Bayes Model, Random Forest Model, SVM Model, Cancerclass Model, Adaboost Model, Logitboost Model, Logistics Model.
  We use Naive-Bayes Model to make an example:
```
#building model
train_set <- list(MEL_GSE91061, MEL_phs000452, RCC_Braun_2020)
mymodel <- build_Model('NB', train_set, Stem.Sig, response_NR = TRUE)

#test model
test_set <- list(MEL_GSE78220, MEL_PRJEB23709)
ROC <- test_Model(mymodel, test_set)

#Drawing roc curve and calculating the AUC of roc curver.
library(pROC)
plot(ROC)
auc(ROC)
## Area under the curve: 0.5983
```
![Screenshot](https://github.com/Chengxugorilla/tigeR/raw/main/man/figures/ROC.png)
### 4. Immunotherapy Response
Perform differential expression analysis and survival analysis in certain gene and return the result.
```
Immunotherapy_Response(gene='CD274', MEL_GSE91061)
```
### 5. Visualization
You can use plt_diff and plt_surv to visualize your analysis.
```
plt_diff('CD274',MEL_GSE91061,'Treatment') +
  ggtitle("Treatment vs UnTreatment") +
  theme(plot.title = element_text(hjust = 0.5)) 
```
![Screenshot](https://github.com/Chengxugorilla/tigeR/raw/main/man/figures/Treatment.png)
```
plt_diff('CD274',MEL_GSE91061,'Response') +
  ggtitle("Responder vs Non-Responder") +
  theme(plot.title = element_text(hjust = 0.5))
```
![Screenshot](https://github.com/Chengxugorilla/tigeR/raw/main/man/figures/Response.png)
```
P <- plt_surv('CD274',MEL_GSE91061)
P$plot <- P$plot +
  ggtitle("Survival analysis") +
  theme(plot.title = element_text(hjust = 0.5))
P
```
![Screenshot](https://github.com/Chengxugorilla/tigeR/raw/main/man/figures/Survival.png)
### 6. Cibersort
Cibersort algorithm is embeded in tigeR package. Users can use `CIBERSORT()` function for cell fraction deconvolution.
```
data("LM22",package = "tigeR")

result <- CIBERSORT(LM22,MEL_GSE78220,perm=10, QN=T)
result[[2]]
```
![Screenshot](https://github.com/Chengxugorilla/tigeR/raw/main/man/figures/CIBERSORT.png)

## TIGER web server
http://tiger.canceromics.org/#/
