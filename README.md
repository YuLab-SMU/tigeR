## 1.Built-in Model
![Screenshot](https://github.com/Chengxugorilla/tigeR/raw/main/logo.png)

```
library(tigeR)
library(e1071)
library(pROC)
data(bayesmodel, package = "tigeR")
data(Stem.Sig, package = "tigeR")

test_Expr <- extract_mtr('MEL_GSE78220_exp')
test_Expr <- dataPreprocess(test_Expr, Stem.Sig)
test_response <- extract_label('MEL_GSE78220_meta')
predict_response_R <- predict(bayesmodel, t(test_Expr), type = 'class') == 'R'

#Obtaining the meta informations of patients whose prediction results are 'Response'.
data("MEL_GSE78220_meta")

rc <- MEL_GSE78220_meta[predict_response_R,]
rc$response <- sub('CR|MR|PR|SD|CRPR', 'R', rc$response)
rc$response <- sub('PD', 'NR', rc$response)
rc$response <- as.factor(rc$response)

#Drawing roc curve and calculating the AUC of roc curver.
roc1 <- roc(rc$overall.survival..days., response = rc$vital.status)
plot(roc1)
auc(roc1)

```

## 2.Naive Bayes Model

```
#data preparation
library(tigeR)
Dataloader(c(4,5,13,14,18))
SElist <- list(MEL_GSE91061, MEL_phs000452, RCC_Braun_2020)

#building model
mymodel <- build_Model('NB', SElist, Stem.Sig, response_NR = TRUE)

#read tigeR Built-in datasets
library(magrittr)
extract_mtr('MEL_GSE78220') %>% dataPreprocess(Stem.Sig, turn2HL = TRUE) -> test_Expr1
extract_mtr('MEL_PRJEB23709') %>% dataPreprocess(Stem.Sig, turn2HL = TRUE) -> test_Expr2
feature <- intersect(rownames(test_Expr1),rownames(test_Expr2))

test_Expr <- cbind(test_Expr1[rownames(test_Expr1) %in% feature,], test_Expr2[rownames(test_Expr2) %in% feature,])
response <- extract_label(c('MEL_GSE78220','MEL_PRJEB23709'))
value <- as.numeric(predict(mymodel, t(test_Expr), type = 'raw')[,1])

#Drawing roc curve and calculating the AUC of roc curver.
library(pROC)
ROC <- roc(response, value)
plot(ROC)
auc(ROC)
```
## 3.Random Forest Model

```
#data preparation
library(tigeR)
Dataloader(c(4,5,13,14,18))
SElist <- list(MEL_GSE91061, MEL_phs000452, RCC_Braun_2020)

#building model
mymodel <- build_Model('RF', SElist, Stem.Sig, rmBE = TRUE,response_NR = TRUE)

##testing model
library(pROC)

#read tigeR Built-in datasets
library(magrittr)
selected_gene <- rownames(mymodel$importance)
test_Expr1 <- extract_mtr('MEL_GSE78220')
test_Expr1 <- test_Expr1[rownames(test_Expr1) %in% selected_gene,]
test_Expr2 <- extract_mtr('MEL_PRJEB23709')
test_Expr2 <- test_Expr2[rownames(test_Expr2) %in% selected_gene,]

test_Expr <- cbind(test_Expr1, test_Expr2[rownames(test_Expr2) %in% rownames(test_Expr1),])
response <- extract_label(c('MEL_GSE78220','MEL_PRJEB23709'))
value <- as.numeric(predict(mymodel, t(test_Expr), type = 'vote')[,1])

#Drawing roc curve of patients whose prediction results are 'Responder' and calculating the AUC of roc curver.
ROC <- roc(response, value)
plot(ROC)
auc(ROC)

```
## 4.SVM Model

```
#data preparation
library(tigeR)
Dataloader(c(4,5,13,14,18))
SElist <- list(MEL_GSE91061, MEL_phs000452)

#building model
mymodel <- build_Model('SVM', SElist, Stem.Sig, rmBE = FALSE,response_NR = TRUE)

##testing model
library(pROC)

#read tigeR Built-in datasets
library(magrittr)
selected_gene <- colnames(mymodel$SV)
test_Expr1 <- extract_mtr('MEL_GSE78220')
test_Expr1 <- test_Expr1[rownames(test_Expr1) %in% selected_gene,]
test_Expr2 <- extract_mtr('MEL_PRJEB23709')
test_Expr2 <- test_Expr2[rownames(test_Expr2) %in% selected_gene,]

test_Expr <- cbind(test_Expr1, test_Expr2[rownames(test_Expr2) %in% rownames(test_Expr1),])
response <- extract_label(c('MEL_GSE78220','MEL_PRJEB23709'))
value <- predict(mymodel, t(test_Expr),type='eps-regression')

#Drawing roc curve of patients whose prediction results are 'Responder' and calculating the AUC of roc curver.
ROC <- roc(response, value)
plot(ROC)
auc(ROC)

```
## 5.Cancerclass Model

```
#data preparation
library(tigeR)
Dataloader(c(4,5,13,14,18))
SElist <- list(MEL_GSE91061, MEL_phs000452)

#building model
mymodel <- build_Model('CC', SElist, Stem.Sig, rmBE = TRUE)

exprs <- cbind(extract_mtr('MEL_GSE78220'), extract_mtr('MEL_PRJEB23709'))
exprs <- dataPreprocess(exprs, Stem.Sig, turn2HL = FALSE)
pData <- rbind(MEL_GSE78220@colData, MEL_PRJEB23709@colData)
rownames(pData) <- pData$sample_id
pData$class <- c(extract_label('MEL_GSE78220'), extract_label('MEL_PRJEB23709'))
identical(rownames(pData),colnames(exprs))
metadata <- data.frame(labelDescription = colnames(pData), row.names = colnames(pData))
adf <- new("AnnotatedDataFrame", data = as.data.frame(pData), varMetadata = metadata)

exampleSet2 <- new("ExpressionSet", exprs = exprs, phenoData = adf)
prediction <- cancerclass::predict(mymodel, exampleSet2, positive = "R", dist = "cor")
value <- as.numeric(prediction@prediction[,3])
response <- extract_label(c('MEL_GSE78220','MEL_PRJEB23709'))

library(pROC)
#Drawing roc curve of patients whose prediction results are 'Responder' and calculating the AUC of roc curver.
ROC <- roc(response, value)
plot(ROC)
auc(ROC)

```

## 6. Adaboost Model

```
#data preparation
library(tigeR)
Dataloader(c(4,5,13,14,18))
SElist <- list(MEL_GSE91061, MEL_phs000452)

#building model
mymodel <- build_Model('ADB', SElist, Stem.Sig, rmBE = FALSE)

##testing model
library(pROC)

#read tigeR Built-in datasets
library(magrittr)
selected_gene <- names(mymodel$importance)
test_Expr1 <- extract_mtr('MEL_GSE78220')
test_Expr1 <- test_Expr1[rownames(test_Expr1) %in% selected_gene,]
test_Expr2 <- extract_mtr('MEL_PRJEB23709')
test_Expr2 <- test_Expr2[rownames(test_Expr2) %in% selected_gene,]

library(adabag)
test_Expr <- cbind(test_Expr1, test_Expr2[rownames(test_Expr2) %in% rownames(test_Expr1),])
response <- extract_label(c('MEL_GSE78220','MEL_PRJEB23709'))
pred <- predict.boosting(mymodel, as.data.frame(t(test_Expr)))
value <- pred$votes[,1]

#Drawing roc curve of patients whose prediction results are 'Responder' and calculating the AUC of roc curver.
ROC <- roc(response, value)
plot(ROC)
auc(ROC)

```

## 7. Logitboost Model

```
#data preparation
library(tigeR)
Dataloader(c(4,5,13,14,18))
SElist <- list(MEL_GSE91061, MEL_phs000452)

#building model
mymodel <- build_Model('LGB', SElist, Stem.Sig, rmBE = TRUE)

##testing model
library(pROC)

#read tigeR Built-in datasets
library(magrittr)
extract_mtr('MEL_GSE78220') %>% dataPreprocess(Stem.Sig, turn2HL = FALSE) -> test_Expr1
extract_mtr('MEL_PRJEB23709') %>% dataPreprocess(Stem.Sig, turn2HL = FALSE) -> test_Expr2
feature <- intersect(rownames(test_Expr1),rownames(test_Expr2))

test_Expr <- cbind(test_Expr1[rownames(test_Expr1) %in% feature,], test_Expr2[rownames(test_Expr2) %in% feature,])
response <- extract_label(c('MEL_GSE78220','MEL_PRJEB23709'))
value <- predict(mymodel,t(test_Expr),type = 'raw')[,1]

#Drawing roc curve of patients whose prediction results are 'Responder' and calculating the AUC of roc curver.
ROC <- roc(response, value)
plot(ROC)
auc(ROC)

```
## 8. Logistics Model

```
#data preparation
library(tigeR)
Dataloader(c(4,5,13,14,18))
SElist <- list(MEL_GSE91061, MEL_phs000452)

#building model
mymodel <- build_Model('LGT', SElist, Stem.Sig[1:10], rmBE = FALSE, response_NR = TRUE)

##testing model
library(pROC)

#read tigeR Built-in datasets
library(magrittr)
extract_mtr('MEL_GSE78220') %>% dataPreprocess(Stem.Sig[1:10], turn2HL = FALSE) -> test_Expr1
extract_mtr('MEL_PRJEB23709') %>% dataPreprocess(Stem.Sig[1:10], turn2HL = FALSE) -> test_Expr2

test_Expr <- cbind(test_Expr1, test_Expr2[rownames(test_Expr2) %in% rownames(test_Expr1),])
response <- c(extract_label('MEL_GSE78220'), extract_label('MEL_PRJEB23709'))
df <- data.frame(t(max_min_normalization(test_Expr)))
value <- as.numeric(predict(mymodel, df, type = 'response'))

#Drawing roc curve of patients whose prediction results are 'Responder' and calculating the AUC of roc curver.
ROC <- roc(response, value)
plot(ROC)
auc(ROC)

```
## 9. Immunotherapy Response
```
library(tigeR)

DataLoader(5)
Immunotherapy_Response(gene='CD274', MEL_GSE91061)

```
## 10. Visualization
```
library(tigeR)

DataLoader(5)
plt_diff('CD274',MEL_GSE91061,'Treatment') # Treatment vs UnTreatment
plt_diff('CD274',MEL_GSE91061,'Response') # Responder vs Non-Responder
plt_surv('CD274',MEL_GSE91061) # Survival analysis

```
## 11. Cibersort
```
library(tigeR)

Dataloader(4)
mixture <- as.matrix(MEL_GSE78220_exp[,-1])
rownames(mixture) <- unlist(MEL_GSE78220_exp[,1])
data("LM22",package = "tigeR")

result <- CIBERSORT(LM22,mixture,perm=10, QN=T)

```
