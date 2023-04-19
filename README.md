# tigeR
## 1.Built-in Model

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
#Please load data set in Baidu cloud before running the code！
data(Stem.Sig,package = 'tigeR')

library(SummarizedExperiment)
SE1 <- MEL_GSE91061
SE2 <- MEL_phs000452
SE3 <- RCC_Braun_2020

SElist <- list(SE1, SE2, SE3)

#building model
library(tigeR)
mymodel <- build_NB_model(SElist, Stem.Sig, response_NR = TRUE)

##testing model
library(pROC)

#read tigeR Built-in datasets
library(magrittr)
extract_mtr('MEL_GSE78220_exp') %>% dataPreprocess(Stem.Sig, turn2HL = TRUE) -> test_Expr1
extract_mtr('MEL_PRJEB23709_exp') %>% dataPreprocess(Stem.Sig, turn2HL = TRUE) -> test_Expr2
test_Expr <- cbind(test_Expr1, test_Expr2[rownames(test_Expr2) %in% rownames(test_Expr1),])
test_response <- c(extract_label('MEL_GSE78220_meta'), extract_label('MEL_PRJEB23709_meta'))

#Obtaining the meta informations
rc <- rbind(MEL_GSE78220_meta, MEL_PRJEB23709_meta)
rc$response %<>% sub('CR|MR|PR|CRPR', 'R',.) %>% sub('PD|SD', 'NR',.) %>% as.factor()

value <- as.numeric(predict(mymodel, t(test_Expr), type = 'raw')[,1])
#Drawing roc curve and calculating the AUC of roc curver.
roc1 <- roc(rc$response, value)
plot(roc1)
auc(roc1)

```
## 3.Random Forest Model

```
#Please load data set in Baidu cloud before running the code！
data(Stem.Sig, package = "tigeR")

library(SummarizedExperiment)
SE1 <- MEL_GSE91061
SE2 <- MEL_phs000452
SE3 <- RCC_Braun_2020

SElist <- list(SE1, SE2)

#building model
library(tigeR)
mymodel <- build_RF_model(SElist, Stem.Sig, rmBE = TRUE)

##testing model
library(pROC)

#read tigeR Built-in datasets
library(magrittr)
extract_mtr('MEL_GSE78220_exp') %>% dataPreprocess(Stem.Sig, turn2HL = FALSE) -> test_Expr1
extract_mtr('MEL_PRJEB23709_exp') %>% dataPreprocess(Stem.Sig, turn2HL = FALSE) -> test_Expr2
test_Expr <- cbind(test_Expr1, test_Expr2[rownames(test_Expr2) %in% rownames(test_Expr1),])
test_response <- c(extract_label('MEL_GSE78220_meta'), extract_label('MEL_PRJEB23709_meta'))

#Obtaining the meta informations
rc <- rbind(MEL_GSE78220_meta, MEL_PRJEB23709_meta)
rc$response %<>% sub('CR|MR|PR|CRPR', 'R',.) %>% sub('PD|SD', 'NR',.) %>% as.factor()

value <- as.numeric(predict(mymodel, t(test_Expr), type = 'vote')[,1])
#Drawing roc curve of patients whose prediction results are 'Responder' and calculating the AUC of roc curver.
roc1 <- roc(rc$response, value)
plot(roc1)
auc(roc1)


```

## 5.Cancerclass Model

```
#Please load data set in Baidu cloud before running the code！
data(Stem.Sig, package = "tigeR")

#standardization of expression matrix
exp1 <- as.matrix(MEL_GSE91061_exp[,-1])
exp1 <- apply(exp1, 2, as.numeric)
rownames(exp1) <- MEL_GSE91061_exp[,1]

exp2 <- as.matrix(MEL_phs000452_exp[,-1])
exp2 <- apply(exp2, 2, as.numeric)
rownames(exp2) <- MEL_phs000452_exp[,1]

meta1 <- MEL_GSE91061_meta
meta2 <- MEL_phs000452_meta

library(tigeR)
response1 <- meta1$response
response1 <- response_standardize(response1)
response2 <- meta2$response
response2 <- response_standardize(response2)

#obtain index of NA samples
filt1 <- grep('UNK',response1)
filt2 <- grep('UNK',response2)

#remove UNK. exp2 has no UNK
exp1 <- exp1[,-filt1]
response1 <- response1[-filt1]
#exp2 <- exp2[,-filt1]
#response2 <- response2[-filt2]

library(SummarizedExperiment)
colData1 <- DataFrame(response = response1)
SE1 <- SummarizedExperiment(assays = exp1,
                            colData = colData1)

colData2 <- DataFrame(response = response2)
SE2 <- SummarizedExperiment(assays = exp2,
                            colData = colData2)

SElist <- list(SE1, SE2)


data("MEL_GSE78220_meta")
data("MEL_PRJEB23709_meta")
#building model
mymodel <- build_CC_model(SElist, Stem.Sig, rmBE = TRUE)

exprs <- cbind(extract_mtr('MEL_GSE78220_exp'), extract_mtr('MEL_PRJEB23709_exp'))
exprs <- dataPreprocess(exprs, Stem.Sig, turn2HL = FALSE)
pData <- rbind(MEL_GSE78220_meta, MEL_PRJEB23709_meta)
rownames(pData) <- pData$sample_id
pData$class <- c(extract_label('MEL_GSE78220_meta'), extract_label('MEL_PRJEB23709_meta'))
identical(rownames(pData),colnames(exprs))
metadata <- data.frame(labelDescription = colnames(pData), row.names = colnames(pData))
adf <- new("AnnotatedDataFrame", data = as.data.frame(pData), varMetadata = metadata)

exampleSet2 <- new("ExpressionSet", exprs = exprs, phenoData = adf)
prediction <- cancerclass::predict(mymodel, exampleSet2, positive = "R", dist = "cor")
value <- as.numeric(prediction@prediction[,3])

##testing model
library(pROC)

#read tigeR Built-in datasets
test_Expr <- cbind(extract_mtr('MEL_GSE78220_exp'), extract_mtr('MEL_PRJEB23709_exp'))
test_Expr <- dataPreprocess(test_Expr, Stem.Sig, turn2HL = FALSE)
test_response <- c(extract_label('MEL_GSE78220_meta'), extract_label('MEL_PRJEB23709_meta'))

#Obtaining the meta informations of patients whose prediction results are 'Response'.
data("MEL_GSE78220_meta")
data("MEL_PRJEB23709_meta")

rc <- rbind(MEL_GSE78220_meta, MEL_PRJEB23709_meta)
rc$response <- sub('CR|MR|PR|CRPR', 'R', rc$response)
rc$response <- sub('PD|SD', 'NR', rc$response)
rc$response <- as.factor(rc$response)

#Drawing roc curve of patients whose prediction results are 'Responder' and calculating the AUC of roc curver.
roc <- roc(rc$response, value)
plot(roc)
auc(roc)

```

## 6. Adaboost Model

```
#Please load data set in Baidu cloud before running the code！
data(Stem.Sig, package = "tigeR")

#standardization of expression matrix
exp1 <- as.matrix(MEL_GSE91061_exp[,-1])
exp1 <- apply(exp1, 2, as.numeric)
rownames(exp1) <- MEL_GSE91061_exp[,1]

exp2 <- as.matrix(MEL_phs000452_exp[,-1])
exp2 <- apply(exp2, 2, as.numeric)
rownames(exp2) <- MEL_phs000452_exp[,1]

meta1 <- MEL_GSE91061_meta
meta2 <- MEL_phs000452_meta

library(tigeR)
response1 <- meta1$response
response1 <- response_standardize(response1)
response2 <- meta2$response
response2 <- response_standardize(response2)

#obtain index of NA samples
filt1 <- grep('UNK',response1)
filt2 <- grep('UNK',response2)

#remove UNK. exp2 has no UNK
exp1 <- exp1[,-filt1]
response1 <- response1[-filt1]
#exp2 <- exp2[,-filt1]
#response2 <- response2[-filt2]

library(SummarizedExperiment)
colData1 <- DataFrame(response = response1)
SE1 <- SummarizedExperiment(assays = exp1,
                            colData = colData1)

colData2 <- DataFrame(response = response2)
SE2 <- SummarizedExperiment(assays = exp2,
                            colData = colData2)

SElist <- list(SE1, SE2)


data("MEL_GSE78220_meta")
data("MEL_PRJEB23709_meta")
#building model
library(adabag)
mymodel <- build_Adaboost_model(SElist, Stem.Sig, rmBE = FALSE)

exp_test <- cbind(extract_mtr('MEL_GSE78220_exp'), extract_mtr('MEL_PRJEB23709_exp'))
exp_test <- dataPreprocess(exp_test, Stem.Sig, turn2HL = FALSE)

pred <- predict.boosting(mymodel, as.data.frame(t(exp_test)))

value <- pred$votes[,1]

##testing model
library(pROC)

#read tigeR Built-in datasets
test_Expr1 <- extract_mtr('MEL_GSE78220_exp')
test_Expr1 <- dataPreprocess(test_Expr1, Stem.Sig, turn2HL = FALSE)
test_Expr2 <- extract_mtr('MEL_PRJEB23709_exp')
test_Expr2 <- dataPreprocess(test_Expr2, Stem.Sig, turn2HL = FALSE)
test_Expr <- cbind(test_Expr1, test_Expr2[rownames(test_Expr2) %in% rownames(test_Expr1),])
test_response <- c(extract_label('MEL_GSE78220_meta'), extract_label('MEL_PRJEB23709_meta'))

#Obtaining the meta informations of patients whose prediction results are 'Response'.

rc <- rbind(MEL_GSE78220_meta, MEL_PRJEB23709_meta)
rc$response <- sub('CR|MR|PR|CRPR', 'R', rc$response)
rc$response <- sub('PD|SD', 'NR', rc$response)
rc$response <- as.factor(rc$response)

#Drawing roc curve of patients whose prediction results are 'Responder' and calculating the AUC of roc curver.
roc1 <- roc(rc$response, value)
plot(roc1)
auc(roc1)

```

## 7. Logitboost Model

```
#Please load data set in Baidu cloud before running the code！
data(Stem.Sig, package = "tigeR")

#standardization of expression matrix
exp1 <- as.matrix(MEL_GSE91061_exp[,-1])
exp1 <- apply(exp1, 2, as.numeric)
rownames(exp1) <- MEL_GSE91061_exp[,1]

exp2 <- as.matrix(MEL_phs000452_exp[,-1])
exp2 <- apply(exp2, 2, as.numeric)
rownames(exp2) <- MEL_phs000452_exp[,1]

meta1 <- MEL_GSE91061_meta
meta2 <- MEL_phs000452_meta

library(tigeR)
response1 <- meta1$response
response1 <- response_standardize(response1)
response2 <- meta2$response
response2 <- response_standardize(response2)

#obtain index of NA samples
filt1 <- grep('UNK',response1)
filt2 <- grep('UNK',response2)

#remove UNK. exp2 has no UNK
exp1 <- exp1[,-filt1]
response1 <- response1[-filt1]
#exp2 <- exp2[,-filt1]
#response2 <- response2[-filt2]

library(SummarizedExperiment)
colData1 <- DataFrame(response = response1)
SE1 <- SummarizedExperiment(assays = exp1,
                            colData = colData1)

colData2 <- DataFrame(response = response2)
SE2 <- SummarizedExperiment(assays = exp2,
                            colData = colData2)

SElist <- list(SE1, SE2)

#building model
library(tigeR)
mymodel <- build_Logitboost_model(SElist, Stem.Sig, rmBE = TRUE)

##testing model
library(pROC)

#read tigeR Built-in datasets
test_Expr1 <- extract_mtr('MEL_GSE78220_exp')
test_Expr1 <- dataPreprocess(test_Expr1, Stem.Sig, turn2HL = FALSE)
test_Expr2 <- extract_mtr('MEL_PRJEB23709_exp')
test_Expr2 <- dataPreprocess(test_Expr2, Stem.Sig, turn2HL = FALSE)
test_Expr <- cbind(test_Expr1, test_Expr2[rownames(test_Expr2) %in% rownames(test_Expr1),])
test_response <- c(extract_label('MEL_GSE78220_meta'), extract_label('MEL_PRJEB23709_meta'))

#Obtaining the meta informations of patients whose prediction results are 'Response'.

rc <- rbind(MEL_GSE78220_meta, MEL_PRJEB23709_meta)
rc$response <- sub('CR|MR|PR|CRPR', 'R', rc$response)
rc$response <- sub('PD|SD', 'NR', rc$response)
rc$response <- as.factor(rc$response)

#Drawing roc curve of patients whose prediction results are 'Responder' and calculating the AUC of roc curver.
roc1 <- roc(rc$response, value)
plot(roc1)
auc(roc1)

```
