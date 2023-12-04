#' @title Prepare expression matrix for down string analysis
#' @description dataPreprocess() will remove missing genes. Then returns the sub-matrix of the genes whose SYMBOLs are in the signature.
#' @param exp_mtr An expression matrix which rownames are gene SYMBOL and colnames are sample ID.
#' @param Signature The aiming gene set(only Gene SYMBOL allowed).
#' @param turn2HL If TRUE, the expression value of a gene is divided to "HIGH" or "LOW" based on its median expression.
#' @export
#'

dataPreprocess <- function(exp_mtr, Signature = NULL, turn2HL = TRUE){
  if(is.null(Signature))
    Signature <- rownames(exp_mtr)

  exp_mtr[is.na(exp_mtr)] <- 0
  rowname <- rownames(exp_mtr)
  colname <- colnames(exp_mtr)
  exp_mtr <- apply(exp_mtr, 2, as.numeric)

  rownames(exp_mtr) <- rowname
  colnames(exp_mtr) <- colname

  exp_mtr <- exp_mtr[rownames(exp_mtr) %in% unlist(Signature),]
  colname <- colnames(exp_mtr)
  exp_mtr <- t(apply(exp_mtr, 1, zero2na))   #converse lines which full of zero to NA

  rowname <- rownames(exp_mtr)
  exp_mtr <- apply(exp_mtr, 2, as.numeric)

  colnames(exp_mtr) <- colname
  rownames(exp_mtr) <- rowname

  if(all(is.na(exp_mtr)))
    return(exp_mtr)

  filt_NA_mtr <- exp_mtr[!apply(exp_mtr, 1, is.NA_vec),]

  if (nrow(filt_NA_mtr) != nrow(exp_mtr)){
    NA_percentage <- round((1 - nrow(filt_NA_mtr) / nrow(exp_mtr)) * 100, digits = 2)
    message(paste0(NA_percentage,'% genes in Signature are abscent in expression matrix!It may affect model performance!'))
  }

  if(turn2HL){
    count_mean_mtr <- filt_NA_mtr
    count_mean_mtr[is.na(count_mean_mtr)] <- 0
    mean_vec <- apply(count_mean_mtr, 1, mean)

    for (i in 1:length(filt_NA_mtr[,1])) {
      for (j in 1:length(filt_NA_mtr[1,])) {
        if(is.na(filt_NA_mtr[i,j])){
          next
        }else
          if(filt_NA_mtr[i,j] >= mean_vec[i])
            filt_NA_mtr[i,j] <- 'HIGH'
          else
            filt_NA_mtr[i,j] <- 'LOW'
      }
    }
    exp_mtr[!apply(exp_mtr, 1, is.NA_vec),] <- filt_NA_mtr
  }
  else {
    exp_mtr <- t(apply(exp_mtr, 1, zero2na))
  }

  return(exp_mtr[!apply(exp_mtr, 1, is.NA_vec),])
}

#' @title Build naive bayes prediction model for immunotherapy response
#' @description Generate a naive bayes model.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR If TRUE, only use R or NR to represent Immunotherapy response of patients.
#' @importFrom e1071 naiveBayes

build_NB_model <- function(SE, Signature, rmBE = FALSE, response_NR = TRUE){
  data <- dataProcess(SE, Signature, rmBE, response_NR, TRUE)
  model <- naiveBayes(t(data[[1]]), data[[2]]$response, laplace = 1)
}


#' @title Build SVM prediction model for immunotherapy response
#' @description Generate a pport Vector Machine model.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR If TRUE, only use R or NR to represent Immunotherapy response of patients.
#' @importFrom stats na.omit

build_SVM_model <- function(SE, Signature, rmBE = TRUE, response_NR){
  data <- dataProcess(SE, Signature, rmBE, response_NR, FALSE)
  model <- e1071::svm(x = t(na.omit(data[[1]])),
                      y = as.numeric(as.factor(data[[2]]$response)),
                      scale = TRUE,
                      type = 'eps',
                      kernel = 'radial',
                      probability = TRUE)
  return(model)
}

#' @title Build random forest prediction model for immunotherapy response
#' @description Generate a random forest model.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR If TRUE, only use R or NR to represent Immunotherapy response of patients.
#' @importFrom randomForest randomForest
#' @importFrom stats na.omit

build_RF_model <- function(SE, Signature, rmBE = FALSE, response_NR = TRUE){
  data <- dataProcess(SE, Signature, rmBE, response_NR, FALSE)
  model <- randomForest::randomForest(x = t(na.omit(data[[1]])),
                                      y = as.factor(data[[2]]$response),
                                      ntree = 150,
                                      mtry = 9)
}

#' @title Build cancerclass prediction model for immunotherapy response
#' @description Generate a cancerclass model.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR If TRUE, only use R or NR to represent Immunotherapy response of patients.
#' @importFrom cancerclass fit

build_CC_model <- function(SE, Signature, rmBE = TRUE, response_NR = TRUE){
  data <- dataProcess(SE, Signature, rmBE, response_NR, FALSE)
  pData <- data.frame(class = data[[2]]$response, row.names = colnames(data[[1]]))
  metadata <- data.frame(labelDescription = colnames(pData), row.names = colnames(pData))
  adf <- new("AnnotatedDataFrame", data = pData, varMetadata = metadata)
  exampleSet <- methods::new("ExpressionSet", exprs = data[[1]], phenoData = adf)
  model <- cancerclass::fit(exampleSet, method = "welch.test")
}


#' @title Build Adaboost prediction model for immunotherapy response
#' @description Generate a Adaboost model.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR If TRUE, only use R or NR to represent Immunotherapy response of patients.
#' @importFrom adabag boosting

build_Adaboost_model <- function(SE, Signature, rmBE = TRUE, response_NR = TRUE){
  data <- dataProcess(SE, Signature, rmBE, response_NR, FALSE)
  df <- data.frame(class = data[[2]]$response, t(data[[1]]))
  df$class <- factor(df$class)
  model <- adabag::boosting(class ~ ., data = df, mfinal = 10, boos = TRUE)
}



#' @title Build Logitboost prediction model for immunotherapy response
#' @description Generate a LogitBoost model.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR If TRUE, only use R or NR to represent Immunotherapy response of patients.
#' @importFrom stats glm
#' @importFrom stats binomial

build_Logitboost_model <- function(SE, Signature, rmBE = TRUE, response_NR = TRUE){
  data <- dataProcess(SE, Signature, rmBE, response_NR, FALSE)
  model <- LogitBoost(xlearn = t(data[[1]]), ylearn = factor(data[[2]]$response), nIter = 300)
}


#' @title Build Logistics prediction model for immunotherapy response
#' @description Generate a Logistics model.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR If TRUE, only use R or NR to represent Immunotherapy response of patients.
#' @importFrom caTools LogitBoost

build_Logistics_model <- function(SE, Signature, rmBE = FALSE, response_NR = TRUE){
  data <- dataProcess(SE, Signature, rmBE, response_NR, FALSE)
  df <- data.frame(response=ifelse(data[[2]]$response=='R',1,0),t(data[[1]]))
  model <- glm(response ~.,data=df,family = binomial(link = "logit"),control=list(maxit=100))
}


#' @title Build machine learning prediction model for immunotherapy response
#' @description Generate immunotherapy prognosis prediction model.
#' @param Model represents the type of model you want to build. You have several options to choose from: "NB" for Naive Bayes, "SVM" for Support Vector Machine, "RF" for Random Forest, "CC" for Cancerclass, "ADB" for Adaboost, "LGB" for Logitboost, and "LGT" for Logistics.
#' @param SE the dataset you wish to use to build your model. A SummarizedExperiment (SE) object, which can be either a single SE object or a list of SE objects. Note that for each SE object, the colData must contain treatment information under the column name Treatment.
#' @param feature_genes refers to the specific set of genes you wish to use for model construction.
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR If TRUE, only use R or NR to represent Immunotherapy response of patients.
#' @export

build_Model <- function(Model, SE, feature_genes, rmBE = FALSE, response_NR = TRUE){
  if(Model == 'NB')
    model <- build_NB_model(SE, feature_genes, rmBE, response_NR)
  else if(Model == 'RF')
    model <- build_RF_model(SE, feature_genes, rmBE, response_NR)
  else if(Model == 'SVM')
    model <- build_SVM_model(SE, feature_genes, rmBE, response_NR)
  else if(Model == 'CC')
    model <- build_CC_model(SE, feature_genes, rmBE, response_NR)
  else if(Model == 'ADB')
    model <- build_Adaboost_model(SE, feature_genes, rmBE, response_NR)
  else if(Model == 'LGB')
    model <- build_Logitboost_model(SE, feature_genes, rmBE, response_NR)
  else if(Model == 'LGT')
    model <- build_Logistics_model(SE, feature_genes, rmBE, response_NR)
  else
    stop("Please check your parameter! Avaliable value of Model('NB','SVM','RF','CC','ADB','LGB','LGT').")
  return(model)
}
