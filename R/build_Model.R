#' @title Build machine learning prediction model for immunotherapy response
#' @description Generate immunotherapy prognosis prediction model.
#' @param SE the dataset you wish to use to build your model. A SummarizedExperiment (SE) object,
#' which can be either a single SE object or a list of SE objects. Note that for each SE object,
#' the colData must contain treatment information under the column name Treatment.
#' @param mtr the dataset you wish to use to build your model. A SummarizedExperiment (SE) object, which can be either a single SE object or a list of SE objects. Note that for each SE object, the colData must contain treatment information under the column name Treatment.
#' @param meta refers to the specific set of genes you wish to use for model construction.
#' @param Model represents the type of model you want to build. You have several options to choose from: "NB" for Naive Bayes, "SVM" for Support Vector Machine, "RF" for Random Forest, "CC" for Cancerclass, "ADB" for Adaboost, "LGB" for Logitboost, and "LGT" for Logistics.
#' @param feature_genes refers to the specific set of genes you wish to use for model construction.
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR If TRUE, only use R or NR to represent Immunotherapy response of patients.
#' @param ... the arguments
#' @export

build_Model <- function(SE, mtr, meta, Model, feature_genes, rmBE = FALSE, response_NR = TRUE, ...){
  UseMethod("build_Model")
}



#' @title Build machine learning prediction model for immunotherapy response
#' @description Generate immunotherapy prognosis prediction model.
#' @param mtr the dataset you wish to use to build your model. A SummarizedExperiment (SE) object, which can be either a single SE object or a list of SE objects. Note that for each SE object, the colData must contain treatment information under the column name Treatment.
#' @param meta refers to the specific set of genes you wish to use for model construction.
#' @param Model represents the type of model you want to build. You have several options to choose from: "NB" for Naive Bayes, "SVM" for Support Vector Machine, "RF" for Random Forest, "CC" for Cancerclass, "ADB" for Adaboost, "LGB" for Logitboost, and "LGT" for Logistics.
#' @param response_NR If TRUE, only use R or NR to represent Immunotherapy response of patients.
#' @param ... the arguments
#' @export

build_Model.matrix <- function(mtr, meta, Model, response_NR = TRUE, ...){
  SE_obj <- SummarizedExperiment::SummarizedExperiment(assays=S4Vectors::SimpleList(mtr),
                                                       colData=S4Vectors::DataFrame(meta),
                                                       checkDimnames=TRUE)
  build_Model.default(SE_obj, Model, feature_genes = NULL, rmBE = FALSE, response_NR, ...)
}


#' @title Build machine learning prediction model for immunotherapy response
#' @description Generate immunotherapy prognosis prediction model.
#' @param SE the dataset you wish to use to build your model. A SummarizedExperiment (SE) object, which can be either a single SE object or a list of SE objects. Note that for each SE object, the colData must contain treatment information under the column name Treatment.
#' @param Model represents the type of model you want to build. You have several options to choose from: "NB" for Naive Bayes, "SVM" for Support Vector Machine, "RF" for Random Forest, "CC" for Cancerclass, "ADB" for Adaboost, "LGB" for Logitboost, and "LGT" for Logistics.
#' @param feature_genes refers to the specific set of genes you wish to use for model construction.
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR If TRUE, only use R or NR to represent Immunotherapy response of patients.
#' @param ... the arguments
#' @export

build_Model.default <- function(SE, Model, feature_genes, rmBE = FALSE, response_NR = TRUE, ...){
  if(Model == 'NB')
    model <- build_NB_model(SE, feature_genes, rmBE, response_NR, ...)
  else if(Model == 'RF')
    model <- build_RF_model(SE, feature_genes, rmBE, response_NR, ...)
  else if(Model == 'SVM')
    model <- build_SVM_model(SE, feature_genes, rmBE, response_NR, ...)
  else if(Model == 'CC')
    model <- build_CC_model(SE, feature_genes, rmBE, response_NR, ...)
  else if(Model == 'ADB')
    model <- build_Adaboost_model(SE, feature_genes, rmBE, response_NR, ...)
  else if(Model == 'LGB')
    model <- build_Logitboost_model(SE, feature_genes, rmBE, response_NR, ...)
  else if(Model == 'LGT')
    model <- build_Logistics_model(SE, feature_genes, rmBE, response_NR, ...)
  else
    stop("Please check your parameter! Avaliable value of Model('NB','SVM','RF','CC','ADB','LGB','LGT').")
  return(model)
}


#' @title Build naive bayes prediction model for immunotherapy response
#' @description Generate a naive bayes model.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR If TRUE, only use R or NR to represent Immunotherapy response of patients.
#' @param laplace positive double controlling Laplace smoothing. The default (0) disables Laplace smoothing.
#' @param ... the arguments
#' @importFrom e1071 naiveBayes

build_NB_model <- function(SE, Signature, rmBE = FALSE, response_NR = TRUE, laplace=1, ...){
  data <- dataProcess(SE, Signature, rmBE, response_NR, TRUE)
  model <- naiveBayes(t(data[[1]]), data[[2]]$response_NR, laplace, ...)
  model$threshold <- data[[3]]
  model
}


#' @title Build SVM prediction model for immunotherapy response
#' @description Generate a pport Vector Machine model.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR If TRUE, only use R or NR to represent Immunotherapy response of patients.
#' @param type the kernel used in training and predicting.
#' @param probability logical indicating whether the model should allow for probability predictions. description
#' @param ... the arguments
#' @importFrom stats na.omit

build_SVM_model <- function(SE, Signature, rmBE = TRUE, response_NR, type = 'eps', probability = TRUE, ...){
  data <- dataProcess(SE, Signature, rmBE, response_NR, FALSE)
  model <- e1071::svm(x = t(na.omit(data[[1]])),
                      y = as.numeric(as.factor(data[[2]]$response)),
                      probability, type, ...)
  model$features <- Signature
  return(model)
}

#' @title Build random forest prediction model for immunotherapy response
#' @description Generate a random forest model.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR If TRUE, only use R or NR to represent Immunotherapy response of patients.
#' @param ... the arguments
#' @importFrom randomForest randomForest
#' @importFrom stats na.omit

build_RF_model <- function(SE, Signature, rmBE = FALSE, response_NR = TRUE, ...){
  data <- dataProcess(SE, Signature, rmBE, response_NR, FALSE)
  v_Args <- list(...)
  Args <- c(list(x = t(na.omit(data[[1]])),
                 y = as.factor(data[[2]]$response),
                 ntree = ifelse("ntree" %in% names(v_Args),v_Args[["ntree"]],150)),
            v_Args[names(v_Args) != "ntree"])
  model <- do.call(randomForest::randomForest, Args)
}

#' @title Build cancerclass prediction model for immunotherapy response
#' @description Generate a cancerclass model.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR If TRUE, only use R or NR to represent Immunotherapy response of patients.
#' @param ... the arguments
#' @importFrom cancerclass fit

build_CC_model <- function(SE, Signature, rmBE = TRUE, response_NR = TRUE, ...){
  data <- dataProcess(SE, Signature, rmBE, response_NR, FALSE)
  pData <- data.frame(class = data[[2]]$response, row.names = colnames(data[[1]]))
  metadata <- data.frame(labelDescription = colnames(pData), row.names = colnames(pData))
  adf <- new("AnnotatedDataFrame", data = pData, varMetadata = metadata)
  exampleSet <- methods::new("ExpressionSet", exprs = data[[1]], phenoData = adf)
  model <- do.call(cancerclass::fit, c(list(exampleSet),list(...)))
  colnames(model@predictor)  <- c("welch.test","R","N")
  model
}


#' @title Build Adaboost prediction model for immunotherapy response
#' @description Generate a Adaboost model.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR If TRUE, only use R or NR to represent Immunotherapy response of patients.
#' @param ... the arguments

build_Adaboost_model <- function(SE, Signature, rmBE = TRUE, response_NR = TRUE, ...){
  data <- dataProcess(SE, Signature, rmBE, response_NR, FALSE)
  df <- data.frame(class = data[[2]]$response, t(data[[1]]))
  df$class <- factor(df$class)
  v_Args <- list(...)
  Args <- c(list(formula=class ~ .,
                 data=df,
                 mfinal=ifelse("mfinal" %in% names(v_Args),v_Args[["mfinal"]],10)),
            v_Args[names(v_Args) != "mfinal"])
  model <- do.call(adabag::boosting, Args)
  model$features <- rownames(data[[1]])
  model
}



#' @title Build Logitboost prediction model for immunotherapy response
#' @description Generate a LogitBoost model.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR If TRUE, only use R or NR to represent Immunotherapy response of patients.
#' @param ... the arguments

build_Logitboost_model <- function(SE, Signature, rmBE = TRUE, response_NR = TRUE, ...){
  v_Args <- list(...)
  nIter <- ifelse("nIter" %in% names(v_Args),v_Args[["nIter"]],150)
  data <- dataProcess(SE, Signature, rmBE, response_NR, FALSE)
  model <- LogitBoost(xlearn = t(data[[1]]), ylearn = factor(data[[2]]$response), nIter = nIter)
}


#' @title Build Logistics prediction model for immunotherapy response
#' @description Generate a Logistics model.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR If TRUE, only use R or NR to represent Immunotherapy response of patients.
#' @param ... the arguments

build_Logistics_model <- function(SE, Signature, rmBE = FALSE, response_NR = TRUE, ...){
  data <- dataProcess(SE, Signature, rmBE, response_NR, FALSE)
  df <- data.frame(response=ifelse(data[[2]]$response=='R',1,0),t(data[[1]]))

  v_Args <- list(...)
  Args <- c(list(formula=response ~.,
                 data=df,
                 family=stats::binomial(link = "logit"),
                 control=ifelse("control" %in% names(v_Args),v_Args[["control"]],list(maxit=100))),
            v_Args[names(v_Args) != "control"])
  model <- do.call(stats::glm, Args)
}

