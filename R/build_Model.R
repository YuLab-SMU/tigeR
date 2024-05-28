#' @title Build machine learning prediction model for immunotherapy response
#' @description generate immunotherapy response prediction machine-learning model.
#' @param ... the arguments
#' @section S3 methods:
#' \describe{
#'   \item{build_Model.matrix}{\code{\link{build_Model.matrix}}: Plots objects of class "myclass". See \code{?build_Model.matrix} for details.}
#' }
#' @rdname build_Model-build_Model.matrix-build_Model.default
#' @export

build_Model <- function(...){
  UseMethod("build_Model")
}



#' @title Build machine learning prediction model for immunotherapy response
#' @param mtr the dataset you wish to use to build your model. A SummarizedExperiment (SE) object, which can be either a single SE object or a list of SE objects. Note that for each SE object, the colData must contain treatment information under the column name Treatment.
#' @param meta refers to the specific set of genes you wish to use for model construction.
#' @param Model represents the type of model you want to build. You have several options to choose from: "NB" for Naive Bayes, "SVM" for Support Vector Machine, "RF" for Random Forest, "CC" for Cancerclass, "ADB" for Adaboost, "LGB" for Logitboost, and "LGT" for Logistics.
#' @param response_NR If TRUE, classify patients with CR, MR, PR as Responders (R), and those with PD, SD, NR as Non-Responders(NR).
#' @param ... the arguments
#' @rdname build_Model-build_Model.matrix-build_Model.default
#' @export

build_Model.matrix <- function(mtr, meta, Model, response_NR = TRUE, ...){
  SE_obj <- SummarizedExperiment::SummarizedExperiment(assays=S4Vectors::SimpleList(mtr),
                                                       colData=S4Vectors::DataFrame(meta),
                                                       checkDimnames=TRUE)
  build_Model.default(SE_obj, Model, feature_genes = NULL, rmBE = FALSE, response_NR, ...)
}


#' @title Build machine learning prediction model for immunotherapy response
#' @param SE the dataset you wish to use to build your model. A SummarizedExperiment (SE) object, which can be either a single SE object or a list of SE objects. Note that for each SE object, the colData must contain treatment information under the column name Treatment.
#' @param Model represents the type of model you want to build. You have several options to choose from: "NB" for Naive Bayes, "SVM" for Support Vector Machine, "RF" for Random Forest, "CC" for Cancerclass, "ADB" for Adaboost, "LGB" for Logitboost, and "LGT" for Logistics.
#' @param feature_genes refers to the specific set of genes you wish to use for model construction.
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR If TRUE, classify patients with CR, MR, PR as Responders (R), and those with PD, SD, NR as Non-Responders(NR).
#' @param PT_drop If TRUE, only Untreated patient will be use for model training.
#' @param ... the arguments
#' @rdname build_Model-build_Model.matrix-build_Model.default
#' @export

build_Model.default <- function(SE, Model, feature_genes, rmBE = FALSE, response_NR = TRUE, PT_drop=TRUE, ...){
  model <- switch (Model,
                   NB = build_NB_model(SE, feature_genes, rmBE, response_NR, PT_drop=PT_drop, ...),
                   RF = build_RF_model(SE, feature_genes, rmBE, response_NR, PT_drop=PT_drop, ...),
                   SVM = build_SVM_model(SE, feature_genes, rmBE, response_NR,PT_drop=PT_drop, ...),
                   CC = build_CC_model(SE, feature_genes, rmBE, response_NR, PT_drop=PT_drop, ...),
                   ADB = build_Adaboost_model(SE, feature_genes, rmBE, response_NR, PT_drop=PT_drop, ...),
                   LGB = build_Logitboost_model(SE, feature_genes, rmBE, response_NR, PT_drop=PT_drop, ...),
                   LGT = build_Logistics_model(SE, feature_genes, rmBE, response_NR, PT_drop=PT_drop, ...),
                   SURV = build_SURV_Model(SE, feature_genes, rmBE, PT_drop=PT_drop, ...))

  if(is.null(model))
    stop("Please check your parameter! Avaliable value of Model('NB','SVM','RF','CC','ADB','LGB','LGT','SURV').")

  return(model)
}


#' @title Build naive bayes prediction model for immunotherapy response
#' @description Generate a naive bayes model.
#' @param SE a SummarizedExperiment(SE) object or a list consists of multiple SE objects. The colData of the SE object(s) must contain treatment information named Treatment.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data sets using internal Combat method.
#' @param response_NR If TRUE, classify patients with CR, MR, PR as Responders (R), and those with PD, SD, NR as Non-Responders(NR).
#' @param laplace positive double controlling Laplace smoothing. The default (0) disables Laplace smoothing.
#' @param PT_drop If TRUE, only Untreated patient will be use for model training.
#' @param ... the arguments

build_NB_model <- function(SE, Signature, rmBE = FALSE, response_NR = TRUE, laplace=1, PT_drop,...){
  data <- dataProcess(SE, Signature, rmBE, response_NR, TRUE)
  if(PT_drop)
    data <- PT_filter(data)

  model <- e1071::naiveBayes(t(data[[1]]), data[[2]]$response_NR, laplace, ...)
  model$threshold <- data[[3]]
  model
}


#' @title Build SVM prediction model for immunotherapy response
#' @description Generate a pport Vector Machine model.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR If TRUE, classify patients with CR, MR, PR as Responders (R), and those with PD, SD, NR as Non-Responders(NR).
#' @param type the kernel used in training and predicting.
#' @param probability logical indicating whether the model should allow for probability predictions. description
#' @param PT_drop If TRUE, only Untreated patient will be use for model training.
#' @param ... the arguments
#' @importFrom stats na.omit

build_SVM_model <- function(SE, Signature, rmBE = TRUE, response_NR, type = 'eps', probability = TRUE, PT_drop, ...){
  data <- dataProcess(SE, Signature, rmBE, response_NR, FALSE)
  if(PT_drop)
    data <- PT_filter(data)
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
#' @param response_NR If TRUE, classify patients with CR, MR, PR as Responders (R), and those with PD, SD, NR as Non-Responders(NR).
#' @param PT_drop If TRUE, only Untreated patient will be use for model training.
#' @param ... the arguments
#' @importFrom randomForest randomForest
#' @importFrom stats na.omit

build_RF_model <- function(SE, Signature, rmBE = FALSE, response_NR = TRUE, PT_drop, ...){
  data <- dataProcess(SE, Signature, rmBE, response_NR, FALSE)
  if(PT_drop)
    data <- PT_filter(data)
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
#' @param response_NR If TRUE, classify patients with CR, MR, PR as Responders (R), and those with PD, SD, NR as Non-Responders(NR).
#' @param PT_drop If TRUE, only Untreated patient will be use for model training.
#' @param ... the arguments
#' @importFrom cancerclass fit

build_CC_model <- function(SE, Signature, rmBE = TRUE, response_NR = TRUE, PT_drop,...){
  data <- dataProcess(SE, Signature, rmBE, response_NR, FALSE)
  if(PT_drop)
    data <- PT_filter(data)
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
#' @param response_NR If TRUE, classify patients with CR, MR, PR as Responders (R), and those with PD, SD, NR as Non-Responders(NR).
#' @param PT_drop If TRUE, only Untreated patient will be use for model training.
#' @param ... the arguments

build_Adaboost_model <- function(SE, Signature, rmBE = TRUE, response_NR = TRUE, PT_drop, ...){
  data <- dataProcess(SE, Signature, rmBE, response_NR, FALSE)
  if(PT_drop)
    data <- PT_filter(data)
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
#' @param response_NR If TRUE, classify patients with CR, MR, PR as Responders (R), and those with PD, SD, NR as Non-Responders(NR).
#' @param PT_drop If TRUE, only Untreated patient will be use for model training.
#' @param ... the arguments

build_Logitboost_model <- function(SE, Signature, rmBE = TRUE, response_NR = TRUE, PT_drop, ...){
  v_Args <- list(...)
  nIter <- ifelse("nIter" %in% names(v_Args),v_Args[["nIter"]],150)
  data <- dataProcess(SE, Signature, rmBE, response_NR, FALSE)
  if(PT_drop)
    data <- PT_filter(data)
  model <- LogitBoost(xlearn = t(data[[1]]), ylearn = factor(data[[2]]$response), nIter = nIter)
}


#' @title Build Logistics prediction model for immunotherapy response
#' @description Generate a Logistics model.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR If TRUE, classify patients with CR, MR, PR as Responders (R), and those with PD, SD, NR as Non-Responders(NR).
#' @param PT_drop If TRUE, only Untreated patient will be use for model training.
#' @param ... the arguments

build_Logistics_model <- function(SE, Signature, rmBE = FALSE, response_NR = TRUE, PT_drop, ...){
  data <- dataProcess(SE, Signature, rmBE, response_NR, FALSE)
  if(PT_drop)
    data <- PT_filter(data)
  df <- data.frame(response=ifelse(data[[2]]$response=='R',1,0),t(data[[1]]),check.names = FALSE)

  v_Args <- list(...)
  Args <- c(list(formula=response ~.,
                 data=df,
                 family=stats::binomial(link = "logit"),
                 control=ifelse("control" %in% names(v_Args),v_Args[["control"]],list(maxit=100))),
            v_Args[names(v_Args) != "control"])
  model <- do.call(stats::glm, Args)
}


#' @title Build Lasso-cox prediction model for Survival
#' @description Generate a Lasso-cox model.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param PT_drop If TRUE, only Untreated patient will be use for model training.
#' @param lambda the lambda value for Lasso.
#' @param ... the arguments

build_SURV_Model <- function(SE, Signature, rmBE = FALSE, PT_drop, lambda = 1,...){
  data <- dataProcess(SE, Signature, rmBE, FALSE, FALSE)
  if(PT_drop)
    data <- PT_filter(data)

  time <- as.numeric(data[[2]]$overall.survival..days.)
  status <- sub('Dead','1', data[[2]]$vital.status) %>%
    sub('Alive','0',.) %>%
    as.numeric()

  df <- na.omit(data.frame(time, status,t(data[[1]])))
  model <-
  penalized::penalized(Surv(time, status) ~ .,
                       data = df,
                       model = "cox",
                       lambda1 = lambda)
  model@penalized <- model@penalized[model@penalized!=0]
  return(model)
}


#' @title Post Treatment filter
#' @description Post Treatment filter
#' @param data the data

PT_filter <- function(data){
  idx_UT <- which(data[[2]]$Treatment == "PRE")
  if(length(idx_UT) == 0)
    stop("All patients in data set have been treated. Setting the parameter PT_drop to FALSE to run anyway.")
  data[[1]] <- data[[1]][,idx_UT,drop=FALSE]
  data[[2]] <- data[[2]][idx_UT,,drop=FALSE]
  return(data)
}
