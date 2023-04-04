#' @title extract an specific gene subset of the expression matrix
#' @description The function will return an expression matrix which only cotains the expression of user designated genes.
#' @param exp_mtr An expression matrix which rownames are gene symbol and colnames are sample ID.
#' @param Signature The aiming gene set(only Gene SYMBOL).
#' @param turn2HL Whether turn numeric data to 'HIGH' and LOW".
#' @export
#'

dataPreprocess <- function(exp_mtr, Signature, turn2HL = TRUE){
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

#' @title perform naive bayes prediction model.
#' @description Generate a naive bayes model.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @import e1071
#' @import sva
#' @export

build_NB_model <- function(SE, Signature, rmBE = TRUE){
  if (!is.list(SE)){
    if (is.numeric(SummarizedExperiment::assay(SE))){
      exp_mtr <- dataPreprocess(SummarizedExperiment::assay(SE), Signature)
      response <- SE$response
    } else{
      stop("The assay must be numeric!")
    }
  } else if (is.list(SE)){
    if (all(lapply(lapply(SE, SummarizedExperiment::assay), is.numeric) == TRUE)){
      batch_count <- unlist(lapply(SE, ncol))

      batch <- c()
      response <- c()
      for (i in 1:length(batch_count)) {
        batch <- c(batch, rep(paste0('batch', i), batch_count[i]))
        response <- c(response,SE[[i]]$response)
      }
      Expr <- matrix(unlist(lapply(SE, SummarizedExperiment::assay)), nrow = nrow(SummarizedExperiment::assay(SE[[1]])))
      rownames(Expr) <- rownames(SummarizedExperiment::assay(SE[[1]]))
      if(rmBE){
        model <- model.matrix(~as.factor(response))
        Expr <- dataPreprocess(Expr, rownames(Expr), turn2HL = FALSE)
        inte_Expr <- sva::ComBat(dat = Expr,batch = as.factor(batch),mod = model)
      } else {
        inte_Expr <- Expr
      }
      exp_mtr <- dataPreprocess(inte_Expr, Signature)
    } else{
      stop("The matrices in list must be numeric!")
    }
  } else{
    stop("Parameter 'exp' must be matrix or list!")
  }

  model <- e1071::naiveBayes(t(exp_mtr), response, laplace = 1)
  return(model)
}


#' @title perform SVM prediction model.
#' @description Generate a Support Vector Machine model.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @import e1071
#' @import sva
#' @export

build_SVM_model <- function(SE, Signature, rmBE = TRUE){
  if (!is.list(SE)){
    if (is.numeric(SummarizedExperiment::assay(SE))){
      exp_mtr <- dataPreprocess(SummarizedExperiment::assay(SE), Signature, turn2HL = FALSE)
      response <- SE$response
    } else{
      stop("The assay must be numeric!")
    }
  } else if (is.list(SE)){
    if (all(lapply(lapply(SE, SummarizedExperiment::assay), is.numeric) == TRUE)){
      batch_count <- unlist(lapply(SE, ncol))

      batch <- c()
      response <- c()
      for (i in 1:length(batch_count)) {
        batch <- c(batch, rep(paste0('batch', i), batch_count[i]))
        response <- c(response,SE[[i]]$response)
      }
      Expr <- matrix(unlist(lapply(SE, SummarizedExperiment::assay)), nrow = nrow(SummarizedExperiment::assay(SE[[1]])))
      rownames(Expr) <- rownames(SummarizedExperiment::assay(SE[[1]]))
      if(rmBE){
        model <- model.matrix(~as.factor(response))
        inte_Expr <- sva::ComBat(dat = Expr,batch = as.factor(batch),mod = model)
      } else {
        inte_Expr <- Expr
      }
      exp_mtr <- dataPreprocess(inte_Expr, Signature, turn2HL = FALSE)
    } else{
      stop("The matrices in list must be numeric!")
    }
  } else{
    stop("Parameter 'exp' must be matrix or list!")
  }

  model <- e1071::svm(x = t(na.omit(exp_mtr)),
                      y = as.factor(response),
                      scale = TRUE,
                      type = 'C',
                      kernel = 'linear',
                      probability = TRUE)
  return(model)
}

#' @title perform ramdom forest prediction model.
#' @description Generate a random forest model.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @import randomForest
#' @import sva
#' @export

build_RF_model <- function(SE, Signature, rmBE = TRUE){
  if (!is.list(SE)){
    if (is.numeric(SummarizedExperiment::assay(SE))){
      exp_mtr <- dataPreprocess(SummarizedExperiment::assay(SE), Signature, turn2HL = FALSE)
      response <- SE$response
    } else{
      stop("The assay must be numeric!")
    }
  } else if (is.list(SE)){
    if (all(lapply(lapply(SE, SummarizedExperiment::assay), is.numeric) == TRUE)){
      batch_count <- unlist(lapply(SE, ncol))

      batch <- c()
      response <- c()
      for (i in 1:length(batch_count)) {
        batch <- c(batch, rep(paste0('batch', i), batch_count[i]))
        response <- c(response,SE[[i]]$response)
      }
      Expr <- matrix(unlist(lapply(SE, SummarizedExperiment::assay)), nrow = nrow(SummarizedExperiment::assay(SE[[1]])))
      rownames(Expr) <- rownames(SummarizedExperiment::assay(SE[[1]]))
      if(rmBE){
        model <- model.matrix(~as.factor(response))
        inte_Expr <- sva::ComBat(dat = Expr,batch = as.factor(batch),mod = model)
      } else {
        inte_Expr <- Expr
      }
      exp_mtr <- dataPreprocess(inte_Expr, Signature, turn2HL = FALSE)
    } else{
      stop("The matrices in list must be numeric!")
    }
  } else{
    stop("Parameter 'exp' must be matrix or list!")
  }

  model <- randomForest::randomForest(x = t(na.omit(exp_mtr)),
                                      y = as.factor(response),
                                      ntree = 150,
                                      mtry = 9,
                                      cutoff = c(0.80, 0.20))
  return(model)
}

#' @title perform cancerclass prediction model.
#' @description Generate a cancerclass model.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @importFrom cancerclass fit
#' @import sva
#' @export

build_CC_model <- function(SE, Signature, rmBE = TRUE){
  if (!is.list(SE)){
    if (is.numeric(SummarizedExperiment::assay(SE))){
      exp_mtr <- dataPreprocess(SummarizedExperiment::assay(SE), Signature, turn2HL = FALSE)
      response <- SE$response
    } else{
      stop("The assay must be numeric!")
    }
  } else if (is.list(SE)){
    if (all(lapply(lapply(SE, SummarizedExperiment::assay), is.numeric) == TRUE)){
      batch_count <- unlist(lapply(SE, ncol))

      batch <- c()
      response <- c()
      for (i in 1:length(batch_count)) {
        batch <- c(batch, rep(paste0('batch', i), batch_count[i]))
        response <- c(response,SE[[i]]$response)
      }
      Expr <- matrix(unlist(lapply(SE, SummarizedExperiment::assay)), nrow = nrow(SummarizedExperiment::assay(SE[[1]])))
      rownames(Expr) <- rownames(SummarizedExperiment::assay(SE[[1]]))
      if(rmBE){
        model <- model.matrix(~as.factor(response))
        inte_Expr <- sva::ComBat(dat = Expr,batch = as.factor(batch),mod = model)
      } else {
        inte_Expr <- Expr
      }
      exp_mtr <- dataPreprocess(inte_Expr, Signature, turn2HL = FALSE)
    } else{
      stop("The matrices in list must be numeric!")
    }
  } else{
    stop("Parameter 'exp' must be matrix or list!")
  }

  pData <- data.frame(class = response, row.names = colnames(exp_mtr))
  metadata <- data.frame(labelDescription = colnames(pData), row.names = colnames(pData))
  adf <- new("AnnotatedDataFrame", data = pData, varMetadata = metadata)
  exampleSet <- methods::new("ExpressionSet", exprs = exp_mtr, phenoData = adf)
  model <- cancerclass::fit(exampleSet, method = "welch.test")

  return(model)
}


#' @title perform Adaboost prediction model.
#' @description Generate a Adaboost model.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @importFrom adabag boosting
#' @import sva
#' @export

build_Adaboost_model <- function(SE, Signature, rmBE = TRUE){
  if (!is.list(SE)){
    if (is.numeric(SummarizedExperiment::assay(SE))){
      exp_mtr <- dataPreprocess(SummarizedExperiment::assay(SE), Signature, turn2HL = FALSE)
      response <- SE$response
    } else{
      stop("The assay must be numeric!")
    }
  } else if (is.list(SE)){
    if (all(lapply(lapply(SE, SummarizedExperiment::assay), is.numeric) == TRUE)){
      batch_count <- unlist(lapply(SE, ncol))

      batch <- c()
      response <- c()
      for (i in 1:length(batch_count)) {
        batch <- c(batch, rep(paste0('batch', i), batch_count[i]))
        response <- c(response,SE[[i]]$response)
      }
      Expr <- matrix(unlist(lapply(SE, SummarizedExperiment::assay)), nrow = nrow(SummarizedExperiment::assay(SE[[1]])))
      rownames(Expr) <- rownames(SummarizedExperiment::assay(SE[[1]]))
      if(rmBE){
        model <- model.matrix(~as.factor(response))
        inte_Expr <- sva::ComBat(dat = Expr,batch = as.factor(batch),mod = model)
      } else {
        inte_Expr <- Expr
      }
      exp_mtr <- dataPreprocess(inte_Expr, Signature, turn2HL = FALSE)
    } else{
      stop("The matrices in list must be numeric!")
    }
  } else{
    stop("Parameter 'exp' must be matrix or list!")
  }

  df <- data.frame(class = response, t(exp_mtr))
  df$class <- factor(df$class)
  model <- adabag::boosting(class ~ ., data = df, mfinal = 10, boos = TRUE)

  return(model)
}


#' @title perform LogitBoost prediction model.
#' @description Generate a LogitBoost model.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @importFrom caTools LogitBoost
#' @import sva
#' @export

build_Logitboost_model <- function(SE, Signature, rmBE = TRUE){
  if (!is.list(SE)){
    if (is.numeric(SummarizedExperiment::assay(SE))){
      exp_mtr <- dataPreprocess(SummarizedExperiment::assay(SE), Signature, turn2HL = FALSE)
      response <- SE$response
    } else{
      stop("The assay must be numeric!")
    }
  } else if (is.list(SE)){
    if (all(lapply(lapply(SE, SummarizedExperiment::assay), is.numeric) == TRUE)){
      batch_count <- unlist(lapply(SE, ncol))

      batch <- c()
      response <- c()
      for (i in 1:length(batch_count)) {
        batch <- c(batch, rep(paste0('batch', i), batch_count[i]))
        response <- c(response,SE[[i]]$response)
      }
      Expr <- matrix(unlist(lapply(SE, SummarizedExperiment::assay)), nrow = nrow(SummarizedExperiment::assay(SE[[1]])))
      rownames(Expr) <- rownames(SummarizedExperiment::assay(SE[[1]]))
      if(rmBE){
        model <- model.matrix(~as.factor(response))
        inte_Expr <- sva::ComBat(dat = Expr,batch = as.factor(batch),mod = model)
      } else {
        inte_Expr <- Expr
      }
      exp_mtr <- dataPreprocess(inte_Expr, Signature, turn2HL = FALSE)
    } else{
      stop("The matrices in list must be numeric!")
    }
  } else{
    stop("Parameter 'exp' must be matrix or list!")
  }

  model <- caTools::LogitBoost(xlearn = t(exp_mtr), ylearn = factor(response), nIter = 300)

  return(model)
}


