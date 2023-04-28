#' @title judge whether all the elements in vector V are NA
#' @param V must be a vector.

is.NA_vec <- function(V){
  return(all(is.na(V)))
}

#' @title turn the 0 elements in vector V to NA
#' @param V must be a vector.

zero2na <- function(V){
  if(all(is.na(V)))
    return(V)
  if(any(is.na(V)) && all(V[!is.na(V)] == 0)){
    V[!is.na(V)] <- NA
    return(V)
  }else
    if(all(V == 0))
      return(rep(NA, length(V)))
  return(V)
}

#' @title Ranking features in matrix with Gini index
#' @description By calculating the Gini index of different genes, you can get an overview of the classification efficiency of different genes.
#' @param mtr An gene expression matrix which rownames are gene symbols and colnames are sample ID.
#' @param label The classification labels of the samples.
#' @param ascending If ascending = TRUE, the result will be display in ascending order.
#' @export

Gini_rank <- function(mtr, label, ascending = TRUE){
  features_Gini <- apply(mtr, 1, Gini, label = label)
  features_Rank <- names(sort(features_Gini))
  if(ascending == FALSE)
    features_Rank <- rev(features_Rank)
  return(features_Rank)
}

#' @title Ranking features in vector with Gini coefficient
#' @param vec An vector. Usually is one row of an matrix.
#' @param label The classification label of matrix.
#'

Gini <- function(vec, label){
  NR <- grep('NR', label)
  R <- (1:length(label))[!1:length(label) %in% NR]
  Gini_H <- 1 - (length(grep('TRUE',vec[R] == 'HIGH'))/length(vec[vec == 'HIGH'])) ^ 2 - (1 - length(grep('TRUE',vec[R] == 'HIGH'))/length(vec[vec == 'HIGH']))^2
  Gini_L <- 1 - (length(grep('TRUE',vec[R] == 'LOW'))/length(vec[vec == 'LOW'])) ^ 2 - (1 - length(grep('TRUE',vec[R] == 'LOW'))/length(vec[vec == 'LOW']))^2
  Gini_Gene <- (length(vec[vec == 'HIGH'])*Gini_H + length(vec[vec == 'LOW'])*Gini_L)/length(label)
  return(Gini_Gene)
}

#' @title Binding expression matrices from data folder in tigeR together
#' @description Extract expression data in particular data set or data sets from the data folder in tigeR. If there are more than one data set, this function will return an matrix which binds all the expression matrices by column.
#' @param datasetNames the name of data set or data sets you want to use.
#' @export
#'

extract_mtr <- function(datasetNames){
  for (name in datasetNames) {
    #browser()
    if(!exists('inteMatrix', envir = current_env())){
      data(list = name, envir = current_env(), overwrite = TRUE)
      exp <- get(name)
      exp_mtr <- exp[,-1]
      exp_mtr <- as.matrix(exp_mtr)
      rownames(exp_mtr) <- exp[,1]

      inteMatrix <- exp_mtr
      if(length(datasetNames > 1))
        next
    }
    data(list = name, envir = current_env(), overwrite = TRUE)
    exp <- get(name)
    exp_mtr <- exp[,-1]
    exp_mtr <- as.matrix(exp_mtr)
    rownames(exp_mtr) <- exp[,1]

    inteMatrix <- cbind(inteMatrix, exp_mtr)
  }
  return(inteMatrix)
}


#' @title Binding response data from data folder in tigeR together
#' @description Extract response data in particular data set or data sets from the data folder in tigeR. If there are more than one data set, this function will return an vector which contains the response data of every data sets.
#' @param datasetNames the name of data set or data sets you want to use.
#' @export
#'

extract_label <-function(datasetNames){
  for (name in datasetNames) {
    if(!exists('inteVector', envir = current_env())){
      data(list = name, envir = current_env(), overwrite = TRUE)
      meta <- get(name)
      meta_mtr <- meta$response

      inteVector <- meta_mtr
      if(length(datasetNames > 1))
        next
    }
    data(list = paste0name, envir = current_env(), overwrite = TRUE)
    meta <- get(name)
    meta_mtr <- meta$response
    inteVector <- c(inteVector, meta_mtr)
  }

  inteVector <- sub('CR|MR|PR|SD|CRPR', 'R', inteVector)
  inteVector <- sub('PD', 'NR', inteVector)
  inteVector[inteVector == 'UNK'] <- NA
  return(inteVector)
}

#' @title standardization of response labels
#' @description turn label to 'R' or 'NR'
#' @param V an vector
#' @export

response_standardize <- function(V){
  V <- sub('CR|MR|PR|CRPR', 'R', V)
  V <- sub('PD|SD', 'NR', V)
  return(V)
}


#' @title Max_Min normalization.
#' @description (x - min(x))/(max(x) - min(x))
#' @param exp_mtr an matrix which rows represent genes columns represent samples.
#' @export

max_min_normalization <- function(exp_mtr){
  mini <- apply(exp_mtr, 1, min)
  maxi <- apply(exp_mtr, 1, max)
  interval <- maxi - mini
  exp_mtr <- (exp_mtr - mini) / interval
  return(exp_mtr)
}


#' @title perform naive bayes prediction model.
#' @description Generate a naive bayes model.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param isList whether SE is list
#' @importFrom SummarizedExperiment assay
#' @importFrom magrittr %>%

bind_mtr <- function(SE,isList){
  if (!isList){
    exp_mtr <- assay(SE)
  } else if (isList){
    lapply(SE, assay) %>% unlist() %>% matrix(nrow = nrow(assay(SE[[1]]))) -> exp_mtr
    assay(SE[[1]]) %>% rownames() ->rownames(exp_mtr)
    lapply(SE,colnames) %>% unlist() -> colnames(exp_mtr)
  }
  return(exp_mtr)
}


#' @title perform naive bayes prediction model.
#' @description Generate a naive bayes model.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param isList whether SE is list
#' @importFrom SummarizedExperiment assay
#' @importFrom magrittr %>%

bind_meta <- function(SE,isList){
  if (!isList){
    meta <- as.data.frame(SE@colData)
  } else if (isList){
    meta <- as.data.frame(SE[[1]]@colData)
    meta$batch <- rep('batch1',nrow(meta))
    for (i in 2:length(SE)) {
      paste0('batch', i) %>% rep(nrow(SE[[i]]@colData)) -> batch
      SE[[i]]@colData %>% as.data.frame() %>% cbind(batch) %>% rbind(meta,.) -> meta
    }
  }
  return(meta)
}


#' @title Remove batch effect.
#' @description Generate a naive bayes model.
#' @param mtr an expression matrix.
#' @param meta meta informations.

rmBE <- function(mtr, meta){
  model <- model.matrix(~as.factor(meta$response))
  mtr <- dataPreprocess(mtr, rownames(mtr), turn2HL = FALSE)
  mtr <- sva::ComBat(dat = mtr,batch = as.factor(meta$batch),mod = model)
  return(mtr)
}


#' @title Filting missing response value.
#' @description Generate a naive bayes model.
#' @param response a vector which contains response information.
#' @importFrom magrittr %>%

response_filter <- function(response){
  idx <- grep('NE|UNK',response)
  return(idx)
}

