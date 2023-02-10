#' @title judge whether all the elements in vector V are NA
#'

is.NA_vec <- function(V){
  return(all(is.na(V)))
}

#' @title turn the 0 elements in vector V to NA
#'

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
  features_Gini <- apply(final_Expr, 1, Gini)
  features_Rank <- names(sort(features_Gini))
  if(ascending == FALSE)
    features_Rank <- rev(features_Rank)
  return(features_Rank)
}

#' @title Ranking features in vector with Gini coefficient
#'

Gini <- function(vec){
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
      data(list = paste0(name, '.Response'), envir = current_env())
      exp_mtr <- exp[,-1]
      exp_mtr <- as.matrix(exp_mtr)
      rownames(exp_mtr) <- exp[,1]

      inteMatrix <- exp_mtr
      if(length(datasetNames > 1))
        next
    }
    data(list = paste0(name, '.Response'), envir = current_env())
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
      data(list = paste0(name, '.meta'), envir = current_env())
      meta_mtr <- meta$response

      inteVector <- meta_mtr
      if(length(datasetNames > 1))
        next
    }
    data(list = paste0(name, '.meta'), envir = current_env())
    meta_mtr <- meta$response
    inteVector <- c(inteVector, meta_mtr)
  }

  inteVector <- sub('CR|MR|PR|SD|CRPR', 'R', inteVector)
  inteVector <- sub('PD', 'NR', inteVector)
  inteVector[inteVector == 'UNK'] <- NA
  return(inteVector)
}
