#' @title judge whether all the elements in vector V are NA
#'

is.NA_vec <- function(V){
  return(all(is.na(V)))
}

#' @title
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
#'

Gini_rank <- function(mtr, label, ascending = TRUE){
  features_Gini <- apply(final_Expr, 1, Gini)
  features_Rank <- names(sort(features_Gini))
  if(ascending == FALSE)
    features_Rank <- rev(features_Rank)
  return(features_Rank)
}

#' @title Ranking features in vector with Gini index
#'

Gini <- function(vec){
  NR <- grep('NR', label)
  R <- (1:length(label))[!1:length(label) %in% NR]
  Gini_H <- 1 - (length(grep('TRUE',vec[R] == 'HIGH'))/length(vec[vec == 'HIGH'])) ^ 2 - (1 - length(grep('TRUE',vec[R] == 'HIGH'))/length(vec[vec == 'HIGH']))^2
  Gini_L <- 1 - (length(grep('TRUE',vec[R] == 'LOW'))/length(vec[vec == 'LOW'])) ^ 2 - (1 - length(grep('TRUE',vec[R] == 'LOW'))/length(vec[vec == 'LOW']))^2
  Gini_Gene <- (length(vec[vec == 'HIGH'])*Gini_H + length(vec[vec == 'LOW'])*Gini_L)/length(label)
  return(Gini_Gene)
}

#' @title Binding expression matrixes from data folder together
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


#' @title extract and modify response data for naive bayes model
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
