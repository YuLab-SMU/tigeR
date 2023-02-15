#' @title extract an specific gene subset of the expression matrix
#' @description The function will return an expression matrix which only cotains the expression of user designated genes.
#' @param exp_mtr An expression matrix which rownames are gene symbol and colnames are sample ID.
#' @param Signature The aiming gene set(only Gene SYMBOL).
#' @export
#'

dataPreprocess <- function(exp_mtr, Signature){
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

  return(exp_mtr)
}

#' @title perform naive bayes prediction model.
#' @description Generate a naive bayes model.
#' @param exp a matrix or list consist of matrix. The rownames of every matrix must be GENE SYMBOL. If it is a list, the row num and row names must be the same.
#' @param response an character vector corresponding to samples' clinical status (respond or not,survival or dead etc.)
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @import e1071
#' @import sva
#' @export

predict_Response <- function(exp, response, Signature, rmBE = TRUE){
  if (!all(response %in% c('R', 'NR'))){
    stop("Parameter 'response' must be an character vector which only contains'R' or 'NR'!")
  }

  if (is.matrix(exp)){
    if (is.numeric(exp)){
      exp_mtr <- dataPreprocess(exp, Signature)
    } else{
      stop("The data type must be numeric!")
    }
  } else if (is.list(exp)){
    if (is.numeric(unlist(exp))){
      batch_count <- unlist(lapply(exp, ncol))

      batch <- c()
      for (i in 1:length(batch_count)) {
        batch <- c(batch, rep(paste0('batch', i), batch_count[i]))
      }
      Expr <- matrix(unlist(exp), nrow = nrow(exp[[1]]))
      rownames(Expr) <- rownames(exp[[1]])
      model <- model.matrix(~as.factor(response))
      combat_Expr <- sva::ComBat(dat = Expr,batch = as.factor(batch),mod = model)
      exp_mtr <- dataPreprocess(combat_Expr, Signature)
    } else{
      stop("The matrices in list must be numeric!")
    }
  } else{
    stop("Parameter 'exp' must be matrix or list!")
  }

  model <- e1071::naiveBayes(t(exp_mtr), response, laplace = 0)
  return(model)
}


























