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



