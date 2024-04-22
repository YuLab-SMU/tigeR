#' @title Calculating Signature score of existing immunotherapy response Signature.
#' @description By employing the Signature_calculation() function, you can obtain a comprehensive signature score matrix for the 23 signatures in tigeR. In this matrix, the columns represent the signature scores, and the rows denote the sample names.
#' @param SE a SummarizedExperiment object for which you want to calculate the Signature Score.
#' @param exp_mtr an expression matrix for which you want to calculate the Signature Score.
#' @param meta meta data of samples
#' @param Signature a genes vector represents user-defined signature for Immunotherapy response. If NULL, the function will only calculate 23 built-in signatures in tigeR.
#' @param method the method for calculating gene set scores which has several options: Average_mean, Weighted_mean, or GSVA. The method can be set to NULL if the length of the parameter geneSet is 1. This means that if you are working with only one gene, the specific calculation method may not be applicable or necessary.
#' @param PT_drop If TRUE, only Untreated patient will be use for model training.
#' @export

score_biomk_SE <- function(SE=NULL, exp_mtr=NULL, meta=NULL, Signature=NULL, method="Average_mean",PT_drop=TRUE){
  if(!missing(SE)){
    isList <- is.list(SE)
    exp_mtr <- bind_mtr(SE, isList)
    meta <- bind_meta(SE, isList)
  }

  if(PT_drop){
    idx_UT <- which(meta$Treatment == 'PRE')
    if(length(idx_UT) == 0)
      stop("All patients are treated!")
    exp_mtr <- exp_mtr[,idx_UT,drop=FALSE]
    meta <- meta[idx_UT,,drop=FALSE]
  }

  Average_mean_Sigs <- NULL
  Weighted_mean_Sigs <- NULL
  ZScore_PCA_Sigs <- NULL
  data(Average_mean_Sigs, package = 'tigeR', envir = current_env())
  data(Weighted_mean_Sigs, package = 'tigeR', envir = current_env())
  data(ZScore_PCA_Sigs, package = 'tigeR', envir = current_env())
  df <- data.frame(IRS=IRS_grading(exp_mtr),tGE8=tGE8_grading(exp_mtr))
  for (i in Average_mean_Sigs) {
    df <- cbind(df, average_mean_signature(exp_mtr, i))
  }
  for (i in Weighted_mean_Sigs) {
    df <- cbind(df,weight_mean_signature(exp_mtr,i))
  }
  for (i in ZScore_PCA_Sigs) {
    df <- cbind(df,ZScore_PCA_signature(exp_mtr,i))
  }
  colnames(df) <- c('IRS','tGE8',
                    names(Average_mean_Sigs),
                    names(Weighted_mean_Sigs),
                    names(ZScore_PCA_Sigs))
  if(!is.null(Signature)){
    sig <- Core(exp_mtr, Signature, method)
    df <- cbind(`Customed Signature`=sig,df)
  }
  SummarizedExperiment(assays=S4Vectors::SimpleList(t(df)),
                       colData=S4Vectors::DataFrame(colData(SE)),
                       checkDimnames=TRUE)
}
