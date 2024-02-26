#' @title Calculating Signature score of existing immunotherapy prognosis Signature.
#' @description By employing the Signature_calculation() function, you can obtain a comprehensive signature score matrix for the 23 signatures in TigeR. In this matrix, the columns represent the signature scores, and the rows denote the sample names.
#' @param SE a SummarizedExperiment object for which you want to calculate the Signature Score.
#' @param exp_mtr an expression matrix for which you want to calculate the Signature Score.
#' @export
#'

Signature_calculation <- function(SE=NULL,exp_mtr=NULL){
  if(!missing(SE)){
    isList <- is.list(SE)
    exp_mtr <- bind_mtr(SE, isList)
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
  return(df)
}


#' @title Giving immunotherapy prognosis using IRS.
#' @description The function will return a vector calculated using IRS.
#' @param exp_mtr An expression matrix which rownames are gene symbol and colnames are sample ID.
#' @export
#'

IRS_grading <- function(exp_mtr){
  IRS <- c('EGFR','MAPK1','TFRC','IRF1','ADIPOR2','GBP2','CTSS','THBS1','GBP2','CCN2','PSMD11','SRC','KIR2DL4','NOX4','MAP2K1','ELAVL1')
  Expr <- dataPreprocess(exp_mtr, IRS, turn2HL = FALSE)
  rownames_identify <- rownames(Expr)
  for (gene in IRS) {
    if(!(gene %in% rownames(Expr))){
      vec <- rep(0, ncol(Expr))
      Expr <- rbind(Expr, vec)
    }
  }


  idx <- which(!IRS %in% rownames(Expr))
  if(length(idx) != 0)
    rownames(Expr) <- c(rownames_identify, IRS[idx])
  if(nrow(Expr) > 1){
    result  <- 0.599 * ifelse(Expr[rownames(Expr) == IRS[1],]>Expr[rownames(Expr) == IRS[2],],1,0) +
      0.579 * ifelse(Expr[rownames(Expr) == IRS[3],]>Expr[rownames(Expr) == IRS[4],],1,0) +
      0.376 * ifelse(Expr[rownames(Expr) == IRS[5],]>Expr[rownames(Expr) == IRS[6],],1,0) -
      0.305 * ifelse(Expr[rownames(Expr) == IRS[7],]>Expr[rownames(Expr) == IRS[8],],1,0) -
      0.472 * ifelse(Expr[rownames(Expr) == IRS[9],]>Expr[rownames(Expr) == IRS[10],],1,0) +
      0.469 * ifelse(Expr[rownames(Expr) == IRS[11],]>Expr[rownames(Expr) == IRS[12],],1,0) -
      0.645 * ifelse(Expr[rownames(Expr) == IRS[13],]>Expr[rownames(Expr) == IRS[14],],1,0) +
      0.495 * ifelse(Expr[rownames(Expr) == IRS[15],]>Expr[rownames(Expr) == IRS[16],],1,0)
  } else{
    warning("There absence in IRS genes. Please check your data!")
    return()
  }
  return(result)
}


#' @title Giving immunotherapy prognosis using tGE8.
#' @description The function will return a vector calculated using tGE8.
#' @param exp_mtr An expression matrix which rownames are gene symbol and colnames are sample ID.
#' @importFrom stats sd
#' @importFrom stats median
#' @export

tGE8_grading <- function(exp_mtr){
  tGE8 <- c('IFNG','CXCL9','CD8A','GZMA','GZMB','CXCL10','PRF1','TBX21')
  Expr <- dataPreprocess(exp_mtr, tGE8, turn2HL = FALSE)
  if(nrow(Expr) == 8){
    Expr <- dataPreprocess(exp_mtr, tGE8, turn2HL = FALSE)
    average <- apply(Expr, 1, mean)
    standard_error <- apply(Expr, 1, sd)
    ZScore <- (Expr - average) / standard_error
    result <- apply(ZScore, 2, median)
    names(result) <- colnames(ZScore)
  } else{
    warning("There absence in tGE8 genes. Please check your data!")
    return()
  }
  return(result)
}
