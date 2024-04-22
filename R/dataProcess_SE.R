#' @title Process data before running machine learning algorithm
#' @description Process data before running machine learning algorithm
#' @param SE a SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR If TRUE, only use R or NR to represent Immunotherapy response of patients.
#' @param turn2HL If TRUE, the expression value of a gene is divided to "HIGH" or "LOW" based on its median expression.
#' @export

dataProcess_SE <- function(SE, Signature, rmBE, response_NR, turn2HL){
  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)
  meta <- bind_meta(SE, isList)

  if(response_NR)
    meta$response <- response_standardize(meta$response)

  if(rmBE){
    exp_mtr <- rmBE(exp_mtr,meta)
    colnames(exp_mtr) <- rownames(meta)
  }

  idx <- response_filter(meta$response)
  idx_all <- seq_along(meta[,1])
  if(length(idx)==0){
    idx <- idx_all
  }else{
    idx <- idx_all[!idx_all %in% idx]
  }

  f <- dataPreprocess(exp_mtr, Signature, turn2HL, meta)
  if(turn2HL){
    exp_mtr <- f[[1]][,idx]
  }else{
    exp_mtr <- f[,idx,drop=FALSE]
  }
  meta <- meta[idx,]

  absent <- meta$response_NR=="UNK"

  if(turn2HL){
    thres <- f[[2]]
    return(list(exp_mtr[,!absent,drop=FALSE],meta[!absent,,drop=FALSE],thres))
  }

  SE_final <-
    SummarizedExperiment::SummarizedExperiment(assays=S4Vectors::SimpleList(exp_mtr[,!absent,drop=FALSE]),
                                               colData=S4Vectors::DataFrame(meta[!absent,,drop=FALSE]),
                                               checkDimnames=TRUE)
  return(SE_final)
}
