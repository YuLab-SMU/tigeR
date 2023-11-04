#' @title Assessment of gene set
#' @description Perform differential expression analysis and survival analysis in certain gene and return the result.
#' @param SE The dataset you want to use to test your Signature.A SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain treatment information names Treatment.
#' @param Signature The gene which you wanted.
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR If TRUE, only use R or NR to represent Immunotherapy response of patients.
#' @export

Signature_assessment <- function(SE, Signature, rmBE, response_NR){
  if(is.character(Signature)){
    names(Signature) <- Signature
    Signature <- rep(1, length(Signature))
  }

  data <- dataProcess(SE, names(Signature), rmBE, response_NR, FALSE)
  value <- weight_mean_signature(data[[1]], Signature)
  ROC <- roc(data[[2]]$response, value)
}
