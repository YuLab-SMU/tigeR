#' @title Assessment of gene set
#' @description Perform differential expression analysis and survival analysis in certain gene and return the result.
#' @param SE the dataset you wish to use to test your Signature. A SummarizedExperiment (SE) object, which can be either a single SE object or a list of SE objects. Note that for each SE object, the colData must contain treatment information under the column name Treatment.
#' @param Signature The gene which you wanted.
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR If TRUE, only use R or NR to represent Immunotherapy response of patients.
#' @export

Signature_assessment <- function(SE, Signature, rmBE, response_NR){
  if(is.character(Signature)){
    names(Signature) <- Signature
    Signature[] <- rep(1, length(Signature))
  }

  data <- dataProcess(SE, names(Signature), rmBE, response_NR, FALSE)
  value <- weight_mean_signature(data[[1]], Signature)
  ROC <- roc(data[[2]]$response, value)
}

