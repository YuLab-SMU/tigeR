#' @title Assessing the performance of signature on predicting immunotherapy response
#' @description generate a Receiver Operating Characteristic (ROC) object and a curve to assess the predictive performance.
#' @param SE the dataset you use to assess your signature. A SummarizedExperiment (SE) object, which can be either a single SE object or a list of SE objects. Note that for each SE object, the colData must contain treatment information under the column named Treatment.
#' @param Signature a gene vector represents the built-in or user-defined signature.
#' @param method the method for calculating gene set scores which has several options: "Average_mean", "Weighted_mean", or "GSVA". The method can be set to NULL if the length of the parameter geneSet is 1. This means that if you are working with only one gene, the specific calculation method may not be applicable or necessary.
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR If TRUE, classify patients with CR, MR, PR as Responders (R), and those with PD, SD, NR as Non-Responders(NR).
#' @param PT_drop If TRUE, only untreated patient will be use for model training.
#' @param auc.pos the position of the AUC value
#' @param auc.round the decimal places you want to keep for auc value
#' @param textcol the color of the text in the plot.
#' @param panelcol the color of the panel border and ticks in the plot.
#' @return
#'   \describe{
#'   Return a list contain the following elements which including:
#'     \item{\code{ROC}}{a receiver operating characteristic (ROC) object.}
#'     \item{\code{figure}}{a roc curve to show the predictive performance.}
#'   }
#' @examples
#' sig_roc <-
#' roc_biomk(MEL_GSE93157,
#'          Weighted_mean_Sigs$Tcell_inflamed_GEP,
#'          method = "Weighted_mean",
#'          rmBE=TRUE,
#'          response_NR=TRUE)
#' sig_roc
#' @export

roc_biomk <- function(SE, Signature, method = NULL, rmBE=FALSE, response_NR=TRUE, PT_drop=TRUE, auc.pos=c(0.3,0.42), auc.round=3, textcol="black",panelcol="black"){
  data <- dataProcess(SE, names(Signature), rmBE, response_NR, FALSE)
  if(PT_drop)
    data <- PT_filter(data)
  value <- Core(data[[1]], Signature, method)
  ROC <- pROC::roc(data[[2]]$response, value)
  figure <- plt_roc(ROC,auc.pos=auc.pos,auc.round=auc.round,textcol=textcol,panelcol=panelcol)

  return(list(ROC,figure))
}

