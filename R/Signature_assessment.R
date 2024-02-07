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
  ROC <- pROC::roc(data[[2]]$response, value)
  figure <- pROC::ggroc(ROC, color = "black", size = 1) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = 1,
                                       y = 1, yend = 0), color = "#646464", size = 0.5,
                          linetype = "solid") +
    ggplot2::annotate("text",x = 0.3, y = 0.42, label = paste0("AUC=", round(ROC$auc, digits = 4)), size = 4.5,color = "#646464") +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = element_blank(),
                   plot.title = element_text(face = "bold",
                                             hjust = 0.5,
                                             size = "14", color = "#646464"),
                   axis.title = element_text(face = "bold", size = "12", color = "#646464"),
                   axis.text = element_text(face = "bold", size = "9", color = "#646464"))
  return(list(ROC,figure))
}

