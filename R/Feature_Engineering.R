#' @title Merge SummarizedExperiment object
#' @description Merge SummarizedExperiment object
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @export

merge_SE <- function(SE){
  browser()
  isList <- is.list(SE)
  expr_mtr <- bind_mtr(SE,isList)
  col_data <- bind_meta(SE,isList)

  SE_obj <- SummarizedExperiment::SummarizedExperiment(assays=S4Vectors::SimpleList(expr_mtr),
                                                       colData=S4Vectors::DataFrame(col_data),
                                                       checkDimnames=TRUE)
}
