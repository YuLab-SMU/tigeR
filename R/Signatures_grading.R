#' @title Giving immunotherapy prognosis using MEMTS.
#' @description The function will return a vector calculated by Metastasis Related Epithelial-Mesenchymal Transition Signature.
#' @param exp_mtr An expression matrix which rownames are gene symbol and colnames are sample ID.
#' @export
#'

MEMTS_grading <- function(exp_mtr){
  MEMTS <- c('SNAI2','PFN2','NOTCH2','NID2','MEST','MATN2','LAMA1','ITGB3','GPX7','FBN2','ECM2','DPYSL3','BDNF')
  Expr <- dataPreprocess(exp_mtr,MEMTS,turn2HL = FALSE)
  result <- apply(Expr, 2, mean)
  names(result) <- colnames(Expr)
  return(result)
}
