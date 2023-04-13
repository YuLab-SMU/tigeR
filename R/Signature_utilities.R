#' @title Calculate Signature score with average mean.
#' @description Calculate Signature score with average mean.
#' @param exp_mtr An expression matrix which rownames are gene symbol and colnames are sample ID.
#' @param Signature Gene Signature.
#' @export
#'

average_mean_signature <- function(exp_mtr, Signature){
  Expr <- dataPreprocess(exp_mtr, Signature, turn2HL = FALSE)
  rownames_identify <- rownames(Expr)
  for (gene in Signature) {
    if(!(gene %in% rownames(Expr))){
      vec <- rep(0, ncol(Expr))
      Expr <- rbind(Expr, vec)
    }
  }


  idx <- which(!Signature %in% rownames(Expr))
  if(length(idx) != 0)
    rownames(Expr) <- c(rownames_identify, Signature[idx])

  result <- apply(Expr, 2, mean)
  names(result) <- colnames(Expr)
  return(result)
}


#' @title Calculate Signature score with weighted mean.
#' @description Calculate Signature score with weighted mean.
#' @param exp_mtr An expression matrix which rownames are gene symbol and colnames are sample ID.
#' @param Signature Gene Signature.
#' @export
#'

weight_mean_signature <- function(exp_mtr, Signature){
  if(length(names(Signature)) != length(Signature)){
    warning('The length of Signature must equal to the length of weight!')
  }

  Expr <- dataPreprocess(exp_mtr, names(Signature), turn2HL = FALSE)
  rownames_identify <- rownames(Expr)
  for (gene in names(Signature)) {
    if(!(gene %in% rownames(Expr))){
      vec <- rep(0, ncol(Expr))
      Expr <- rbind(Expr, vec)
    }
  }


  idx <- which(!names(Signature) %in% rownames(Expr))
  if(length(idx) != 0)
    rownames(Expr) <- c(rownames_identify, names(Signature[idx]))


  for (gene in names(Signature)) {
    Expr[rownames(Expr) == gene,] <- Expr[rownames(Expr) == gene,] * Signature[names(Signature) == gene]
  }

  result <- apply(Expr, 2, sum)
  names(result) <- colnames(Expr)
  return(result)
}
