#' @title Calculate Signature score with average mean.
#' @description Calculate Signature score with average mean.
#' @param exp_mtr An expression matrix which rownames are gene symbol and colnames are sample ID.
#' @param Signature Gene Signature.

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

  Expr <- Expr[names(Signature),]

  i <- 0
  Expr <- t(apply(Expr, 1, function(x){
    i <- i+1
    x*as.numeric(Signature[i])}))

  result <- apply(Expr, 2, sum)
  names(result) <- colnames(Expr)
  return(result)
}


#' @title Calculate Signature score with ZScore and PCA.
#' @description Calculate Signature score with ZScore and PCA.
#' @param exp_mtr An expression matrix which rownames are gene symbol and colnames are sample ID.
#' @param Signature Gene Signature.

ZScore_PCA_signature <- function(exp_mtr, Signature){
  Expr <- dataPreprocess(exp_mtr, Signature, turn2HL = FALSE)
  if(all(is.na(Expr)))
    return(rep(0,ncol(Expr)))

  average <- apply(Expr, 2, mean)
  standard_error <- apply(Expr, 2, stats::sd)

  ZScore <- t(apply(Expr,1,function(x){
    (x - average) / standard_error
  }))

  if(all(is.na(ZScore)))
    return(rep(0,ncol(Expr)))

  ZScore <-
    apply(ZScore, 2, function(x){
      if(all(is.na(x)))
        x <- rep(0, length(x))
      x
    })
  result  <- stats::prcomp(na.omit(ZScore), center = FALSE, scale = FALSE)$rotation[,1]
  return(result)
}
