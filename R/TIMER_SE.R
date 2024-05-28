#' @title TIMER deconvolution
#' @description use TIMER to predict TME
#' @param SE a SummarizedExperiment object contains the bulk RNA-seq dataset that you want to use for deconvolution and obtaining its cell fraction.
#' @param type the cancer type of data.
#' @export

TIMER_SE <- function(SE,type="SKCM"){
  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)
  TIMER.Immune <- readRDS(system.file("extdata", "TIMER.Immune.rds", package = "tigeR", mustWork = TRUE))

  co_genes <- intersect(rownames(exp_mtr),rownames(TIMER.Immune[[1]]))
  pre_rmBE <- cbind(exp_mtr[co_genes,],TIMER.Immune[[1]][co_genes,])
  batch <- as.factor(c(rep("Tumor",ncol(exp_mtr)),rep("Immune",ncol(TIMER.Immune[[1]]))))
  post_rmBE <- sva::ComBat(pre_rmBE, batch)

  tumor_exp <- post_rmBE[,seq_along(exp_mtr[1,])]
  immune_exp <- post_rmBE[,(ncol(exp_mtr) + 1):ncol(post_rmBE)]

  feature_matrix <- as.data.frame(
    lapply(TIMER.Immune[[2]], function(x){
      apply(immune_exp[,x],1,median)
    }))

  g <- unique(as.vector(apply(feature_matrix, 2, function(x){
    rownames(feature_matrix)[tail(order(x),ceiling(length(x)/100))]
  })))
  feature_matrix <- feature_matrix[!rownames(feature_matrix) %in% g,]

  TIMER.Markers <- readRDS(system.file("extdata", "TIMER.Markers.rds", package = "tigeR", mustWork = TRUE))

  selected_genes <- intersect(TIMER.Markers[[type]],rownames(feature_matrix))

  cancer.expression <- tumor_exp[selected_genes,]
  feature.expression <- feature_matrix[selected_genes,]

  fraction_matrix <-
    apply(cancer.expression, 2, function(x) {
      fraction <- stats::lsfit(feature.expression,x,intercept=FALSE)$coefficients

      drop <- c()
      while(any(fraction<0)){
        drop <- c(drop,which.min(fraction))
        fraction <- rep(0,length(fraction))
        fraction[-drop] <- stats::lsfit(feature.expression[,-drop],x,intercept=FALSE)$coefficients
      }
      fraction
    })
  rownames(fraction_matrix) <- colnames(feature_matrix)
  SummarizedExperiment::SummarizedExperiment(assays=S4Vectors::SimpleList(fraction_matrix),
                                             colData=S4Vectors::DataFrame(SummarizedExperiment::colData(SE)),
                                             checkDimnames=TRUE)
}
