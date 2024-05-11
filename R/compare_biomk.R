#' @title Calculating Signature score of existing immunotherapy response Signature.
#' @description By employing
#' @param SE a SummarizedExperiment object for which you want to calculate the Signature Score.
#' @param Signature a genes vector represents user-defined signature for Immunotherapy response. If NULL, the function will only calculate 23 built-in signatures in tigeR.
#' @param method the method for calculating gene set scores which has several options: Average_mean, Weighted_mean, or GSVA. The method can be set to NULL if the length of the parameter geneSet is 1. This means that if you are working with only one gene, the specific calculation method may not be applicable or necessary.
#' @param PT_drop If TRUE, only Untreated patient will be use for model training.
#' @param show.val If TRUE, the value will be show in the heatplot.
#' @param val.size an integer represents the size of AUC value.
#' @importFrom SummarizedExperiment colData
#' @export

compare_biomk <- function(SE=NULL, Signature=NULL, method="Average_mean",PT_drop=TRUE,show.val=TRUE,val.size=2){
  if(inherits(SE,"SummarizedExperiment"))
    SE <- list(SE)

  auc.l <-
    lapply(SE, function(x){
      Sc <- score_biomk(SE=x,Signature = Signature,method = method, PT_drop = PT_drop)
      auc.v <-
        apply(Sc, 2, function(y){
          response <- colData(x)[rownames(Sc),]$response_NR
          idx_UNK <- which(response=="UNK")
          if(length(idx_UNK)!=0){
            response <- response[-idx_UNK]
            y <- y[-idx_UNK]
          }
          if(all(is.na(y)))
            return(0)

          pROC::roc(response, y)$auc
        })
    })

  cnames <- unlist(
    lapply(SE, function(x){
      unique(x$dataset_group)
    }))
  auc.df <- do.call(cbind, auc.l)
  colnames(auc.df) <- cnames

  R.df <- reshape::melt(auc.df)
  colnames(R.df) <- c("Signature", "Data_set", "AUC")

  R.plt <-
    ggplot(R.df, aes(y = .data$Signature, x = .data$Data_set, fill = .data$AUC)) +
    geom_tile() +
    scale_fill_gradient2(low="blue",mid="white", high="red",midpoint = 0.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(color = "black",angle = 90,
                                     vjust = 0.5, hjust=1),
          axis.text.y = element_text(color = "black"))
  if(show.val)
    plt <- R.plt + geom_text(aes(label = sprintf("%.2f", .data$AUC)), size = val.size)
  return(plt)
}
