#' @title Calculating Signature score of existing immunotherapy response Signature.
#' @description By employing the Signature_calculation() function, you can obtain a comprehensive signature score matrix for the 23 signatures in tigeR. In this matrix, the columns represent the signature scores, and the rows denote the sample names.
#' @param SE a SummarizedExperiment object for which you want to calculate the Signature Score.
#' @param exp_mtr an expression matrix for which you want to calculate the Signature Score.
#' @param meta meta data of samples
#' @param threshold the threshold for signature to discretize patients into R and N group
#' @param positive description
#' @param Signature a genes vector represents user-defined signature for Immunotherapy response. If NULL, the function will only calculate 23 built-in signatures in tigeR.
#' @param method the method for calculating gene set scores which has several options: Average_mean, Weighted_mean, or GSVA. The method can be set to NULL if the length of the parameter geneSet is 1. This means that if you are working with only one gene, the specific calculation method may not be applicable or necessary.
#' @param PT_drop If TRUE, only Untreated patient will be use for model training.
#' @export

response_prediction <- function(SE=NULL, exp_mtr=NULL, meta=NULL, threshold,positive="R",Signature=NULL, method="Average_mean",PT_drop=TRUE){
  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)
  meta <- bind_meta(SE, isList)

  sig <- Core(exp_mtr, Signature, method)
  negative <- ifelse(positive=="R","N","R")
  group <- ifelse(sig >= threshold,positive,negative)
  df <- data.frame(sample=names(sig),value=sig,group=group)
  df <- dplyr::arrange(df,dplyr::desc("value"))
  df$sample <- factor(df$sample,levels = df$sample)

  ggplot(data = df, aes(x = .data$sample, y = .data$value, fill = .data$group)) +
    geom_bar(stat = "identity")  +
    theme(plot.title = element_text(hjust = 0,size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"),
          axis.title.x = element_text(size = 12,face="bold"),
          axis.title.y = element_text(size = 12,face="bold"),
          axis.text.x = element_text(size = 8,face="bold",angle = 90),
          axis.line = element_blank())
}

#' @title Calculating Signature score of existing immunotherapy response Signature.
#' @description By employing the Signature_calculation() function, you can obtain a comprehensive signature score matrix for the 23 signatures in tigeR. In this matrix, the columns represent the signature scores, and the rows denote the sample names.
#' @param SE a SummarizedExperiment object for which you want to calculate the Signature Score.
#' @param exp_mtr an expression matrix for which you want to calculate the Signature Score.
#' @param meta meta data of samples
#' @param threshold the threshold for signature to discretize patients into R and N group
#' @param positive description
#' @param Signature a genes vector represents user-defined signature for Immunotherapy response. If NULL, the function will only calculate 23 built-in signatures in tigeR.
#' @param method the method for calculating gene set scores which has several options: Average_mean, Weighted_mean, or GSVA. The method can be set to NULL if the length of the parameter geneSet is 1. This means that if you are working with only one gene, the specific calculation method may not be applicable or necessary.
#' @param PT_drop If TRUE, only Untreated patient will be use for model training.
#' @param show.val If TRUE, the value will be show in the heatplot.
#' @param ZS dd
#' @import patchwork
#' @export

pred_response <- function(SE=NULL, exp_mtr=NULL, meta=NULL, threshold=0.8,
                     positive="R", Signature=NULL, method="Average_mean",
                     PT_drop=TRUE,show.val=TRUE, sort_by="rankscore",
                     group_by="rankscore", show.real=TRUE, text_col="black",
                     rankscore=TRUE,ZS=TRUE){
  if(!missing(SE)){
    isList <- is.list(SE)
    exp_mtr <- bind_mtr(SE, isList)
    meta <- bind_meta(SE, isList)
  }

  if(PT_drop){
    idx_UT <- which(meta$Treatment == 'PRE')
    if(length(idx_UT) == 0)
      stop("Only Untreated patients can be use to perform survival analysis!")
    exp_mtr <- exp_mtr[,idx_UT,drop=FALSE]
    meta <- meta[idx_UT,,drop=FALSE]
  }

  Average_mean_Sigs <- NULL
  Weighted_mean_Sigs <- NULL
  ZScore_PCA_Sigs <- NULL
  data(Average_mean_Sigs, package = 'tigeR', envir = current_env())
  data(Weighted_mean_Sigs, package = 'tigeR', envir = current_env())
  data(ZScore_PCA_Sigs, package = 'tigeR', envir = current_env())
  df <- data.frame(IRS=IRS_grading(exp_mtr),tGE8=tGE8_grading(exp_mtr))
  for (i in Average_mean_Sigs) {
    df <- cbind(df, average_mean_signature(exp_mtr, i))
  }
  for (i in Weighted_mean_Sigs) {
    df <- cbind(df,weight_mean_signature(exp_mtr,i))
  }
  for (i in ZScore_PCA_Sigs) {
    df <- cbind(df,ZScore_PCA_signature(exp_mtr,i))
  }
  colnames(df) <- c('IRS','tGE8',
                    names(Average_mean_Sigs),
                    names(Weighted_mean_Sigs),
                    names(ZScore_PCA_Sigs))
  if(ZS){
    if(!is.null(Signature)){
      sig <- Core(exp_mtr, Signature, method)
      df <- cbind(Customed.Signature=sig,df)
    }
  }

  negative <- ifelse(positive=="R","N","R")

  vl <- as.matrix(
  apply(df,2,function(x){
    (x - mean(x))/stats::sd(x)
  }))
  if(!ZS){
    sig <- Core(exp_mtr, Signature, method)
    vl <- cbind(Customed.Signature=sig,vl)
  }
  vl <- t(na.omit(t(vl)))
  rank_mtr <-
    apply(vl, 2, function(x){
      order(x)
    })
  rownames(rank_mtr) <- rownames(vl)
  rank_list <-
  lapply(seq_along(rank_mtr[1,]),function(i){
    rank_mtr[,i]
  })

  names(rank_list) <- colnames(rank_mtr)
  result <- RobustRankAggreg::aggregateRanks(rank_list,exact = TRUE)
  result$Name <- as.numeric(result$Name)
  result <- dplyr::arrange(result,"Name")
  Score <- (result$Score - mean(result$Score))/stats::sd(result$Score)

  right.mtr <- data.frame(Sample=rownames(vl),rankscore=Score,vl) %>%
    dplyr::arrange(!!sym(sort_by))
  if(!rankscore)
    right.mtr <- right.mtr[,-2]
  right.mtr$Sample <- factor(right.mtr$Sample,levels = right.mtr$Sample)
  group <- ifelse(right.mtr[,group_by]>=threshold,positive,negative)
  lv <- reshape2::melt(right.mtr,id.vars="Sample")
  colnames(lv) <- c("Sample", "Signature","Score")
  plt.r <-
    ggplot(lv, aes(x = .data$Signature, y = .data$Sample, fill = .data$Score)) +
    geom_tile()  +
    scale_fill_gradient(low = "white", high = "red") +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size=10,color = text_col, angle = 45, vjust = 1, hjust=1),
          axis.text.y = element_text(size=10,color = text_col, margin = margin(r=10)))
  if(show.val)
    plt.r <- plt.r + geom_text(aes(label = sprintf("%.2f", Score)), size = 2)

  df.l <- data.frame(Sample=factor(rownames(right.mtr),levels = right.mtr$Sample),
                     prediction=group)
  if(show.real)
    df.l$real <- meta[rownames(right.mtr),]$response_NR
  df.l <- reshape2::melt(df.l,id.var="Sample")
  colnames(df.l) <- c("Sample","Group","Value")
  plt.l <-
  ggplot(df.l, aes(x = .data$Group, y = .data$Sample, fill = .data$Value)) +
    geom_tile() +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "left",
          legend.title = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=10,color=text_col,angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_blank())
  return(plt.l + plt.r + plot_layout(widths = c(1,30)))
}
