#' @title Calculating differential expression score between responder and non_responder.
#' @description Differential expression score was alculated using the following formula: −𝑆𝐼𝐺𝑁(𝑙𝑜𝑔2(𝐹𝐶)) × 𝑙𝑜𝑔10(𝑃), where 𝐹𝐶 represents the fold change and 𝑃 represents the P value derived from the Wilcoxon rank-sum test
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information names response.
#' @param gene The gene which you wanted.
#' @importFrom SummarizedExperiment assay

DEA_Response <- function(SE, gene){
  browser()
  exp_mtr <- bind_mtr(SE, 0)[gene,]
  meta <- bind_meta(SE, 0)

  idx_R <- which(meta$response_NR == 'R')
  idx_N <- which(meta$response_NR == 'N')

  #log2FC <- log2(apply(exp_mtr[,idx_R], 1, mean)/apply(exp_mtr[,idx_N], 1, mean))
  #P <- apply(exp_mtr, 1, matrix_t.test, P=idx_R, N=idx_N)
  log2FC <- log2(mean(exp_mtr[idx_R])/mean(exp_mtr[idx_N]))
  P <- matrix_t.test(exp_mtr,P=idx_R,N=idx_N)
  Score <- -sign(log2FC) * log10(P)

  result <- data.frame(log2FC,P_value=P,DEA_Score=Score)
  return(result)
}


#' @title Calculating differential expression score between treated and untreated patients.
#' @description Differential expression score was alculated using the following formula: −𝑆𝐼𝐺𝑁(𝑙𝑜𝑔2(𝐹𝐶)) × 𝑙𝑜𝑔10(𝑃), where 𝐹𝐶 represents the fold change and 𝑃 represents the P value derived from the Wilcoxon rank-sum test
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain treatment information names Treatment.
#' @param gene The gene which you wanted.
#' @importFrom SummarizedExperiment assay
#'

DEA_Treatment <- function(SE, gene){
  exp_mtr <- bind_mtr(SE, 0)[gene,]
  meta <- bind_meta(SE, 0)

  idx_Pre <- which(meta$Treatment == 'PRE')
  idx_Post <- which(meta$Treatment %in% c('POST','ON'))

  #log2FC <- log2(apply(exp_mtr[,idx_Pre], 1, mean)/apply(exp_mtr[,idx_Post], 1, mean))
  #P <- apply(exp_mtr, 1, matrix_t.test, P=idx_Pre, N=idx_Post)
  log2FC <- log2(mean(exp_mtr[idx_Pre])/mean(exp_mtr[idx_Post]))
  P <- matrix_t.test(exp_mtr,P=idx_Pre,N=idx_Post)
  Score <- -sign(log2FC) * log10(P)

  result <- data.frame(log2FC,P_value=P,DEA_Score=Score)
  return(result)
}

#' @title Calculating survival score of patients.
#' @description Survival score was calculated using the following formula: −𝑆𝐼𝐺𝑁(𝑙𝑜𝑔2(𝐻𝑅)) × 𝑙𝑜𝑔10(𝑃), where 𝐻𝑅 represents the hazard ratio and 𝑃 represents the P value derived from univariate Cox regression analysis.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain event time(names time), event(names status).
#' @param gene The gene which you wanted.
#' @importFrom magrittr %<>%
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment assay
#' @importFrom stats median

survival_Score <- function(SE, gene){
  meta <- as.data.frame(SE@colData)[,c(1,9,10)]
  mtr <- assay(SE)[rownames(assay(SE)) == gene,]
  if(length(gene) == 1){
    M <- median(mtr)
    for (i in 1:length(mtr)) {
      mtr[i] <- ifelse(mtr[i] >= M,1,0)
    }
    result <- matrix_cox(mtr,meta)
    names(result) <- c('HR', 'P', 'Score')
  }
  else{
    M <- apply(mtr, 1, median)
    for (i in 1:nrow(mtr)) {
      mtr[i,] <- ifelse(mtr[i,] >= M[i],1,0)
    }
    result <- t(apply(mtr, 1, matrix_cox,meta=meta))
    colnames(result) <- c('HR', 'P', 'Score')
  }

  return(result)
}


#' @title Perform differential expression analysis and survival analysis.
#' @description Perform differential expression analysis and survival analysis in certain gene and return the result.
#' @param SE a SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain treatment information names Treatment.
#' @param gene The gene which you wanted.
#' @export

Immunotherapy_Response <- function(SE, gene=NULL){
  if(is.null(gene))
    gene <- rownames(assay(SE))

  R_vs_NR <- DEA_Response(SE, gene)
  names(R_vs_NR) <- c('log2(FC)','P','Score')
  Pre_vs_Post <- DEA_Treatment(SE, gene)
  names(Pre_vs_Post) <- c('log2(FC)','P','Score')
  Survival <- survival_Score(SE, gene)
  names(Survival) <- c('HR','P','Score')

  result <- list(R_vs_NR, Pre_vs_Post, Survival)
  names(result) <- c('Response vs Non-Response',
                     'Pre-Therapy vs Post-Therapy',
                     'Survival')
  return(result)
}
