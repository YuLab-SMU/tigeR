#' @title Perform differential expression analysis and survival analysis.
#' @description Perform differential expression analysis and survival analysis in certain gene and return the result.
#' @param SE a SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain treatment information names Treatment.
#' @param geneSet The geneSet which you wanted.
#' @param method the method for calculating gene set scores. Can be NULL if the length of parameter gene is 1.
#' @export

Immunotherapy_Response <- function(SE, geneSet=NULL, method){
  if(is.null(geneSet))
    geneSet <- rownames(assay(SE))

  R_vs_NR <- DEA_Response(SE, geneSet, method)
  names(R_vs_NR) <- c('log2(FC)','P','Score')
  Pre_vs_Post <- DEA_Treatment(SE, geneSet, method)
  names(Pre_vs_Post) <- c('log2(FC)','P','Score')
  Survival <- survival_Score(SE, geneSet, method)
  names(Survival) <- c('HR','P','Score')

  result <- list(R_vs_NR, Pre_vs_Post, Survival)
  names(result) <- c('Response vs Non-Response',
                     'Pre-Therapy vs Post-Therapy',
                     'Survival')
  return(result)
}


#' @title Calculating differential expression score between responder and non_responder.
#' @description Differential expression score was alculated using the following formula: −𝑆𝐼𝐺𝑁(𝑙𝑜𝑔2(𝐹𝐶)) × 𝑙𝑜𝑔10(𝑃), where 𝐹𝐶 represents the fold change and 𝑃 represents the P value derived from the Wilcoxon rank-sum test
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information names response.
#' @param geneSet The geneSet which you wanted.
#' @param method the method for calculating gene set scores. Can be NULL if the length of parameter gene is 1.
#' @importFrom SummarizedExperiment assay
#' @importFrom stats t.test

DEA_Response <- function(SE, geneSet, method){
  browser()
  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)
  meta <- bind_meta(SE, isList)

  idx_R <- which(meta$response_NR == 'R')
  idx_N <- which(meta$response_NR == 'N')

  Score <- Core(exp_mtr, geneSet, method)
  FC <- abs(mean(Score[idx_R])/mean(Score[idx_N]))
  log2FC <- log2(FC)
  P <- t.test(Score[idx_R],Score[idx_N])$p.value
  Score <- -sign(log2FC) * log10(P)

  result <- data.frame(log2FC,P_value=P,DEA_Score=Score)
  return(result)
}


#' @title Calculating differential expression score between treated and untreated patients.
#' @description Differential expression score was alculated using the following formula: −𝑆𝐼𝐺𝑁(𝑙𝑜𝑔2(𝐹𝐶)) × 𝑙𝑜𝑔10(𝑃), where 𝐹𝐶 represents the fold change and 𝑃 represents the P value derived from the Wilcoxon rank-sum test
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain treatment information names Treatment.
#' @param geneSet The geneSet which you wanted.
#' @param method the method for calculating gene set scores. Can be NULL if the length of parameter gene is 1.
#' @importFrom SummarizedExperiment assay
#' @importFrom stats t.test

DEA_Treatment <- function(SE, geneSet, method){
  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)
  meta <- bind_meta(SE, isList)

  idx_Pre <- which(meta$Treatment == 'PRE')
  idx_Post <- which(meta$Treatment %in% c('POST','ON'))

  Score <- Core(exp_mtr, geneSet, method)
  FC <- abs(mean(Score[idx_Pre])/mean(Score[idx_Post]))
  log2FC <- log2(FC)
  P <- t.test(Score[idx_Pre],Score[idx_Post])$p.value
  Score <- -sign(log2FC) * log10(P)

  result <- data.frame(log2FC,P_value=P,DEA_Score=Score)
  return(result)
}

#' @title Calculating survival score of patients.
#' @description Survival score was calculated using the following formula: −𝑆𝐼𝐺𝑁(𝑙𝑜𝑔2(𝐻𝑅)) × 𝑙𝑜𝑔10(𝑃), where 𝐻𝑅 represents the hazard ratio and 𝑃 represents the P value derived from univariate Cox regression analysis.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain event time(names time), event(names status).
#' @param geneSet The geneSet which you wanted.
#' @param method the method for calculating gene set scores. Can be NULL if the length of parameter gene is 1.
#' @importFrom magrittr %<>%
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment assay
#' @importFrom stats median

survival_Score <- function(SE, geneSet, method){
  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)
  meta <- bind_meta(SE, isList)[,c(1,9,10)]

  if(length(geneSet) > 1){
    Sc <- Core(exp_mtr, geneSet, method)
  }else{
    Sc <- exp_mtr[geneSet,]
  }

  M <- median(Sc)

  for (i in 1:length(exp_mtr)) {
    Score <- ifelse(Sc >= M,1,0)
  }
  result <- matrix_cox(Score, meta)
  names(result) <- c('HR', 'P', 'Score')

  return(result)
}


#' @title cox regression
#' @description wait to write
#' @param V the expression vector of a gene
#' @param meta the meta information.
#' @import survival
#' @importFrom stats na.omit

matrix_cox <- function(V,meta){
  df <- na.omit(cbind(meta,V))
  colnames(df) <- c('ID','time','status','exp')
  df$status %<>% {sub('Dead','1',.)} %>% {sub('Alive','0',.)} %>% as.numeric()

  cox <- summary(coxph(Surv(time, status) ~ exp, data = df))
  HR <- signif(cox$coefficients[2])
  P <- cox$coefficients[5]
  Score <- -sign(log2(HR)) * log10(P)
  return(c(HR,P,Score))
}


#' @title count geneset score by different method
#' @description wait to write
#' @param exp_mtr an expression matrix.
#' @param geneSet The geneSet which you wanted.
#' @param method the method for calculating gene set scores. Can be NULL if the length of parameter gene is 1.
#' @import survival

Core <- function(exp_mtr, geneSet, method){
  if(length(geneSet)==1)
    return(exp_mtr[geneSet,])

  exp_mtr <- stats::na.omit(t(apply(exp_mtr, 1, function(x){
    if(all(x==0)){
      return(rep(NA,length(x)))
    }else{
      return(x)
    }
  })))

  Score <-
  switch(method,
         Average_mean = apply(exp_mtr[geneSet,], 2, mean),
         GSVA = GSVA::gsvaParam(exp_mtr, geneSets=list(geneSet)) %>% GSVA::gsva(),
         Weighted_mean = weight_mean_signature(exp_mtr, geneSet))

  return(as.vector(Score))
}
