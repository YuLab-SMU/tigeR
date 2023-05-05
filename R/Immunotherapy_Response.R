#' @title Calculating differential expression score between responder and non_responder.
#' @description Differential expression score was alculated using the following formula: −𝑆𝐼𝐺𝑁(𝑙𝑜𝑔2(𝐹𝐶)) × 𝑙𝑜𝑔10(𝑃), where 𝐹𝐶 represents the fold change and 𝑃 represents the P value derived from the Wilcoxon rank-sum test
#' @param gene The gene which you wanted.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information names response.
#' @importFrom SummarizedExperiment assay

DEA_Response <- function(gene='CD274',SE){
  meta <- SE@colData
  exp <- assay(SE)[rownames(assay(SE)) == gene,]

  responder <- meta[meta$response_NR == 'R',]$sample_id
  non_responder <- meta[meta$response_NR == 'N',]$sample_id

  R_exp <- exp[names(exp) %in% responder]
  NR_exp <- exp[names(exp) %in% non_responder]

  if(length(R_exp) == 0 ||length(NR_exp) == 0)
    stop('There must be two different groups in you data!')

  FC <- mean(R_exp)/mean(NR_exp)
  P <- wilcox.test(R_exp,NR_exp)$p.value
  Score <- -sign(log2(FC)) * log10(P)

  result <- c(FC, P, Score)
  return(result)
}


#' @title Calculating differential expression score between treated and untreated patients.
#' @description Differential expression score was alculated using the following formula: −𝑆𝐼𝐺𝑁(𝑙𝑜𝑔2(𝐹𝐶)) × 𝑙𝑜𝑔10(𝑃), where 𝐹𝐶 represents the fold change and 𝑃 represents the P value derived from the Wilcoxon rank-sum test
#' @param gene The gene which you wanted.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain treatment information names Treatment.
#' @importFrom SummarizedExperiment assay
#'

DEA_Treatment <- function(gene='CD274',SE){
  meta <- SE@colData
  exp <- assay(SE)[rownames(assay(SE)) == gene,]

  Pre <- meta[meta$Treatment == 'PRE',]$sample_id
  Post <- meta[meta$Treatment %in% c('POST','ON'),]$sample_id

  Pre_exp <- exp[names(exp) %in% Pre]
  Post_exp <- exp[names(exp) %in% Post]

  if(length(Pre_exp) == 0 || length(Post_exp) == 0)
    stop('There must be two different groups in you data!')

  FC <- mean(Post_exp)/mean(Pre_exp)
  P <- wilcox.test(Pre_exp,Post_exp)$p.value
  Score <- -sign(log2(FC)) * log10(P)

  result <- c(FC, P, Score)
  return(result)
}

#' @title Calculating survival score of patients.
#' @description Survival score was calculated using the following formula: −𝑆𝐼𝐺𝑁(𝑙𝑜𝑔2(𝐻𝑅)) × 𝑙𝑜𝑔10(𝑃), where 𝐻𝑅 represents the hazard ratio and 𝑃 represents the P value derived from univariate Cox regression analysis.
#' @param gene The gene which you wanted.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain event time(names time), event(names status).
#' @import survival
#' @importFrom magrittr %<>%
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment assay
#'

survival_Score <- function(gene='CD274',SE){
  meta <- as.data.frame(SE@colData)[,c(1,9,10)]
  mtr <- assay(SE)[rownames(assay(SE)) == gene,]
  mtr %<>% {ifelse(.>=median(mtr),1,0)}

  df <- na.omit(cbind(meta,mtr))
  colnames(df) <- c('ID','time','status','exp')
  df$status %<>% {sub('Dead','1',.)} %>% {sub('Alive','0',.)} %>% as.numeric()

  cox <- summary(coxph(Surv(time, status) ~ exp, data = df))
  HR <- signif(cox$coefficients[2])
  P <- cox$coefficients[5]
  Score <- -sign(log2(HR)) * log10(P)

  result <- c(HR, P, Score)
  return(result)
}


#' @title Perform differential expression analysis and survival analysis.
#' @description Perform differential expression analysis and survival analysis in certain gene and return the result.
#' @param gene The gene which you wanted.
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain treatment information names Treatment.
#' @export

Immunotherapy_Response <- function(gene, SE){
  R_vs_NR <- DEA_Response(gene, SE)
  names(R_vs_NR) <- c('FC','P','Score')
  Pre_vs_Post <- DEA_Treatment(gene, SE)
  names(Pre_vs_Post) <- c('FC','P','Score')
  Survival <- survival_Score(gene, SE)
  names(Survival) <- c('HR','P','Score')

  result <- list(R_vs_NR, Pre_vs_Post, Survival)
  names(result) <- c('Response vs Non-Response',
                     'Pre-Therapy vs Post-Therapy',
                     'Survival')
  return(result)
}