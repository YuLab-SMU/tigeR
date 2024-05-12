#' @title Perform differential expression analysis and survival analysis.
#' @description perform differential expression analysis and survival analysis.
#' @param SE a SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain treatment information names Treatment.
#' @param geneSet The geneSet which you wanted.
#' @param method the method for calculating gene set scores which has several options: Average_mean, Weighted_mean, or GSVA. The method should be set to NULL if the length of the parameter geneSet is 1. This means that if you are working with only one gene, the specific calculation method may not be applicable or necessary.
#' @param PT_drop If TRUE, only Untreated patient will be use for model training.
#' @export

integrate_analysis <- function(SE, geneSet=NULL, method=NULL, PT_drop=TRUE){
  if(is.null(geneSet))
    geneSet <- rownames(assay(SE))

  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)
  meta <- bind_meta(SE, isList)
  Score <- Core(exp_mtr, geneSet, method)

  R_vs_NR <- DEA_Response(Score, meta, PT_drop)
  names(R_vs_NR) <- c('log2(FC)','P','Score')
  Pre_vs_Post <- DEA_Treatment(Score, meta)
  names(Pre_vs_Post) <- c('log2(FC)','P','Score')
  Survival <- survival_Score(Score, meta[,c(1,5,9,10)], PT_drop)
  names(Survival) <- c('HR','P','Score')

  result <- list(R_vs_NR, Pre_vs_Post, Survival)
  names(result) <- c('Response vs Non-Response',
                     'Pre-Therapy vs Post-Therapy',
                     'Survival')
  return(result)
}


#' @title Perform batch differential expression analysis and survival analysis.
#' @description Perform batch differential expression analysis and survival analysis in certain gene and return the result.
#' @param SE a SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain treatment information names Treatment.
#' @param geneSet The genes you want to use for anaylsis.
#' @param method the method for calculating gene set scores which has several options: Average_mean, Weighted_mean, or GSVA. The method should be set to NULL if the length of the parameter geneSet is 1. This means that if you are working with only one gene, the specific calculation method may not be applicable or necessary.
#' @param PT_drop If TRUE, only Untreated patient will be use for model training.
#' @export

Immunotherapy_Response_Batch <- function(SE, geneSet=NULL, method=NULL, PT_drop=TRUE){
  lapply(geneSet,
         function(x) integrate_analysis(SE,geneSet=x,method, PT_drop))
}

#' @title Calculating differential expression score between responder and non_responder.
#' @description Differential expression score was alculated using the following formula: −𝑆𝐼𝐺𝑁(𝑙𝑜𝑔2(𝐹𝐶)) × 𝑙𝑜𝑔10(𝑃), where 𝐹𝐶 represents the fold change and 𝑃 represents the P value derived from the Wilcoxon rank-sum test
#' @param Score an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information names response.
#' @param meta The geneSet which you wanted.
#' @param PT_drop If TRUE, only Untreated patient will be use for model training.

DEA_Response <- function(Score, meta, PT_drop){
  idx_R <- which(meta$response_NR == 'R')
  idx_N <- which(meta$response_NR == 'N')

  if(PT_drop){
    idx_UT <- which(meta$Treatment == 'PRE')
    idx_R <- intersect(idx_R,idx_UT)
    idx_N <- intersect(idx_N,idx_UT)
  }

  if(length(idx_R)==0||length(idx_N)==0){
    warning("The data set must have both Responder and Non-Responder.")
    return(data.frame(log2FC=NA,P_value=NA,DEA_Score=NA))
  }

  FC <- abs(mean(Score[idx_R])/mean(Score[idx_N]))
  log2FC <- log2(FC)

  if(length(idx_R)<2||length(idx_N)<2){
    P <- NA
    message("The p value of differential analysis between groups N and R is not availiable. (Both N and R groups must contain more than 2 samples)")
  }else{
    P <- stats::t.test(Score[idx_R],Score[idx_N])$p.value
  }

  Score <- -sign(log2FC) * log10(P)

  result <- data.frame(log2FC,P_value=P,DEA_Score=Score)
  return(result)
}


#' @title Calculating differential expression score between treated and untreated patients.
#' @description Differential expression score was alculated using the following formula: −𝑆𝐼𝐺𝑁(𝑙𝑜𝑔2(𝐹𝐶)) × 𝑙𝑜𝑔10(𝑃), where 𝐹𝐶 represents the fold change and 𝑃 represents the P value derived from the Wilcoxon rank-sum test
#' @param Score an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information names response.
#' @param meta The geneSet which you wanted.

DEA_Treatment <- function(Score, meta){
  idx_Pre <- which(meta$Treatment == 'PRE')
  idx_Post <- which(meta$Treatment %in% c('POST','ON'))
  if(length(idx_Pre)==0||length(idx_Post)==0){
    warning("The data set must have both Pre-Treatment and Post-Treatment samples.")
    return(data.frame(log2FC=NA,P_value=NA,DEA_Score=NA))
  }

  FC <- abs(mean(Score[idx_Pre])/mean(Score[idx_Post]))
  log2FC <- log2(FC)

  if(length(idx_Pre)<2||length(idx_Post)<2){
    P <- NA
    message("The p value of differential analysis between groups Pre and Post is not availiable. (Both Pre and Post groups must contain more than 2 samples)")
  }else{
    P <- stats::t.test(Score[idx_Pre],Score[idx_Post])$p.value
  }

  Score <- -sign(log2FC) * log10(P)

  result <- data.frame(log2FC,P_value=P,DEA_Score=Score)
  return(result)
}

#' @title Calculating survival score of Untreated patients.
#' @description Survival score was calculated using the following formula: −𝑆𝐼𝐺𝑁(𝑙𝑜𝑔2(𝐻𝑅)) × 𝑙𝑜𝑔10(𝑃), where 𝐻𝑅 represents the hazard ratio and 𝑃 represents the P value derived from univariate Cox regression analysis.
#' @param Score an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information names response.
#' @param meta The geneSet which you wanted.
#' @param PT_drop If TRUE, only Untreated patient will be use for model training.

survival_Score <- function(Score, meta, PT_drop){
  if(PT_drop){
    idx_UT <- which(meta$Treatment == "PRE")
  }else{
    idx_UT <- seq_along(meta$Treatment)
  }

  result <- matrix_cox(ifelse(Score>=stats::median(Score),1,0)[idx_UT],
                       meta[idx_UT, -2])
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
  df$status <- as.numeric(sub('Alive','0',sub('Dead','1',df$status)))

  cox <- summary(coxph(Surv(time, status) ~ exp, data = df))
  HR <- signif(cox$coefficients[2])
  P <- cox$coefficients[5]
  Score <- -sign(log2(HR)) * log10(P)
  return(c(HR,P,Score))
}
