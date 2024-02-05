#' @title plot differential result(Responder vs Non-Responder or Pre-Treatment vs Post-Treatment)
#' @description plot differential result(Responder vs Non-Responder or Pre-Treatment vs Post-Treatment)
#' @param SE SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param gene is the Gene or Gene set you are interested in.
#' @param type 'Treatment' or 'Response'.the type of analysis you want to perform(Responder vs Non-Responder or Pre-Treatment vs Post-Treatment)
#' @param method the method for calculating gene set scores. Can be NULL if the length of parameter gene is 1.
#' @export

plt_diff <- function(SE, gene, type, method='Average_mean'){
  type <- match.arg(type, c('Response','Treatment'))
  method <- match.arg(method, c('Average_mean','GSVA','Weighted_mean'))

  if(type == 'Response'){
    df <- plt_Preprocess(gene, SE, method, 'R vs NR')
    plt <- plt_style(df) + ggplot2::ggtitle("Responder vs Non-Responder")
  }else{
    df <- plt_Preprocess(gene, SE, method, 'T vs UT')
    plt <- plt_style(df) + ggplot2::ggtitle("Treatment vs UnTreatment")
  }
  return(plt +
           ggplot2::labs(title=NULL,x=NULL,y=paste0('log2(', method,' + 1)')) +
           ggplot2::theme(plot.title = element_text(hjust = 0.5)))
}


#' @title plot differential result(Pre-Treatment vs Post-Treatment)
#' @description The association between gene expression and overall survival in the immunotherapy data was calculated using univariate Cox regression analysis.
#' @param SE SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param gene is the Gene or Gene set you are interested in.
#' @param method the method for calculating gene set scores. Can be NULL if the length of parameter gene is 1.
#' @importFrom SummarizedExperiment assay
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom survival survfit
#' @importFrom survival Surv
#' @importFrom survminer ggsurvplot
#' @importFrom ggplot2 theme_bw
#' @importFrom stats median
#' @importFrom stats na.omit
#' @export

plt_surv <- function(SE, gene, method='Average_mean'){
  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)
  meta <- bind_meta(SE, isList)

  Sc <- Core(exp_mtr, gene, method)
  Score <- as.numeric(ifelse(Sc>=median(Sc),1,0))

  time <- as.numeric(meta$overall.survival..days.)
  status <- sub('Dead','1', meta$vital.status) %>% sub('Alive','0',.)

  df <- data.frame(time,status,Score) %>% na.omit() %>% lapply(as.numeric) %>% as.data.frame()

  fit <- survfit(Surv(time, status) ~ Score, data = df)

  P <-
  ggsurvplot(fit,
             data = df,
             pval = TRUE,
             conf.int = TRUE,
             legend.title = ifelse(length(gene)==1,gene,"Score"),
             legend.labs = c('Low','HIGH'),
             risk.table = TRUE,
             risk.table.title = 'Number at risk',
             risk.table.col = "strata",
             linetype = "strata",
             surv.median.line = "hv",
             ggtheme = theme_bw())
  P$plot <- P$plot +
    ggtitle("Survival analysis") +
    theme(plot.title = element_text(hjust = 0.5))
  if(length(gene) > 1)
    P$table$labels$y <- method

  return(P)
}
