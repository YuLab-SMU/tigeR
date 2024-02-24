#' @title plot differential result(Responder vs Non-Responder or Pre-Treatment vs Post-Treatment)
#' @description plot differential result(Responder vs Non-Responder or Pre-Treatment vs Post-Treatment)
#' @param SE SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param gene is the Gene or Gene set you are interested in.
#' @param type 'Treatment' or 'Response'.the type of analysis you want to perform(Responder vs Non-Responder or Pre-Treatment vs Post-Treatment)
#' @param method the method for calculating gene set scores which has several options: Average_mean, Weighted_mean, or GSVA. The method can be set to NULL if the length of the parameter geneSet is 1. This means that if you are working with only one gene, the specific calculation method may not be applicable or necessary.
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
#' @param method the method for calculating gene set scores which has several options: Average_mean, Weighted_mean, or GSVA. The method can be set to NULL if the length of the parameter geneSet is 1. This means that if you are working with only one gene, the specific calculation method may not be applicable or necessary.
#' @param style the ploting style. c("raw,"elegant","brief")
#' @param conf.int logical value. If TRUE, plots confidence interval.
#' @importFrom SummarizedExperiment assay
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom survival survfit
#' @importFrom survival Surv
#' @importFrom survminer ggsurvplot
#' @importFrom ggplot2 theme_bw
#' @importFrom stats median
#' @export

plt_surv <- function(SE, gene, method='Average_mean', style='elegant', conf.int=FALSE){
  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)
  meta <- bind_meta(SE, isList)

  Sc <- Core(exp_mtr, gene, method)
  Score <- as.numeric(ifelse(Sc>=median(Sc),1,0))

  time <- as.numeric(meta$overall.survival..days.)
  status <- sub('Dead','1', meta$vital.status) %>% sub('Alive','0',.)

  df <- data.frame(time,status,Score) %>%
    stats::na.omit() %>%
    lapply(as.numeric) %>%
    as.data.frame()

  return(surv_styling(df, style, conf.int, gene, method))
}


#' @title the style of survival plot
#' @description the style of survival plot
#' @param df the data source
#' @param style plotting style
#' @param conf.int logical value. If TRUE, plots confidence interval.
#' @param gene the gene you interested in.
#' @param method the method
#' @import ggplot2


surv_styling <- function(df, style, conf.int, gene, method){
  fit <- survfit(Surv(time, status) ~ Score, data = df)

  if(style == 'elegant'){
    cox_md <- coxph(Surv(time, status) ~ Score, data = df)
    summary_cox <- summary(cox_md)
    HR <- round(summary_cox$conf.int[,1],2)
    P_val_cox <- round(summary_cox$coefficients[,5],2)
    LCI <- round(summary_cox$conf.int[,3],2)
    UCI <- round(summary_cox$conf.int[,4],2)

    summary_KM <- summary(fit)
    idx_0 <- which(summary_KM$strata == 'Score=0')
    idx_1 <- which(summary_KM$strata == 'Score=1')
    median_0 <- summary_KM$time[which.min(summary_KM$surv[idx_0] - 0.5)]
    median_1 <- summary_KM$time[which.min(summary_KM$surv[idx_1] - 0.5) + length(idx_0)]
    LCI_0 <- summary_KM$time[which(summary_KM$lower[idx_0]<=0.5)[1]]
    LCI_1 <- summary_KM$time[which(summary_KM$lower[idx_1]<=0.5)[1] + length(idx_0)]
    UCI_0 <- summary_KM$time[which(summary_KM$upper[idx_0]<=0.5)[1]]
    UCI_1 <- summary_KM$time[which(summary_KM$upper[idx_1]<=0.5)[1]]
    P_val_KM <- round(survdiff(Surv(time, status) ~ Score, data = df)$pvalue,2)
  }

  lg.labs <- switch (style,
                     raw = c('Low','HIGH'),
                     elegant = c(paste0("Low median: ",median_0,
                                        " days (95% CI ",LCI_0,"-",UCI_0,")"),
                                 paste0("High median: ",median_1,
                                        " days (95% CI ",LCI_1,"-",UCI_1,")")),
                     brief = NULL)

  P <- ggsurvplot(fit, data = df, pval = FALSE, conf.int = conf.int,
                  palette = c('#15C4C8','#F28A72'),
                  legend.title = ifelse(length(gene) == 1,gene,"Score"),
                  legend.labs = lg.labs,
                  risk.table = TRUE,
                  risk.table.title = "Number at risk",
                  risk.table.col = "strata", linetype = "solid",
                  surv.median.line = "hv",
                  size = 1.2)

  if(style == 'raw'){
    P$plot <- P$plot +
      ggtitle("Survival analysis") +
      theme(plot.title = element_text(hjust = 0.5))

    if(length(gene) > 1)
      P$table$labels$y <- method
  } else{
    P$plot <-
      P$plot +
      theme(plot.title = element_text(hjust = 0,size = 12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.title = element_blank(),
            legend.position = switch (style,
                                      elegant = c(0.6,0.9),
                                      brief='none'),
            legend.background = element_rect(fill = "transparent"),
            axis.title.x = element_text(size = 12,face="bold"),
            axis.title.y = element_text(size = 12,face="bold"),
            axis.line = element_line(size = 1)) +
      labs(x="Overall survival (days)")

    P$table <-
      P$table +
      ggtitle("Numer at risk")+
      theme(plot.title = element_text(hjust = -0.1,vjust = 0.5,size = 12,face="bold"),
            text = element_text(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()) +
      scale_y_discrete(labels=c("Low","High"))
    P$table$theme$axis.text.y$hjust <- 0
  }

  if(style == 'elegant'){
    P$plot <-
      P$plot +
      annotate("text", x = 20, y = 0.2, label = paste0("HR ",HR),size = 3.5,hjust=0) +
      annotate("text", x = 20, y = 0.15, label = paste0("95% CI ",LCI,"-",UCI),size = 3.5,hjust=0) +
      annotate("text", x = 20, y = 0.1, label = paste0("P ",P_val_cox),size = 3.5,hjust=0)

    P$table <-
      P$table +
      theme(axis.line = element_line(size = 1))
  }

  if(style == 'brief'){
    P$table <-
      P$table +
      theme(axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_blank())
  }

  return(P)
}

