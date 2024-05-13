#' @title perform survival analysis
#' @description calculates hazard ratios, confidence intervals and P value of cox-ph analysis as well as draw KM curve.
#' @param SE SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param gene is the Gene or Gene set you are interested in.
#' @param method the method for calculating gene set scores which has several options: Average_mean, Weighted_mean, or GSVA. The method can be set to NULL if the length of the parameter geneSet is 1. This means that if you are working with only one gene, the specific calculation method may not be applicable or necessary.
#' @param style the ploting style. c("raw,"elegant","brief")
#' @param conf.int logical value. If TRUE, plots confidence interval.
#' @param val.pos the position of annotation value.
#' @param lg dd
#' @param lg.pos the position of legend. When lg.pos=c(0,0), the legend will be placed at the leftdown of the plot.
#' @param lg.text c("precise","all"). If all, show the 0.95 CI of median.
#' @param p.round the decimal places you want to keep for p value
#' @param PT_drop If TRUE, only Untreated patient will be use for model training.
#' @param drop_zero If TRUE, 0 will not be consider when calculate median.
#' @import ggplot2
#' @importFrom SummarizedExperiment assay
#' @importFrom magrittr %>%
#' @importFrom survival survfit
#' @importFrom survival Surv
#' @importFrom survminer ggsurvplot
#' @importFrom ggplot2 theme_bw
#' @importFrom stats median
#' @export

surv_biomk <- function(SE, gene, method='Average_mean', style='elegant', conf.int=FALSE, val.pos=c(0.03,0.2), lg=NULL,lg.pos=c(0.6,0.9), lg.text="precise",p.round=4,PT_drop=TRUE,drop_zero=FALSE){
  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)
  meta <- bind_meta(SE, isList)

  idx_UT <- seq_along(meta[,1])
  if(PT_drop)
    idx_UT <- which(meta$Treatment == 'PRE')

  if(length(idx_UT) == 0)
    stop("Only Untreated patients can be use to perform survival analysis!")

  exp_mtr <- exp_mtr[,idx_UT,drop=FALSE]
  meta <- meta[idx_UT,,drop=FALSE]

  Sc <- Core(exp_mtr, gene, method)
  if(drop_zero)
    Score <- as.numeric(ifelse(Sc>=median(Sc[Sc!=0]),1,0))
  else{
    Score <- as.numeric(ifelse(Sc>=median(Sc),1,0))
  }

  time <- as.numeric(meta$overall.survival..days.)
  status <- sub('Dead','1', meta$vital.status) %>% sub('Alive','0',.)

  df <- data.frame(time,status,Score) %>%
    stats::na.omit() %>%
    lapply(as.numeric) %>%
    as.data.frame()

  return(surv_styling(df, style, conf.int, gene, method, val.pos,lg, lg.pos, lg.text,p.round))
}


#' @title the style of survival plot
#' @description the style of survival plot
#' @param df the data source
#' @param style plotting style
#' @param conf.int logical value. If TRUE, plots confidence interval.
#' @param gene the gene you interested in.
#' @param method the method
#' @param val.pos the position of annotation value.
#' @param lg dd
#' @param lg.pos the position of legend.
#' @param lg.text c("precise","all"). If all, show the 0.95 CI of median.
#' @param p.round the decimal places you want to keep for p value
#' @import ggplot2
#' @export

surv_styling <- function(df, style, conf.int, gene, method, val.pos, lg,lg.pos, lg.text, p.round){
  fit <- survfit(Surv(time, status) ~ Score, data = df)

  if(style == 'elegant'){
    val.pos <- c(val.pos[1]*1000, val.pos[2]-0.1)
    cox_md <- coxph(Surv(time, status) ~ Score, data = df)
    summary_cox <- summary(cox_md)
    HR <- sprintf("%.2f",summary_cox$conf.int[,1])
    P_val_cox <- sprintf(paste0("%.",p.round,"f"),summary_cox$coefficients[,5])
    LCI <- round(summary_cox$conf.int[,3],2)
    UCI <- round(summary_cox$conf.int[,4],2)

    summary_KM <- summary(fit)
    idx_0 <- which(summary_KM$strata == 'Score=0')
    idx_1 <- which(summary_KM$strata == 'Score=1')
    median_0 <- round(summary_KM$time[which.min(summary_KM$surv[idx_0] - 0.5)],1)
    median_1 <- round(summary_KM$time[which.min(summary_KM$surv[idx_1] - 0.5) + length(idx_0)],1)
    LCI_0 <- round(summary_KM$time[which(summary_KM$lower[idx_0]<=0.5)[1]],1)
    LCI_1 <- round(summary_KM$time[which(summary_KM$lower[idx_1]<=0.5)[1] + length(idx_0)],1)
    UCI_0 <- round(summary_KM$time[which(summary_KM$upper[idx_0]<=0.5)[1]],1)
    UCI_1 <- round(summary_KM$time[which(summary_KM$upper[idx_1]<=0.5)[1]],1)

    P_val_KM <- round(survival::survdiff(Surv(time, status) ~ Score, data = df)$pvalue,p.round)
    if(P_val_cox < 0.001)
      P_val_cox <- "< 0.001"
    else{
      P_val_cox <- paste0("= ",P_val_cox)
    }
  }

  if(lg.text == "precise")
    lg.title <- c(paste0("Low median: ",median_0," days"),
                  paste0("High median: ",median_1," days"))
  if(lg.text == "all")
    lg.title <- c(paste0("Low median: ",median_0,
                         " days (95% CI ",LCI_0,"-",UCI_0,")"),
                  paste0("High median: ",median_1,
                         " days (95% CI ",LCI_1,"-",UCI_1,")"))
  if(lg.text == "specific")
    lg.title <- c(paste0(ifelse(length(gene)==1,gene,lg)," Low"),
                  paste0(ifelse(length(gene)==1,gene,lg)," High"))

  lg.labs <- switch (style,
                     raw = c('Low','High'),
                     elegant = lg.title,
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
      theme(plot.title = element_text(face = "bold", size = "14",
                                      color = "black", hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.title = element_blank(),
            legend.position = switch(style,
                                     elegant = lg.pos,
                                     brief='none'),
            legend.background = element_rect(fill = "transparent"),
            axis.title.x = element_text(size = 12,face="bold"),
            axis.title.y = element_text(size = 12,face="bold"),
            axis.line = element_line(linewidth = 1)) +
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
      annotate("text", x=val.pos[1], y=val.pos[2]+0.2, label = paste0("HR ",HR),size = 3.5,hjust=0) +
      annotate("text", x=val.pos[1], y=val.pos[2]+0.1, label = paste0("95% CI ",LCI,"-",UCI),size = 3.5,hjust=0) +
      annotate("text", x=val.pos[1], y=val.pos[2], label = paste0("P ",P_val_cox),size = 3.5,hjust=0)

    P$table <-
      P$table +
      theme(axis.line = element_line(linewidth = 1))
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

