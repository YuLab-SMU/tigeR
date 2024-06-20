#' @title Visualize differential analysis result (Responder vs NonResponder or Pre-Treatment vs Post-Treatment).
#' @description Visualize differential analysis result (Responder vs NonResponder or Pre-Treatment vs Post-Treatment).
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param gene the Gene or Gene set you are interested in.
#' @param type the comparison group you want to choose, 'Treatment' (Pre-Treatment vs Post-Treatment) or 'Response' (Responder vs Non-Responder ).
#' @param method the method for calculating gene set scores which has several options: Average_mean, Weighted_mean, or GSVA. The method can be set to NULL if the length of the parameter geneSet is 1. This means that if you are working with only one gene, the specific calculation method may not be applicable or necessary.
#' @param PT_drop If TRUE, only Untreated patient will be use for model training.
#' @param p.pos the position of the P value. When p.pos=c(0,0), the legend will be placed in the bottom left corner of the plot.
#' @param p.round the decimal places you want to keep for p value
#' @param log_sc if TRUE, log(value + 1)
#' @param textcol the color of the text in the plot.
#' @param panelcol the color of the panel border and ticks in the plot.
#' @return
#'   \describe{
#'   Return a bar plot visualizing the differential analysis (Responder vs NonResponder or Pre-Treatment vs Post-Treatment)}
#' @examples
#' diff_biomk(SE=MEL_GSE78220,gene='CD274',type='Treatment')
#' @export

diff_biomk <- function(SE, gene, type, method='Average_mean', PT_drop=TRUE, p.pos=c(0.2,0.7),p.round=2,log_sc=TRUE, textcol="black", panelcol="black"){
  p.pos[2] <- p.pos[2]*10
  p.pos <- p.pos + c(0.4,0)
  type <- match.arg(type, c('Response','Treatment'))
  method <- match.arg(method, c('Average_mean','GSVA','Weighted_mean'))

  if(type == 'Response'){
    df <- plt_Preprocess(gene, SE, method, 'R vs NR', PT_drop,log_sc)
    P <- stats::wilcox.test(df[df[,1] == "Responder",3],df[df[,1] == "Non-Responder",3])$p.value
    plt <- plt_style(df,textcol = textcol, panelcol = panelcol) + ggplot2::ggtitle("Responder vs Non-Responder")
  }else{
    df <- plt_Preprocess(gene, SE, method, 'T vs UT', FALSE,log_sc)
    P <- stats::wilcox.test(df[df[,1] == "Post-Therapy",3],df[df[,1] == "Pre-Therapy",3])$p.value
    plt <- plt_style(df,textcol = textcol, panelcol = panelcol) + ggplot2::ggtitle("Treatment vs UnTreatment")
  }
  if(P < 0.1^p.round)
    P <- paste0("< ",format(0.1^p.round,scientific=FALSE))
  else{
    P <- paste0("= ",sprintf(paste0("%.",p.round,"f"),P))
  }
  return(plt +
           ggplot2::labs(title=NULL,x=NULL,y=ifelse(length(gene)==1,
                                                    paste0(gene," expression"),
                                                    "Signature Score")) +
           ggplot2::theme(plot.title = element_text(hjust = 0.5)) +
           annotate("text", x=p.pos[1], y=p.pos[2], label = paste0("P ",P),size = 3.5,hjust=0,color=textcol))
}
