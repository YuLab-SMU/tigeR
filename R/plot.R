#' @title Calculate the correlations between user defined gene or gene set and TCGA cancerous samples.
#' @description Generating an Heatmap showing the pearson correlation coefficient.
#' @param gene is the Gene or Gene set you are interested in.
#' @import utils
#' @import rlang
#' @importFrom stats cor
#'

geneCorr <- function(gene='CD274'){
  projects <- c('TCGA-ESCA','TCGA-SARC','TCGA-CESC','TCGA-UCEC','TCGA-LAML','TCGA-TGCT','TCGA-CHOL','TCGA-MESO','TCGA-ACC','TCGA-DLBC','TCGA-PCPG','TCGA-KICH','TCGA-THCA','TCGA-THYM','TCGA-UCS','TCGA-GBM','TCGA-UVM','TCGA-COAD','TCGA-STAD','TCGA-OV','TCGA-KIRC','TCGA-LGG','TCGA-HNSC','TCGA-BLCA','TCGA-LUAD','TCGA-PRAD','TCGA-LIHC','TCGA-LUSC','TCGA-SKCM','TCGA-KIRP','TCGA-BRCA','TCGA-PAAD','TCGA-READ')
  signatures <- 0
  data(signatures, package = 'tigeR', envir = current_env())
  final_cor_list <- list()
  for (project in projects) {
    cancertype <- strsplit(project, '-')[[1]][2]
    exp_mtr <- 0
    data(list = paste0('TCGA_', cancertype, '_exp'), package = 'tigeR', envir = current_env(), overwrite = TRUE)

    if (!any(gene %in% rownames(exp_mtr))) {
      final_cor_list[[project]] <- rep(0, length(signatures))
      next
    } else
      genes_exp <- exp_mtr[rownames(exp_mtr) %in% gene, ]

    if (!is.null(nrow(genes_exp))) {
      genes_exp <- unlist(apply(genes_exp, 2, mean))
    }
    project_cor <- c()
    for (signature in names(signatures)) {
      geneset <- unlist(signatures[[signature]])
      geneset_mtr <- exp_mtr[rownames(exp_mtr) %in% c(geneset), ]

      if (length(as.vector(geneset_mtr)) == 0) {
        #if no element in matrix
        project_cor[length(project_cor) + 1] <- 0
        next
      } else if (!is.null(dim(geneset_mtr))) {
        #if more than one row in matrix
        sig_exp <- apply(geneset_mtr, MARGIN =  2, mean)
        tmp_mtr <- t(rbind(genes_exp, sig_exp))
        project_cor[length(project_cor) + 1] <- cor(tmp_mtr)[2, 1]
      } else{
        #if only one row in matrix
        sig_exp <- geneset_mtr
        tmp_mtr <- t(rbind(genes_exp, sig_exp))
        project_cor[length(project_cor) + 1] <- cor(tmp_mtr)[2, 1]
      }
    }
    final_cor_list[[project]] <- project_cor
  }

  final_cor <- rbind(final_cor_list[[1]], final_cor_list[[2]])
  for (i in 3:length(final_cor_list)) {
    final_cor <- rbind(final_cor, final_cor_list[[i]])
  }
  final_cor <- cbind(rep(1, length(projects)), final_cor)

  rownames(final_cor) <- names(final_cor_list)
  colnames(final_cor) <- c('Customed Geneset', names(signatures))
  if (all(final_cor[, -1] == 0))
    message(
      'The gene or gene set you input may not exist in the TCGA expression matrix! Please enter correct gene or gene set!'
    )
  return(final_cor)
}


#' @title plot differential result(Responder vs Non-Responder or Pre-Treatment vs Post-Treatment)
#' @description plot differential result(Responder vs Non-Responder or Pre-Treatment vs Post-Treatment)
#' @param gene is the Gene or Gene set you are interested in.
#' @param SE SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param type 'Treatment' or 'Response'.the type of analysis you want to perform(Responder vs Non-Responder or Pre-Treatment vs Post-Treatment)
#' @export

plt_diff <- function(gene='CD274',SE,type){
  if(type == 'Response')
    df <- plt_Preprocess(gene,SE,'R vs NR')
  if(type == 'Treatment')
    df <- plt_Preprocess(gene,SE,'T vs UT')

  plt_style(df)
}


#' @title plot differential result(Pre-Treatment vs Post-Treatment)
#' @description The association between gene expression and overall survival in the immunotherapy data was calculated using univariate Cox regression analysis.
#' @param gene is the Gene or Gene set you are interested in.
#' @param SE SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
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

plt_surv <- function(gene='CD274',SE){
  exp <- assay(SE)[rownames(SE) == gene,]
  exp <- as.numeric(ifelse(exp>=median(exp),1,0))

  time <- as.numeric(SE@colData$overall.survival..days.)
  sub('Dead','1',SE@colData$vital.status) %>% sub('Alive','0',.) -> status

  data.frame(time,status,exp) %>% na.omit() %>% lapply(as.numeric) %>% as.data.frame() -> df

  fit <- survfit(Surv(time, status) ~ exp, data = df)

  ggsurvplot(fit,
             data = df,
             pval = TRUE,
             conf.int = TRUE,
             legend.title = gene,
             legend.labs = c('Low','HIGH'),
             risk.table = TRUE,
             risk.table.title = 'Number at risk',
             risk.table.col = "strata",
             linetype = "strata",
             surv.median.line = "hv",
             ggtheme = theme_bw())
}
