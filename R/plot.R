#' @title Calculate the correlations between user defined gene or gene set and TCGA cancerous samples.
#' @description Generating an Heatmap showing the pearson correlation coefficient.
#' @param gene is the Gene or Gene set you are interested in.
#' @import utils
#' @import rlang
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


#' @title Calculate the association between gene expression and overall survival in the immunotherapy data.
#' @description The association between gene expression and overall survival in the immunotherapy data was calculated using univariate Cox regression analysis.
#' @param gene is the Gene or Gene set you are interested in.
#' @param type is the analysis you want to perform.'cox' or 'auc'
#' @import survival
#' @import stats
#' @import utils
#' @import rlang
#' @importFrom pROC roc
#'

geneSurv <- function(gene='CD274',type='cox') {
  if(!type %in% c('cox', 'auc')){
    message("Incorrect parameter 'type'. 'type' should be 'cox' or 'auc'!")
    return()
  }
  datasets <- c("GBM-PRJNA482620",
                "HNSC-GSE93157",
                "LGG_E-MTAB-6270",
                "LUSC-GSE93157",
                "Melanoma-GSE100797",
                "Melanoma-GSE106128",
                "Melanoma-GSE115821",
                "Melanoma-GSE78220",
                "Melanoma-GSE91061",
                "Melanoma-GSE93157",
                "Melanoma-GSE96619",
                "Melanoma-Nathanson_2017",
                "Melanoma-phs000452",
                "Melanoma-PRJEB23709",
                "Melanoma_GSE145996",
                "nonsqNSCLC-GSE93157",
                "NSCLC_GSE126044",
                "NSCLC_GSE135222",
                "RCC-Braun_2020",
                "RCC-GSE67501",
                "STAD-PRJEB25780")
  signatures <- 0
  data(signatures, package = 'tigeR', envir = current_env(), overwrite = TRUE)

  final_sur_list <- list()
  for (dataset in datasets) {
    exp <- 0
    data(list = as.character(paste0(dataset, '.Response')), envir = current_env(), overwrite = TRUE)
    exp_mtr <- matrix(data = unlist(exp),
                      ncol = length(exp),
                      byrow = FALSE)[, -1]
    exp_mtr <- matrix(data = as.numeric(exp_mtr), nrow = nrow(exp_mtr))
    rownames(exp_mtr) <- exp[[1]]
    colnames(exp_mtr) <- names(exp)[-1]
    exp_mtr <- na.omit(exp_mtr)

    meta <- 0
    data(list = paste0(dataset, '.meta'), package = 'tigeR', envir = current_env(), overwrite = TRUE)
    meta <- meta[meta$sample_id %in% colnames(exp_mtr), ]
    colnames(meta)[colnames(meta) == "overall.survival..days."] <- 'time'
    colnames(meta)[colnames(meta) == "vital.status"] <- 'event'
    meta$event <- ifelse(meta$event == 'Dead', 1, 0)
    if (any(is.na(meta$event))) {
      final_sur_list[[dataset]] <- rep(0, 12)
      next
    }

    Sur_score <- c()
    if (!any(gene %in% rownames(exp_mtr))) {
      Sur_score[length(Sur_score) + 1] <- 0
    } else{
      genes_exp <- exp_mtr[rownames(exp_mtr) %in% gene,]
      if (!is.null(nrow(genes_exp)))
        genes_score <- unlist(apply(genes_exp, 2, mean))
      else
        genes_score <- genes_exp
      if(type == 'cox'){
        meta$gene <- ifelse(genes_score > median(genes_score), 'High', 'Low')
        res.cox <- coxph(Surv(time, event) ~ gene, data = meta)
        res.cox.sum <- summary(res.cox)
        Sur_score[length(Sur_score) + 1] <- -sign(log2(exp(res.cox.sum$coefficients[, 1]))) * log10(res.cox.sum$coefficients[, 5])
      }else if(type == 'auc'){
        meta$gene <- as.numeric(genes_score)
        res.auc <- roc(meta$event, meta$gene, smooth = F, ci = T, auc = T)
        Sur_score[length(Sur_score) + 1] <- res.auc$auc
      }
    }

    for (signature in names(signatures)) {
      Geneset <- unlist(signatures[[signature]])
      Geneset_mtr <- exp_mtr[rownames(exp_mtr) %in% Geneset, ]
      if (length(as.vector(Geneset_mtr)) == 0) {
        Sur_score[length(Sur_score) + 1] <- 0
        next
      }
      else if (!is.null(dim(Geneset_mtr))) {
        Sig_score <- apply(Geneset_mtr, MARGIN =  2, mean)
      } else{
        Sig_score <- Geneset_mtr
      }
      if(type == 'cox'){
        meta$gene <- ifelse(Sig_score > median(Sig_score), 'High', 'Low')
        res.cox <- coxph(Surv(time, event) ~ gene, data = meta)
        res.cox.sum <- summary(res.cox)
        Sur_score[length(Sur_score) + 1] <- -sign(log2(exp(res.cox.sum$coefficients[, 1]))) * log10(res.cox.sum$coefficients[, 5])
      }else if(type == 'auc'){
        meta$gene <- as.numeric(Sig_score)
        res.auc <- pROC::roc(meta$event, meta$gene, smooth = F, ci = T, auc = T)
        Sur_score[length(Sur_score) + 1] <- res.auc$auc
      }
    }
    final_sur_list[[dataset]] <- Sur_score
  }

  final_mtr <- cbind(final_sur_list[[1]], final_sur_list[[2]])
  for (i in 3:length(final_sur_list)) {
    final_mtr <- cbind(final_mtr, final_sur_list[[i]])
  }

  emptycol <-c()
  for (i in 1:dim(final_mtr)[2]) {
    if(all(final_mtr[,i] == 0))
      emptycol <- c(emptycol, i)
  }
  final_mtr <- final_mtr[,-emptycol]

  rownames(final_mtr) <- c('Customed Geneset', names(signatures))
  colnames(final_mtr) <- names(final_sur_list)[-emptycol]
  if (all(final_mtr == 0))
    message(
      'The gene or gene set you input may not exist in the TCGA expression matrix! Please enter correct gene or gene set!'
    )
  return(final_mtr)
}


#' @title plot differential result(responder vs nonresponder)
#' @description The association between gene expression and overall survival in the immunotherapy data was calculated using univariate Cox regression analysis.
#' @param gene is the Gene or Gene set you are interested in.
#' @param SE SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @importFrom SummarizedExperiment assay
#' @import ggplot2
#' @import ggpubr
#' @export

plt_RvsNR <- function(gene='CD274',SE){
  exp <- assay(SE)[rownames(SE) == gene,]
  exp <- log2(exp + 1)

  group <- as.vector(SE@colData$response_NR)
  df <- data.frame(group,exp)
  idx <- which(df$group == 'UNK')
  df <- df[-idx,]
  df$group %<>% sub('N','Non-Responder(NR)',.) %>% sub('R','Responder(R)',.)

  mytheme <- ggplot2::theme(plot.title=element_text(face='bold',
                                           size='14',color='black'),
                   axis.title=element_text(face='bold',
                                           size='14',color='black'),
                   axis.text=element_text(face='bold',
                                          size='9',color='black'),
                   panel.background=element_rect(fill='white',color='black',
                                                 size=1.3),
                   legend.position='right',
                   legend.title =element_text(face='bold',
                                              size='14',color='black'))

  plot <- ggplot(df, aes(x=group,y=exp,color=group)) +
    geom_boxplot()+
    geom_jitter(aes(fill=group),width =0.2,shape = 21,size=1)+
    mytheme+
    labs(title='ALL',x=NULL,y='Gene Expression(log2(FPKM + 1))')
  return(plot)
}


#' @title plot differential result(Pre-Treatment vs Post-Treatment)
#' @description The association between gene expression and overall survival in the immunotherapy data was calculated using univariate Cox regression analysis.
#' @param gene is the Gene or Gene set you are interested in.
#' @param SE SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @importFrom SummarizedExperiment assay
#' @import ggplot2
#' @import ggpubr
#' @export

plot_TvsUT <- function(gene='CD274',SE){
  exp <- assay(MEL_GSE91061)[rownames(MEL_GSE91061) == gene,]
  exp <- log2(exp + 1)

  group <- as.vector(SE@colData$Treatment)
  df <- data.frame(group,exp)
  df$group %<>% sub('PRE','Pre-Therapy',.) %>% sub('ON','Post-Therapy',.)

  mytheme <- theme(plot.title=element_text(face='bold',
                                           size='14',color='black'),
                   axis.title=element_text(face='bold',
                                           size='14',color='black'),
                   axis.text=element_text(face='bold',
                                          size='9',color='black'),
                   panel.background=element_rect(fill='white',color='black',
                                                 size=1.3),
                   legend.position='right',
                   legend.title =element_text(face='bold',
                                              size='14',color='black'))

  plot <- ggplot(df, aes(x=group,y=exp,color=group)) +
    geom_boxplot() +
    geom_jitter(aes(fill=group),width =0.2,shape = 21,size=1) +
    mytheme +
    labs(title='ALL',x=NULL,y='Gene Expression(log2(FPKM + 1))')
  return(plot)
}
