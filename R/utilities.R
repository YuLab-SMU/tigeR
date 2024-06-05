#' @title Process data before running machine learning algorithm
#' @description Process data before running machine learning algorithm
#' @param SE a SummarizedExperiment(SE) object or a list consists of multiple SE objects. The colData of the SE object(s) must contain treatment information named Treatment.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR if TRUE, classify patients with CR, MR, PR as Responders (R), and those with PD, SD, NR as Non-Responders(NR).
#' @param turn2HL if TRUE, the expression value of a gene is divided to "HIGH" or "LOW" based on its median expression.
#' @export

dataProcess <- function(SE, Signature, rmBE, response_NR, turn2HL){
  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)
  meta <- bind_meta(SE, isList)

  if(response_NR)
    meta$response <- response_standardize(meta$response)

  if(rmBE){
    exp_mtr <- rmBE(exp_mtr,meta)
    colnames(exp_mtr) <- rownames(meta)
  }

  idx <- response_filter(meta$response)
  idx_all <- seq_along(meta[,1])
  if(length(idx)==0){
    idx <- idx_all
  }else{
    idx <- idx_all[!idx_all %in% idx]
  }

  f <- dataPreprocess(exp_mtr, Signature, turn2HL, meta)
  if(turn2HL){
    exp_mtr <- f[[1]][,idx]
  }else{
    exp_mtr <- f[,idx,drop=FALSE]
  }
  meta <- meta[idx,]

  absent <- meta$response_NR=="UNK"

  if(turn2HL){
    thres <- f[[2]]
    return(list(exp_mtr[,!absent,drop=FALSE],meta[!absent,,drop=FALSE],thres))
  }

  return(list(exp_mtr[,!absent,drop=FALSE],meta[!absent,,drop=FALSE]))
}


#' @title Prepare expression matrix for down string analysis
#' @description dataPreprocess will remove missing genes. Then returns the sub-matrix of the genes whose SYMBOLs are in the signature.
#' @param exp_mtr An expression matrix which rownames are gene SYMBOL and colnames are sample ID.
#' @param Signature The aiming gene set(only Gene SYMBOL allowed).
#' @param turn2HL If TRUE, the expression value of a gene is divided to "HIGH" or "LOW" based on its median expression.
#' @param meta refers to the specific set of genes you wish to use for model construction.
#' @export

dataPreprocess <- function(exp_mtr, Signature = NULL, turn2HL = TRUE, meta = NULL){
  if(is.null(Signature))
    Signature <- rownames(exp_mtr)

  genes <- S4Vectors::intersect(Signature, rownames(exp_mtr))

  absent_genes <- length(Signature) - length(genes)
  if(length(Signature)>length(genes))
    message(paste0(absent_genes," Signature genes are not found in expression matrix. The function can execute properly, but the performance of the model may be compromised."))

  exp_mtr <- exp_mtr[genes,,drop=FALSE]
  rowname <- rownames(exp_mtr)
  colname <- colnames(exp_mtr)
  exp_mtr <- apply(exp_mtr, 2, as.numeric)
  if(is.vector(exp_mtr))
    exp_mtr <- matrix(exp_mtr, nrow = 1)

  rownames(exp_mtr) <- rowname

  idx_R <- which(meta$response_NR=="R")
  idx_N <- which(meta$response_NR=="N")

  exp_mtr <- t(apply(exp_mtr, 1, function(x, Batch) {
    for (i in unique(Batch)) {
      tmp <- x[which(Batch == i)]
      isNA <- is.na(tmp)
      if (all(isNA))
        next
      if (all(tmp[!isNA] == 0))
        tmp[!isNA] <- NA
      x[which(Batch == i)] <- tmp
    }
    x
  }, Batch = meta$batch))

  if(is.vector(exp_mtr))
    exp_mtr <- matrix(exp_mtr, nrow = 1)

  colnames(exp_mtr) <- colname

  if(turn2HL){
    threshold <-
      apply(exp_mtr,1,function(x){
        mean(mean(x[idx_R],na.rm=TRUE),mean(x[idx_N],na.rm=TRUE),na.rm=TRUE)
      })
    exp_mtr <- t(
      apply(exp_mtr,1,function(x){
      thres <- mean(mean(x[idx_R],na.rm=TRUE),mean(x[idx_N],na.rm=TRUE),na.rm=TRUE)
      x <- ifelse(x>=thres,'HIGH','LOW')

    }))
    if(is.vector(exp_mtr))
      exp_mtr <- matrix(exp_mtr, nrow = 1)
  }else{
    exp_mtr[is.na(exp_mtr)] <- 0
  }

  if(turn2HL)
    return(list(exp_mtr, threshold))

  return(exp_mtr)
}


#' @title judge whether all the elements in vector V are NA
#' @param V must be a vector.

is.NA_vec <- function(V){
  return(all(is.na(V)))
}

#' @title turn the 0 elements in vector V to NA
#' @param V must be a vector.

zero2na <- function(V){
  if(all(is.na(V)))
    return(V)
  if(any(is.na(V)) && all(V[!is.na(V)] == 0)){
    V[!is.na(V)] <- NA
    return(V)
  }else
    if(sum(V != 0) < 2)
      return(rep(NA, length(V))) #if only 1 or 2 element is non-zero number, return NA
  return(V)
}

#' @title Ranking features in matrix with Gini index
#' @description calculating the Gini index and get an overview of the classification efficiency of genes.
#' @param SE a SummarizedExperiment(SE) object or a list consists of multiple SE objects. The colData of the SE object(s) must contain treatment information named Treatment.
#' @importFrom stats setNames
#' @export

gini_gene <- function(SE){
  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)
  mtr <- dataPreprocess(exp_mtr,rownames(exp_mtr), turn2HL = TRUE)
  label <- bind_meta(SE, isList)$response_NR

  genes <- rownames(exp_mtr)
  features_Gini <- apply(mtr, 1, Gini, label = label)
  Gini <- setNames(features_Gini[genes], genes)

  return(Gini)
}

#' @title Ranking features in vector with Gini coefficient
#' @param vec An vector. Usually is one row of an matrix.
#' @param label The classification label of matrix.
#'

Gini <- function(vec, label){
  NR <- grep('NR|N', label)
  R <- grep('R', label)

  Gini_H <- Gini_internal(vec, R, "HIGH")
  Gini_L <- Gini_internal(vec, R, "LOW")
  Gini_Gene <- (sum(vec == 'HIGH')*Gini_H + sum(vec == 'LOW')*Gini_L)/length(label)
  return(Gini_Gene)
}

Gini_internal <- function(vec, index, category) {
  ii <- sum(vec[index] == category)
  tt <- sum(vec == category)

  1 - (ii/tt)^2 - (1 - ii/tt)^2
}

#' @title differential gene
#' @description return differential expression gene between Responder and Non-Responder.
#' @param SE a SummarizedExperiment(SE) object or a list consists of multiple SE objects. The colData of the SE object(s) must contain treatment information named Treatment.
#' @export

diff_gene <- function(SE){
  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)
  meta <- bind_meta(SE, isList)

  idx_R <- which(meta$response_NR == 'R')
  idx_N <- which(meta$response_NR == 'N')

  log2FC1 <- log2(apply(exp_mtr[,idx_R], 1, mean)/apply(exp_mtr[,idx_N], 1, mean)+1)
  P <- apply(exp_mtr, 1, matrix_t.test, P=idx_R, N=idx_N)
  Q <- -10 * log10(P)

  result <- data.frame(`log2(FC+1)`=log2FC1,p_value=P,q_value=Q)
  return(result)
}

#' @title differential gene
#' @description Return differential expression gene between Responder and Non-Responder
#' @param V the expression vector of a gene
#' @param P the index of Positive samples in vector V
#' @param N the index of Negative samples in vector V
#' @importFrom stats t.test

matrix_t.test <- function(V, P, N){
  return(t.test(V[P], V[N])$p.value)
}


#' @title Binding expression matrices from data folder in tigeR together
#' @description Extract expression data in particular data set or data sets from the data folder in tigeR. If there are more than one data set, this function will return an matrix which binds all the expression matrices by column.
#' @param datasetNames the name of data set or data sets you want to use.
#' @importFrom magrittr %>%
#' @importFrom rlang current_env
#' @export
#'

extract_mtr <- function(datasetNames){
  for (name in datasetNames) {
    if(!exists('inteMatrix', envir = current_env())){
      exp_mtr <- get(name) %>% assay()

      inteMatrix <- exp_mtr
      if(length(datasetNames > 1))
        next
    }
    exp_mtr <- get(name) %>% assay()

    inteMatrix <- cbind(inteMatrix, exp_mtr)
  }
  return(inteMatrix)
}


#' @title Binding response data from data folder in tigeR together
#' @description Extract response data in particular data set or data sets from the data folder in tigeR. If there are more than one data set, this function will return an vector which contains the response data of every data sets.
#' @param datasetNames the name of data set or data sets you want to use.
#' @importFrom magrittr %$%
#' @importFrom rlang current_env
#' @importFrom SummarizedExperiment colData
#' @export

extract_label <-function(datasetNames){
  for (name in datasetNames) {
    if(!exists('inteVector', envir = current_env())){
      get(name) %$% colData()$response_NR -> response

      inteVector <- response
      if(length(datasetNames > 1))
        next
    }
    get(name) %$% colData()$response_NR ->response

    inteVector <- c(inteVector, response)
  }

  inteVector[inteVector == 'UNK'] <- NA
  return(inteVector)
}

#' @title standardization of response labels
#' @description turn label to 'R' or 'NR'
#' @param V an vector
#' @export

response_standardize <- function(V){
  V <- sub('CR|MR|PR|CRPR', 'R', V)
  V <- sub('PD|SD|NR', 'N', V)
  return(V)
}


#' @title Max_Min normalization.
#' @description (x - min(x))/(max(x) - min(x))
#' @param exp_mtr an matrix which rows represent genes columns represent samples.
#' @export

max_min_normalization <- function(exp_mtr){
  mini <- apply(exp_mtr, 1, min)
  maxi <- apply(exp_mtr, 1, max)
  interval <- maxi - mini
  exp_mtr <- (exp_mtr - mini) / interval
  return(exp_mtr)
}


#' @title perform naive bayes prediction model.
#' @description Generate a naive bayes model.
#' @param SE a SummarizedExperiment(SE) object or a list consists of multiple SE objects. The colData of the SE object(s) must contain treatment information named Treatment.
#' @param isList whether SE is list
#' @importFrom SummarizedExperiment assay
#' @export

bind_mtr <- function(SE,isList){
  if (!isList){
    exp_mtr <- assay(SE)
  } else if (isList){
    exp_mtr <- matrix(unlist(lapply(SE, assay)),nrow = nrow(assay(SE[[1]])))
    rownames(exp_mtr) <- rownames(assay(SE[[1]]))
    colnames(exp_mtr) <- unlist(lapply(SE,colnames))
  }
  return(exp_mtr)
}


#' @title perform naive bayes prediction model.
#' @description Generate a naive bayes model.
#' @param SE a SummarizedExperiment(SE) object or a list consists of multiple SE objects. The colData of the SE object(s) must contain treatment information named Treatment.
#' @param isList whether SE is list
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment colData
#' @export

bind_meta <- function(SE,isList){
  if (!isList){
    meta <- as.data.frame(colData(SE))
  } else if (isList){
    if(length(SE) == 1)
      return(as.data.frame(colData(SE[[1]])))
    meta <- as.data.frame(colData(SE[[1]]))
    meta$batch <- rep('batch1',nrow(meta))
    for (i in 2:length(SE)) {
      paste0('batch', i) %>% rep(nrow(colData(SE[[i]]))) -> batch
      colData(SE[[i]]) %>% as.data.frame() %>% cbind(batch) %>% rbind(meta,.) -> meta
    }
  }
  return(meta)
}


#' @title Remove batch effect.
#' @description Generate a naive bayes model.
#' @param mtr an expression matrix.
#' @param meta meta informations.
#' @importFrom stats model.matrix

rmBE <- function(mtr, meta){
  batch <- unique(meta$dataset_id)
  if("batch" %in% colnames(meta)){
    model <- model.matrix(~as.factor(meta$response))
    mtr <- dataPreprocess(mtr, rownames(mtr), turn2HL = FALSE)
    mtr <- sva::ComBat(dat = mtr,batch = as.factor(meta$batch),mod = model)
  }else if(length(batch)>1){
    model <- model.matrix(~as.factor(meta$response))
    mtr <- dataPreprocess(mtr, rownames(mtr), turn2HL = FALSE)
    mtr <- sva::ComBat(dat = mtr,batch = as.factor(meta$dataset_id),mod = model)
  }
  return(mtr)
}


#' @title Filting missing response value.
#' @description return the index of NE or UNK in response vector.
#' @param response a vector which contains response information.
#' @importFrom magrittr %>%
#' @export

response_filter <- function(response){
  idx <- grep('NE|UNK',response)
  return(idx)
}


#' @title Build plot theme
#' @description return the ploting theme
#' @param df a dataframe
#' @param textcol the color of the text in the plot.
#' @param panelcol the color of the panel border and ticks in the plot.
#' @import ggplot2
#' @importFrom rlang .data

plt_style <- function(df, textcol, panelcol){
  diff_theme <- theme(plot.background = element_rect(color="transparent",
                                                     fill="transparent"),
                      plot.title = element_text(face = "bold",size = "14", color = textcol),
                      axis.title = element_text(face = "bold", size = "12", color = textcol),
                      axis.text = element_text(face = "bold", size = "10", color = textcol),
                      axis.ticks = element_line(color = panelcol),
                      panel.background = element_rect(fill = "transparent"),
                      panel.border = element_rect(linewidth = 1.5, fill="transparent", color = panelcol),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      legend.background = element_rect(fill="transparent"),
                      legend.position = "right",
                      legend.title = element_text(face = "bold", size = "12",color = textcol),
                      legend.text = element_text(face='bold', size='10',color=textcol),
                      legend.key = element_blank(),
                      aspect.ratio = 1)
  df$group <- sub("Non-Responder","NR",df$group)
  df$group <- sub("Responder","R",df$group)
  df$group <- sub("Post-Therapy","Post",df$group)
  df$group <- sub("Pre-Therapy","Pre",df$group)
  mycolor <- c("#5f96e8","#ee822f")
  names(mycolor) <- unique(df$group)
  if(all(df$group %in% c("Pre","Post")))
    df$group <- factor(df$group,levels = c("Pre","Post"))

  ggplot(df, aes(x = .data$group,
                 y = .data$Score,
                 color = .data$group)) +
    scale_color_manual(values = mycolor) +
    geom_boxplot(lwd=1.2,outliers=FALSE) +
    geom_jitter(aes(fill = .data$group),width = 0.2, size = 1.5) +
    diff_theme
}


#' @title Prepare data for plot
#' @description Preparing data for ploting.
#' @param gene is the Gene or Gene set you are interested in.
#' @param SE a SummarizedExperiment(SE) object or a list consists of multiple SE objects. The colData of the SE object(s) must contain treatment information named Treatment.
#' @param method the method for calculating gene set scores. Can be NULL if the length of parameter gene is 1.
#' @param type the type of information
#' @param PT_drop If TRUE, only Untreated patient will be use for model training.
#' @param log_sc if TRUE, log(value + 1)
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%

plt_Preprocess <- function(gene, SE, method, type, PT_drop, log_sc){
  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)
  meta <- bind_meta(SE, isList)

  if(type == 'R vs NR'){
    if(PT_drop){
      idx_UT <- which(meta$Treatment == 'PRE')
      if(length(idx_UT) == 0)
        stop("All patients in data set have been treated. Setting the parameter PT_drop to FALSE to run anyway.")
      meta <- meta[idx_UT,,drop=FALSE]
      exp_mtr <- exp_mtr[,idx_UT,drop=FALSE]
    }

    group <- as.vector(meta$response_NR)
  }
  if(type == 'T vs UT')
    group <- as.vector(meta$Treatment)

  Sc <- Core(exp_mtr, gene, method)
  Score <- Sc
  if(log_sc)
    Score <- log2(Sc + 1)

  df <- data.frame(group,Score,Sc)
  idx <- response_filter(df$group)
  if(length(idx) != 0)
    df <- df[-idx,]

  if(type == 'R vs NR')
    df$group %<>% sub('R','Responder',.) %>% sub('N','Non-Responder',.)
  if(type == 'T vs UT')
    df$group %<>% sub('PRE','Pre-Therapy',.) %>% sub('ON|EDT|POST*','Post-Therapy',.)

  return(df)
}


#' @title Prepare data for plot
#' @description Preparing data for ploting.
#' @param ROC the ROC object
#' @param auc.pos the position of the AUC value
#' @param auc.round the decimal places you want to keep for auc value
#' @param textcol the color of the text in the plot.
#' @param panelcol the color of the panel border and ticks in the plot.

plt_roc <- function(ROC,auc.pos,auc.round,textcol,panelcol){
    pROC::ggroc(ROC, color = "black", size = 1) +
    ggplot2::annotate("segment", x = 0, xend = 1, y = 1, yend = 0,
                      color = "#646464", size = 0.5, linetype = "solid") +
    ggplot2::annotate("text",x = auc.pos[1], y = auc.pos[2],
                      label = ifelse(ROC$auc < 0.1^auc.round,
                                     paste0("AUC < ",format(0.1^auc.round,scientific=FALSE)),
                                     paste0("AUC = ",sprintf(paste0("%.",auc.round,"f"),ROC$auc))),
                      size = 4.5,color = textcol) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.background = element_rect(color="transparent",
                                                  fill="transparent"),
                   plot.title = element_text(face = "bold",size = "14", color = textcol, hjust = 0.5),
                   axis.title = element_text(face = "bold", size = "12", color = textcol),
                   axis.text = element_text(face = "bold", size = "9", color = textcol),
                   axis.ticks = element_line(color = panelcol),
                   panel.background = element_rect(fill = "transparent"),
                   panel.border = element_rect(fill = "transparent", color = panelcol, linewidth = 1.5),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   legend.position = "right",
                   legend.title = element_text(face = "bold", size = "12",color = textcol),
                   legend.text = element_text(face='bold',
                                              size='9',color=textcol),
                   legend.key = element_blank(),
                   aspect.ratio = 1)
}

#' @title count geneset score by different method
#' @description wait to write
#' @param exp_mtr an expression matrix.
#' @param geneSet the gene or geneset which you wanted to investigate.
#' @param method the method for calculating gene set scores. Can be NULL if the length of parameter gene is 1.

Core <- function(exp_mtr, geneSet, method){
  if(is.null(method)){
    if(is.numeric(geneSet)){
      method <- "Weighted_mean"
    }else if(is.character(geneSet)){
      method <- "Average_mean"
    }else{return(exp_mtr[geneSet,])}
  }


  if(length(geneSet)==1){
    exist_genes <- intersect(rownames(exp_mtr), geneSet)
    if(length(exist_genes)==1)
      return(exp_mtr[geneSet,])
    else
      stop(paste0(geneSet, " does not present in expression matrix."))
  }

  if(method == "Weighted_mean"){
    if(!is.numeric(geneSet)){
      stop("If argument 'method' is 'Weighted_mean', the Signature gene set must be a numeric vector with gene names.")
    }
    geneSet0 <- geneSet
    geneSet <- names(geneSet)
  }

  exist_genes <- intersect(rownames(exp_mtr), geneSet)

  if(length(geneSet) != length(exist_genes)){
    str <- ""
    for (g in geneSet[!geneSet %in% exist_genes]) {
      str <- paste(str, g, sep = " ")
    }
    message(str, " does not exist in expression matrix.")
  }

  cn <- colnames(exp_mtr)
  exp_mtr <- stats::na.omit(t(apply(exp_mtr, 1, function(x){
    if(all(is.na(x)|x == 0)){
      return(rep(NA,length(x)))
    }else{
      return(x)
    }
  })))
  colnames(exp_mtr) <- cn

  exist_genes <- intersect(rownames(exp_mtr), geneSet)

  Score <-
    switch(method,
           Average_mean = apply(exp_mtr[exist_genes,,drop=FALSE], 2, mean),
           GSVA = GSVA::gsvaParam(exp_mtr, geneSets=list(exist_genes)) %>% GSVA::gsva(),
           Weighted_mean = weight_mean_signature(exp_mtr, geneSet0[names(geneSet0) %in% exist_genes]))

  Score <- as.vector(Score)
  names(Score) <- colnames(exp_mtr)
  return(Score)
}

globalVariables(".")
