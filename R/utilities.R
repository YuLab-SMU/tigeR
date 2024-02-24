#' @title Process data before running machine learning algorithm
#' @description Process data before running machine learning algorithm
#' @param SE a SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param Signature an gene set you interested in
#' @param rmBE whether remove batch effect between different data set using internal Combat method
#' @param response_NR If TRUE, only use R or NR to represent Immunotherapy response of patients.
#' @param turn2HL If TRUE, the expression value of a gene is divided to "HIGH" or "LOW" based on its median expression.
#' @export

dataProcess <- function(SE, Signature, rmBE, response_NR, turn2HL){
  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)
  meta <- bind_meta(SE, isList)

  if(response_NR)
    meta$response <- response_standardize(meta$response)

  if(rmBE && isList){
    exp_mtr <- rmBE(exp_mtr,meta)
    colnames(exp_mtr) <- rownames(meta)
  }

  idx <- response_filter(meta$response)
  idx_all <- 1:nrow(meta)
  if(length(idx)==0){
    idx <- idx_all
  }else{
    idx <- idx_all[!idx_all %in% idx]
  }

  f <- dataPreprocess(exp_mtr, Signature, turn2HL, meta)
  if(turn2HL){
    exp_mtr <- f[[1]][,idx]
  }else{
    exp_mtr <- f[,idx]
  }
  meta <- meta[idx,]

  absent <- meta$response_NR=="UNK"

  if(turn2HL){
    thres <- f[[2]]
    return(list(exp_mtr[,!absent],meta[!absent,],thres))
  }

  return(list(exp_mtr[,!absent],meta[!absent,]))
}


#' @title Prepare expression matrix for down string analysis
#' @description dataPreprocess() will remove missing genes. Then returns the sub-matrix of the genes whose SYMBOLs are in the signature.
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
    message(paste0(absent_genes," Signature genes are not found in expression matrix."))

  exp_mtr <- exp_mtr[genes,]
  rowname <- rownames(exp_mtr)
  colname <- colnames(exp_mtr)
  exp_mtr <- apply(exp_mtr, 2, as.numeric)
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
#' @description By calculating the Gini index of different genes, you can get an overview of the classification efficiency of different genes.
#' @param SE a SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @importFrom stats setNames
#' @export

Gini_gene <- function(SE){
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
  ## selected that belongs to category
  ii <- sum(vec[index] == category)
  ## all that belongs to category
  tt <- sum(vec == category)

  1 - (ii/tt)^2 - (1 - ii/tt)^2
}

#' @title differential gene
#' @description Return differential expression gene between Responder and Non-Responder
#' @param SE a SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
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
#' @export

extract_label <-function(datasetNames){
  for (name in datasetNames) {
    if(!exists('inteVector', envir = current_env())){
      get(name) %$% .@colData$response_NR ->response

      inteVector <- response
      if(length(datasetNames > 1))
        next
    }
    get(name) %$% .@colData$response_NR ->response

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
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
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
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param isList whether SE is list
#' @importFrom SummarizedExperiment assay
#' @importFrom magrittr %>%
#' @export

bind_meta <- function(SE,isList){
  if (!isList){
    meta <- as.data.frame(SE@colData)
  } else if (isList){
    meta <- as.data.frame(SE[[1]]@colData)
    meta$batch <- rep('batch1',nrow(meta))
    for (i in 2:length(SE)) {
      paste0('batch', i) %>% rep(nrow(SE[[i]]@colData)) -> batch
      SE[[i]]@colData %>% as.data.frame() %>% cbind(batch) %>% rbind(meta,.) -> meta
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
  model <- model.matrix(~as.factor(meta$response))
  mtr <- dataPreprocess(mtr, rownames(mtr), turn2HL = FALSE)
  mtr <- sva::ComBat(dat = mtr,batch = as.factor(meta$batch),mod = model)
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
#' @import ggplot2
#' @importFrom rlang .data

plt_style <- function(df){
  diff_theme <- theme(plot.title = element_text(face = "bold",
                                                size = "14", color = "#646464"),
                      axis.title = element_text(face = "bold", size = "12", color = "#646464"),
                      axis.text = element_text(face = "bold", size = "9", color = "#646464"),
                      panel.background = element_rect(fill = "white",color = "#646464", size = 1.3),
                      legend.position = "right",
                      legend.title = element_text(face = "bold", size = "14",
                                                  color = "#646464"), panel.grid.major = element_blank(),
                      legend.text = element_text(face='bold',
                                                 size='8.5',color='#646464'),
                      panel.grid.minor = element_blank(),
                      aspect.ratio = 1)
  df$group <- sub("Non-Responder","N",df$group)
  df$group <- sub("Responder","R",df$group)
  df$group <- sub("Post-Therapy","Post",df$group)
  df$group <- sub("Pre-Therapy","Pre",df$group)
  mycolor <- c("#5f96e8","#ee822f")
  names(mycolor) <- unique(df$group)
  ggplot(df, aes(x = .data$group,
                 y = .data$Score,
                 color = .data$group)) +
    scale_color_manual(values = mycolor) +
    geom_boxplot(lwd=1.2) + geom_jitter(aes(fill = .data$group),
                                      width = 0.2, size = 1.5) + diff_theme
}


#' @title Prepare data for plot
#' @description Preparing data for ploting.
#' @param gene is the Gene or Gene set you are interested in.
#' @param SE SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param method the method for calculating gene set scores. Can be NULL if the length of parameter gene is 1.
#' @param type the type of information
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%

plt_Preprocess <- function(gene, SE, method, type){
  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)
  meta <- bind_meta(SE, isList)

  Sc <- Core(exp_mtr, gene, method)
  Score <- log2(Sc + 1)

  if(type == 'R vs NR')
    group <- as.vector(meta$response_NR)
  if(type == 'T vs UT')
    group <- as.vector(meta$Treatment)

  df <- data.frame(group,Score)
  idx <- response_filter(df$group)
  if(length(idx) != 0)
    df <- df[-idx,]

  if(type == 'R vs NR')
    df$group %<>% sub('R','Responder',.) %>% sub('N','Non-Responder',.)
  if(type == 'T vs UT')
    df$group %<>% sub('PRE','Pre-Therapy',.) %>% sub('ON|EDT','Post-Therapy',.)

  return(df)
}


#' @title count geneset score by different method
#' @description wait to write
#' @param exp_mtr an expression matrix.
#' @param geneSet The geneSet which you wanted.
#' @param method the method for calculating gene set scores. Can be NULL if the length of parameter gene is 1.
#' @import survival

Core <- function(exp_mtr, geneSet, method){
  if(is.null(method)){
    return(exp_mtr[geneSet,])
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
    warning(str, "does not exist in expression matrix.")
  }

  exp_mtr <- stats::na.omit(t(apply(exp_mtr, 1, function(x){
    if(all(is.na(x))||all(x == 0)){
      return(rep(NA,length(x)))
    }else{
      return(x)
    }
  })))

  exist_genes <- intersect(rownames(exp_mtr), geneSet)

  Score <-
    switch(method,
           Average_mean = apply(exp_mtr[exist_genes,], 2, mean),
           GSVA = GSVA::gsvaParam(exp_mtr, geneSets=list(exist_genes)) %>% GSVA::gsva(),
           Weighted_mean = weight_mean_signature(exp_mtr, geneSet0[names(geneSet0) %in% exist_genes]))

  return(as.vector(Score))
}

globalVariables(".")
