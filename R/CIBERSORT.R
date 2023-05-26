#' @title core algorithm
#' @description core algorithm of CIBERSORT
#' @param X cell-specific gene expression
#' @param Y mixed expression per sample

CoreAlg <- function(X, Y){
  #try different values of nu
  svn_itor <- 3

  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-e1071::svm(X,Y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }

  if(Sys.info()['sysname'] == 'Windows')
    out <- parallel::mclapply(1:svn_itor, res, mc.cores=1)
  else
    out <- parallel::mclapply(1:svn_itor, res, mc.cores=svn_itor)

  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)

  #do cibersort
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - Y)^2)))
    corrv[t] <- cor(k, Y)
    t <- t + 1
  }

  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]

  #get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  w <- (q/sum(q))

  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]

  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
}

#' @title do permutations
#' @description do permutations
#' @param perm Number of permutations
#' @param X cell-specific gene expression
#' @param Y mixed expression per sample

doPerm <- function(perm, X, Y){
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()

  while(itor <= perm){
    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])

    #standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)

    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr)

    mix_r <- result$mix_r

    #store correlation
    if(itor == 1) {dist <- mix_r}
    else {dist <- rbind(dist, mix_r)}

    itor <- itor + 1
  }
  newList <- list("dist" = dist)
}

#' @title Cibersort functions
#' @description Cibersort functions which perform deconvolution to bulk RNA-seq data. And return the cell fraction.
#' @param sig_matrix file path to gene expression from isolated cells
#' @param mix_matrix heterogenous mixed expression
#' @param perm Number of permutations
#' @param QN Perform quantile normalization or not (TRUE/FALSE)

Ciber <- function(sig_matrix, mix_matrix, perm=0, QN=TRUE){
  #read in data
  X <- data.matrix(sig_matrix)
  Y <- data.matrix(mix_matrix)

  #order
  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]

  P <- perm #number of permutations

  #anti-log if max < 50 in mixture file
  if(max(Y) < 50) {Y <- 2^Y}

  #quantile normalization of mixture file
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- preprocessCore::normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }

  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]

  #standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))

  #empirical null distribution of correlation coefficients
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y)$dist)}

  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")

  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999

  #iterate through mixtures
  while(itor <= mixtures){

    y <- Y[,itor]

    #standardize mixture
    y <- (y - mean(y)) / sd(y)

    #run SVR core algorithm
    result <- CoreAlg(X, y)

    #get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse

    #calculate p-value
    if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}

    #print output
    out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    if(itor == 1) {output <- out}
    else {output <- rbind(output, out)}

    itor <- itor + 1

  }

  #return matrix object containing all results
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")
  obj
}


#' @title Cibersort functions
#' @description Cibersort functions which perform deconvolution to bulk RNA-seq data. And return the a list which first element is cell fraction and second is a box plot.
#' @param sig_matrix file path to gene expression from isolated cells
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param perm Number of permutations
#' @param QN Perform quantile normalization or not (TRUE/FALSE)
#' @importFrom magrittr %>%
#' @importFrom reshape2 melt
#' @importFrom dplyr desc
#' @importFrom dplyr pull
#' @importFrom dplyr arrange
#' @importFrom dplyr summarise
#' @importFrom ggpubr group_by
#' @import ggplot2
#' @export

CIBERSORT <- function(sig_matrix, SE, perm=0, QN=TRUE){
  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)

  result <- Ciber(sig_matrix,exp_mtr,perm,QN)

  TME_data <- as.data.frame(result[,1:22])
  TME_data$group <- bind_meta(SE, isList)$response_NR
  TME_data$sample <- rownames(TME_data)

  TME_New <- melt(TME_data)

  colnames(TME_New) <- c("Group","Sample","Celltype","Composition")

  plot_order <- TME_New[TME_New$Group=="R",] %>%
    group_by(.data$Celltype) %>%
    summarise(m = median(.data$Composition)) %>%
    arrange(desc(.data$m)) %>%
    pull(.data$Celltype)

  TME_New$Celltype = factor(TME_New$Celltype,levels = plot_order)

  ciber_theme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 10,color ="black"),
                   axis.text = element_text(size= 10,color = "black"),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12))

  box_TME <- ggplot(TME_New, aes(x = .data$Celltype, y = .data$Composition))+
    labs(y="Cell composition",x= NULL,title = "TME Cell composition")+
    geom_boxplot(aes(fill = .data$Group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+
    scale_fill_manual(values = c("#99CCFF", "#CCCC00"))+
    theme_classic() + ciber_theme +
    stat_compare_means(aes(group =  .data$Group),
                       label = "p.signif",
                       method = "wilcox.test",
                       hide.ns = T)
  list(result, box_TME)
}
