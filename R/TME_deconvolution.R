#' @title Cibersort functions
#' @description Cibersort functions which perform deconvolution to bulk RNA-seq data. And return the a list which first element is cell fraction and second is a box plot.
#' @param SE an SummarizedExperiment object contains the bulk RNA-seq dataset that you want to use for deconvolution and obtaining its cell fraction.
#' @param method the TME analysis method you want to apply for.
#' @param ... other argument
#' @export

deconv_TME <- function(SE, method, ...){
  switch(method,
         TIMER = TIMER(SE, ...),
         CIBERSORT = CIBERSORT(SE=SE, ...),
         MCPCounter = MCPCounter(SE, ...),
         xCell = xCell(SE, ...),
         IPS = IPS(SE, ...),
         epic = epic(SE, ...),
         ESTIMATE = ESTIMATE(SE, ...),
         ABIS = ABIS(SE, ...),
         ConsensusTME = ConsensusTME(SE, ...),
         quanTIseq = quanTIseq(SE),
         default = stop("The value of parameter 'method' is not available."))
}


#' @title TIMER deconvolution
#' @description use TIMER to predict TME
#' @param SE an SummarizedExperiment object contains the bulk RNA-seq dataset that you want to use for deconvolution and obtaining its cell fraction.
#' @param type the cancer type of data.
#' @export

TIMER <- function(SE,type="SKCM"){
  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)
  TIMER.Immune <- readRDS(system.file("extdata", "TIMER.Immune.rds", package = "tigeR", mustWork = TRUE))

  co_genes <- intersect(rownames(exp_mtr),rownames(TIMER.Immune[[1]]))
  pre_rmBE <- cbind(exp_mtr[co_genes,],TIMER.Immune[[1]][co_genes,])
  batch <- as.factor(c(rep("Tumor",ncol(exp_mtr)),rep("Immune",ncol(TIMER.Immune[[1]]))))
  post_rmBE <- sva::ComBat(pre_rmBE, batch)

  tumor_exp <- post_rmBE[,seq_along(exp_mtr[1,])]
  immune_exp <- post_rmBE[,(ncol(exp_mtr) + 1):ncol(post_rmBE)]

  feature_matrix <- as.data.frame(
    lapply(TIMER.Immune[[2]], function(x){
      apply(immune_exp[,x],1,median)
    }))

  g <- unique(as.vector(apply(feature_matrix, 2, function(x){
    rownames(feature_matrix)[tail(order(x),ceiling(length(x)/100))]
  })))
  feature_matrix <- feature_matrix[!rownames(feature_matrix) %in% g,]

  TIMER.Markers <- readRDS(system.file("extdata", "TIMER.Markers.rds", package = "tigeR", mustWork = TRUE))

  selected_genes <- intersect(TIMER.Markers[[type]],rownames(feature_matrix))

  cancer.expression <- tumor_exp[selected_genes,]
  feature.expression <- feature_matrix[selected_genes,]

  fraction_matrix <-
    apply(cancer.expression, 2, function(x) {
      fraction <- stats::lsfit(feature.expression,x,intercept=FALSE)$coefficients

      drop <- c()
      while(any(fraction<0)){
        drop <- c(drop,which.min(fraction))
        fraction <- rep(0,length(fraction))
        fraction[-drop] <- stats::lsfit(feature.expression[,-drop],x,intercept=FALSE)$coefficients
      }
      fraction
    })
  rownames(fraction_matrix) <- colnames(feature_matrix)
  fraction_matrix
}


#' @title Cibersort functions
#' @description Cibersort functions which perform deconvolution to bulk RNA-seq data. And return the a list which first element is cell fraction and second is a box plot.
#' @param SE an SummarizedExperiment object contains the bulk RNA-seq dataset that you want to use for deconvolution and obtaining its cell fraction.
#' @param sig_matrix gene expression matrix from isolated cells.
#' @param perm the number of permutations.
#' @param QN whether perform quantile normalization or not (TRUE/FALSE).
#' @param style the ploting style. c("raw,"elegant","brief")
#' @param group_color the color of Responder and Non_Responder.
#' @importFrom magrittr %>%
#' @importFrom stats wilcox.test
#' @export

CIBERSORT <- function(SE, sig_matrix, perm=0, QN=TRUE, style='elegant', group_color=c("#5f96e8CC", "#ee822fCC")){
  if(missing(sig_matrix)){
    LM22 <- NULL
    data(LM22,package = "tigeR", envir = current_env())
    sig_matrix <- LM22
  }

  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)

  result <- Ciber(sig_matrix,exp_mtr,perm,QN)

  TME_data <- as.data.frame(result[,1:22])
  TME_data$group <- bind_meta(SE, isList)$response_NR
  TME_data$sample <- rownames(TME_data)

  TME_New <- reshape2::melt(TME_data)

  colnames(TME_New) <- c("Group","Sample","Celltype","Composition")

  plot_order <- TME_New[TME_New$Group=="R",] %>%
    dplyr::group_by(.data$Celltype) %>%
    dplyr::summarise(m = stats::median(.data$Composition)) %>%
    dplyr::arrange(dplyr::desc(.data$m)) %>%
    dplyr::pull(.data$Celltype)

  TME_New$Celltype <- factor(TME_New$Celltype,levels=plot_order)
  TME_New<-TME_New[!TME_New$Group=='UNK',]
  rs <- c()

  warning_status <- 0
  for (i in levels(TME_New$Celltype)) {
    m <- TME_New[TME_New$Celltype==i,]
    R_S <- m[m$Group == 'R', 4]
    N_S <- m[m$Group == "N", 4]
    if(any(R_S%in%N_S))
      warning_status <- 1
    rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
  }
  if(warning_status)
    message("There are identical relative abundance values in groups R and N for the '",i,"'. The p value of the Wilcoxon signed-rank test may not be precise due to ties in the data.")

  #selected_cells <- levels(TME_New$Celltype)[which(rs < 0.05)]
  selected_cells <- levels(TME_New$Celltype)

  if(length(selected_cells) < 5){
    names(rs) <- seq_along(rs)
    selected_cells <- levels(TME_New$Celltype)[as.numeric(names(sort(rs)[1:5]))]
  }

  ciber_theme <- ggplot2::theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                                axis.title = element_text(size = 10,color ="black"),
                                axis.text = element_text(size= 10,color = "black"),
                                axis.text.x = element_text(angle = 45, hjust = 1 ),
                                legend.position = "top",
                                legend.text = element_text(size= 12),
                                legend.title= element_text(size= 12))
  if(style == 'raw'){
    box_TME <-
      ggplot2::ggplot(TME_New[TME_New$Celltype%in%selected_cells,], aes(x = .data$Celltype, y = .data$Composition)) +
      ggplot2::labs(y="Cell composition",x= NULL,title = "TME Cell composition") +
      ggplot2::geom_boxplot(aes(fill = .data$Group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0) +
      ggplot2::scale_fill_manual(values = group_color) +
      ggplot2::theme_classic() + ciber_theme
    y_max <- max(ggplot_build(box_TME)$data[[1]]$ymax)
    box_TME <-
      box_TME +
      ggpubr::stat_compare_means(ggplot2::aes(group = .data$Group),
                                 label = "p.signif",
                                 method = "wilcox.test",
                                 hide.ns = TRUE,
                                 label.y = y_max*1.1) +
      coord_cartesian(ylim = c(0, y_max*1.1))
  }
  if(style == 'elegant'){
    box_TME <- ggplot(TME_New,aes(x=.data$Celltype,y=.data$Composition,fill=.data$Group)) +
      stat_boxplot(data=TME_New,
                   geom = "errorbar",width = 1, color = "black",linetype = "solid",
                   position = position_dodge(0.8),linewidth = 0.7) +
      stat_boxplot(geom = "boxplot", color = "black",linetype = "solid",
                   position = position_dodge(0.8),linewidth = 0.7,
                   width = 0.8, outlier.shape= 19) +
      ggpubr::theme_classic2() + ciber_theme +
      scale_fill_manual(values = group_color)
    y_max <- max(ggplot_build(box_TME)$data[[1]]$ymax)
    box_TME <- box_TME +
      ggpubr::stat_compare_means(ggplot2::aes(group = .data$Group),
                                 label = "p.signif", method = "wilcox.test", hide.ns = FALSE,
                                 label.y = y_max * 1.1,size = 3) + coord_cartesian(ylim = c(0, y_max * 1.1))
  }

  list(t(result), box_TME)
}


#' @title Cibersort functions
#' @description Cibersort functions which perform deconvolution to bulk RNA-seq data. And return the cell fraction.
#' @param sig_matrix file path to gene expression from isolated cells
#' @param mix_matrix heterogenous mixed expression
#' @param perm Number of permutations
#' @param QN Perform quantile normalization or not (TRUE/FALSE)

Ciber <- function(sig_matrix, mix_matrix, perm=0, QN=TRUE){
  #read in data
  X <- sig_matrix
  Y <- mix_matrix

  #order
  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]
  Y <- na.omit(Y)

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
  select_genes <- S4Vectors::intersect(rownames(X),rownames(Y))
  X <- X[select_genes,]
  Y <- Y[select_genes,]

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


#' @title core algorithm
#' @description core algorithm of CIBERSORT
#' @param X cell-specific gene expression
#' @param Y mixed expression per sample
#' @importFrom stats cor

CoreAlg <- function(X, Y){
  #try different values of nu
  svn_itor <- 3

  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-e1071::svm(X,Y,type="nu-regression",kernel="linear",nu=nus,scale=FALSE)
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
    weights <- t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w <- weights/sum(weights)
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
#' @importFrom stats sd

doPerm <- function(perm, X, Y){
  itor <- 1
  Ylist <- as.list(Y)
  dist <- matrix()

  while(itor <= perm){
    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])

    #standardize mixture
    yr <- (yr - mean(yr)) / stats::sd(yr)

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


#' @title MCPCounter deconvolution
#' @description use MCPCounter to predict TME
#' @param SE an SummarizedExperiment object contains the bulk RNA-seq dataset that you want to use for deconvolution and obtaining its cell fraction.
#' @param featuresType type of identifiers for expression features. Defaults to "affy133P2_probesets" for Affymetrix Human Genome 133 Plus 2.0 probesets. Other options are "HUGO_symbols" (Official gene symbols), "ENTREZ_ID" (Entrez Gene ID) or "ENSEMBL_ID" (ENSEMBL Gene ID)
#' @param ... other parameter
#' @export

MCPCounter <- function(SE, featuresType = "HUGO_symbols", ...) {
  featuresType <- match.arg(featuresType,c("affy133P2_probesets","HUGO_symbols","ENTREZ_ID","ENSEMBL_ID"))
  isList <- is.list(SE)
  gene_expression_matrix <- bind_mtr(SE, isList)

  arguments <- rlang::dots_list(gene_expression_matrix, featuresType = featuresType, ..., .homonyms = "last")
  call <- rlang::call2(MCPcounter.estimate, !!!arguments)
  eval(call)
}


#' @title ESTIMATE deconvolution
#' @description use ESTIMATE to predict TME
#' @param expression matrix or data.frame with features in rows and samples in columns
#' @param featuresType type of identifiers for expression features. Defaults to "affy133P2_probesets" for Affymetrix Human Genome 133 Plus 2.0 probesets. Other options are "HUGO_symbols" (Official gene symbols), "ENTREZ_ID" (Entrez Gene ID) or "ENSEMBL_ID" (ENSEMBL Gene ID)
#' @param probesets probs
#' @param genes genes
#' @export

MCPcounter.estimate<-function(expression,featuresType,probesets,genes){
  if(missing(probesets))
    probesets <- read.table(system.file("extdata", "probesets.txt", package = "tigeR", mustWork = TRUE),
                            sep="\t",stringsAsFactors=FALSE,colClasses="character")
  if(missing(genes))
    genes <- read.table(system.file("extdata", "genes.txt", package = "tigeR", mustWork = TRUE),
                        sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)

  if(featuresType == "affy133P2_probesets"){
    features <- probesets
    markers.names <- unique(features[, 2])
    features <- split(features[,1],features[,2])
    features <- lapply(features,intersect,x=rownames(expression))
    features <- features[sapply(features,function(x)length(x)>0)]
    missing.populations <- setdiff(markers.names,names(features))
    features <- features[intersect(markers.names,names(features))]
  } else {
    markersG <- genes
  }

  type <- sub("_"," ",featuresType)
  features <- subset(markersG,markersG[,type]%in%rownames(expression))
  markers.names <- unique(features[, "Cell population"])
  features <- split(features[,type],features[,"Cell population"])
  missing.populations <- setdiff(markers.names,names(features))
  features <- features[intersect(markers.names,names(features))]

  if(length(missing.populations)>0){
    warning(paste("Found no markers for population(s):",paste(missing.populations,collapse=", ")))
  }

  t(as.data.frame(do.call(cbind,
                          lapply(features,function(x){
                            apply(expression[intersect(row.names(expression),x),,drop=FALSE],
                                  2,mean,na.rm=TRUE)
                          }))))
}


#' @title xCell deconvolution
#' @description use xCell to predict TME
#' @param SE an SummarizedExperiment object contains the bulk RNA-seq dataset that you want to use for deconvolution and obtaining its cell fraction.
#' @param signatures a GMT object of signatures.
#' @param genes list of genes to use in the analysis.
#' @param spill the Spillover object for adjusting the scores.
#' @param rnaseq if true than use RNAseq spillover and calibration paramters, else use array parameters.
#' @param file.name string for the file name for saving the scores. Default is NULL.
#' @param scale if TRUE, uses scaling to trnasform scores using fit.vals
#' @param alpha a value to override the spillover alpha parameter. Deafult = 0.5
#' @param save.raw TRUE to save a raw
#' @param parallel.sz integer for the number of threads to use. Default is 4.
#' @param parallel.type Type of cluster architecture when using snow. 'SOCK' or 'FORK'. Fork is faster, but is not supported in windows.
#' @param cell.types.use a character list of the cell types to use in the analysis. If NULL runs xCell with all cell types.
#' @export

xCell <- function(SE, signatures=NULL, genes=NULL, spill=NULL, rnaseq=TRUE, file.name = NULL, scale=TRUE,
                  alpha = 0.5, save.raw = FALSE, parallel.sz = 4, parallel.type = 'SOCK',
                  cell.types.use = NULL){
  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)

  autoload("xCell.data","xCell")
  xCell::xCellAnalysis(exp_mtr, signatures, genes, spill, rnaseq, file.name, scale,
                       alpha, save.raw, parallel.sz, parallel.type, cell.types.use)
}


#' @title IPS deconvolution
#' @description use xCell to predict TME
#' @param SE an SummarizedExperiment object contains the bulk RNA-seq dataset that you want to use for deconvolution and obtaining its cell fraction.
#' @param project project
#' @param plot if TRUE return the plot
#' @export

IPS <- function(SE, project=NULL,plot=FALSE){
  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)
  meta <- bind_meta(SE, isList)

  if(plot){
    my_palette <- grDevices::colorRampPalette(c("blue", "white", "red"))(n = 1000)
    my_palette2 <- grDevices::colorRampPalette(c("black", "white"))(n = 1000)
  }

  IPSG <- readRDS(system.file("extdata", "IPSG.rds", package = "tigeR", mustWork = TRUE))
  IPSG <- IPSG[IPSG$GENE %in% rownames(exp_mtr),]
  unique_ips_genes <- as.vector(unique(IPSG$NAME))

  IPS<-NULL
  MHC<-NULL
  CP<-NULL
  EC<-NULL
  SC<-NULL
  AZ<-NULL

  GVEC <- rownames(exp_mtr)
  VEC <- as.vector(IPSG$GENE)
  ind <- which(is.na(match(VEC,GVEC)))
  MISSING_GENES <- VEC[ind]
  dat <- IPSG[ind,]
  if (length(MISSING_GENES)>0) {
    cat("differently named or missing genes: ",MISSING_GENES,"\n")
  }

  for (i in seq_along(exp_mtr[1,])) {
    GE <- exp_mtr[,i]
    mGE <- mean(GE,na.rm=TRUE)
    sGE <- sd(GE,na.rm=TRUE)
    Z1 <- (exp_mtr[as.vector(IPSG$GENE),i] - mGE)/sGE
    W1 <- IPSG$WEIGHT
    WEIGHT <- NULL
    MIG <- NULL
    k <- 1
    for (gen in unique_ips_genes) {
      MIG[k] <- mean(Z1[which(as.vector(IPSG$NAME)==gen)],na.rm=TRUE)
      WEIGHT[k] <- mean(W1[which (as.vector(IPSG$NAME)==gen)],na.rm=TRUE)
      k <- k+1
    }
    WG <- MIG*WEIGHT
    MHC[i] <- mean(WG[1:10],na.rm=TRUE)
    CP[i] <- mean(WG[11:20],na.rm=TRUE)
    EC[i] <- mean(WG[21:24],na.rm=TRUE)
    SC[i] <- mean(WG[25:26],na.rm=TRUE)
    AZ[i] <- sum(MHC[i],CP[i],EC[i],SC[i],na.rm = TRUE)
    IPS[i] <- ipsmap(AZ[i])

    if (plot) {
      data_a <- data.frame (start = c(0,2.5,5,7.5,10,15,seq(20,39),0,10,20,30),
                            end = c(2.5,5,7.5,10,15,seq(20,40),10,20,30,40), y1=c(rep(2.6,26),rep(0.4,4)),
                            y2=c(rep(5.6,26),rep(2.2,4)),z=c(MIG[c(21:26,11:20,1:10)],EC[i],SC[i],CP[i],MHC[i]),
                            vcol=c(unlist(lapply(MIG[c(21:26,11:20,1:10)],mapcolors,my_palette=my_palette)),
                                   unlist(lapply(c(EC[i],SC[i],CP[i],MHC[i]),mapbw,my_palette2=my_palette2))),
                            label = c(unique_ips_genes[c(21:26,11:20,1:10)],"EC","SC","CP","MHC"))
      data_a$label <- factor(data_a$label, levels=unique(data_a$label))
      plot_a1<-ggplot() + geom_rect(data=data_a,
                                    mapping=aes(xmin=.data$start, xmax=.data$end,
                                                ymin=.data$y1, ymax=.data$y2, fill=.data$label),
                                    size=0.5,color="black", alpha=1) +
        coord_polar() + scale_y_continuous(limits = c(0, 6)) +
        scale_fill_manual(values =as.vector(data_a$vcol),guide=FALSE) +
        theme_bw() + theme(panel.margin = unit(0, 'mm'), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),panel.border = element_blank(),
                           panel.background = element_blank(), axis.line = element_line(colour = "white"),
                           axis.text=element_blank(), axis.ticks= element_blank()) +
        geom_text(aes(x=5, y=1.3, label="EC"), size=4) +
        geom_text(aes(x=15, y=1.3, label="SC"), size=4) +
        geom_text(aes(x=25, y=1.3, label="CP"), size=4) +
        geom_text(aes(x=35, y=1.3, label="MHC"), size=4)

      plot_a2<-plot_a1+geom_text(aes(x=1.25, y=4.1, label="+ Act CD4"), angle=78.75, size=4)+geom_text(aes(x=3.75, y=4.1, label="+ Act CD8"),angle=56.25, size=4)+geom_text(aes(x=6.25, y=4.1, label="+ Tem CD4"), angle=33.75,size=4)+geom_text(aes(x=8.75, y=4.1, label="+ Tem CD8"), angle=11.25,size=4)+geom_text(aes(x=17.5, y=4.1, label="- MDSC"), angle=-67.5,size=4)+geom_text(aes(x=12.5, y=4.1, label="- Treg"), angle=-22.5,size=4)
      plot_a3<-plot_a2+geom_text(aes(x=20.5, y=4.1, label="PD-1 -"), angle=85.5, size=4)+geom_text(aes(x=21.5, y=4.1, label="CTLA4 -"), angle=76.5, size=4)+geom_text(aes(x=22.5, y=4.1, label="LAG3 -"), angle=67.5, size=4)+geom_text(aes(x=23.5, y=4.1, label="TIGIT -"), angle=58.5, size=4)+geom_text(aes(x=24.5, y=4.1, label="TIM3 -"), angle=49.5, size=4)+geom_text(aes(x=25.5, y=4.1, label="PD-L1 -"), angle=40.5, size=4)+geom_text(aes(x=26.5, y=4.1, label="PD-L2 -"), angle=31.5, size=4)+geom_text(aes(x=27.5, y=4.1, label="CD27 +"), angle=22.5, size=4)+geom_text(aes(x=28.5, y=4.1, label="ICOS +"), angle=13.5, size=4)+geom_text(aes(x=29.5, y=4.1, label="IDO1 -"), angle=4.5, size=4)
      plot_a4<-plot_a3+geom_text(aes(x=30.5, y=4.1, label="B2M +"), angle=-4.5, size=4)+geom_text(aes(x=31.5, y=4.1, label="TAP1 +"), angle=-13.5, size=4)+geom_text(aes(x=32.5, y=4.1, label="TAP2 +"), angle=-22.5, size=4)+geom_text(aes(x=33.5, y=4.1, label="HLA-A +"), angle=-31.5, size=4)+geom_text(aes(x=34.5, y=4.1, label="HLA-B +"), angle=-40.5, size=4)+geom_text(aes(x=35.5, y=4.1, label="HLA-C +"), angle=-49.5, size=4)+geom_text(aes(x=36.5, y=4.1, label="HLA-DPA1 +"), angle=-58.5, size=4)+geom_text(aes(x=37.5, y=4.1, label="HLA-DPB1 +"), angle=-67.5, size=4)+geom_text(aes(x=38.5, y=4.1, label="HLA-E +"), angle=-76.5, size=4)+geom_text(aes(x=39.5, y=4.1, label="HLA-F +"), angle=-85.5, size=4)
      plot_a5<-plot_a4+geom_text(aes(x=0, y=6, label=paste("Immunophenoscore: ",IPS[i],sep="")), angle=0,size=6,vjust=-0.5)+ theme(axis.title=element_blank())
      plot_a <-plot_a5 + theme(plot.margin=unit(c(0,0,0,0),"mm")) + geom_text(vjust=1.15,hjust=0,aes(x=25.5, y=6,label="\n\n\n\n   MHC: Antigen Processing                                 EC: Effector Cells\n   CP: Checkpoints | Immunomodulators              SC: Suppressor Cells\n\n", hjust = 0), size=4)

      ## Legend sample-wise (averaged) z-scores
      data_b <- data.frame (start = rep(0,23), end = rep(0.7,23), y1=seq(0,22,by=1), y2=seq(1,23,by=1),z=seq(-3,3,by=6/22),vcol=c(unlist(lapply(seq(-3,3,by=6/22),mapcolors,my_palette=my_palette))), label = LETTERS[1:23])
      data_b_ticks <- data.frame(x = rep(1.2, 7), value = seq(-3,3, by=1), y = seq(0,6, by=1)*(22/6) +0.5)
      legendtheme <- theme(plot.margin = unit(c(2,0,2,0),"inch"),
                           panel.margin = unit(0,"null"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank(),
                           axis.line = element_line(colour = "white"),
                           axis.text=element_blank(),
                           axis.ticks= element_blank(),
                           axis.title.x=element_blank())
      plot_b<-ggplot(hjust=0) +
        geom_rect(data=data_b, mapping=aes(xmin=.data$start, xmax=.data$end, ymin=.data$y1, ymax=.data$y2, fill=.data$label), size=0.5,color="black", alpha=1) +
        scale_x_continuous(limits = c(0, 1.5),expand = c(0,0)) +
        scale_fill_manual(values =as.vector(data_b$vcol),guide=FALSE) +
        geom_text(data=data_b_ticks, aes(x=.data$x, y=.data$y, label=.data$value),hjust="inward", size=4) +
        theme_bw() + legendtheme + ylab("Sample-wise (averaged) z-score")

      ## Legend weighted z-scores
      data_c <- data.frame (start = rep(0,23), end = rep(0.7,23), y1=seq(0,22,by=1),
                            y2=seq(1,23,by=1),z=seq(-2,2,by=4/22),
                            vcol=c(unlist(lapply(seq(-2,2,by=4/22),mapbw,my_palette2=my_palette2))),
                            label = LETTERS[1:23])
      data_c_ticks <- data.frame(x = rep(1.2, 5), value = seq(-2,2, by=1), y = seq(0,4, by=1)*(22/4) +0.5)

      plot_c<-ggplot() + geom_rect(data=data_c, mapping=aes(xmin=.data$start, xmax=.data$end,
                                                            ymin=.data$y1, ymax=.data$y2, fill=.data$label),
                                   size=0.5,color="black", alpha=1) +
        scale_x_continuous(limits = c(0, 1.5),expand = c(0,0)) +
        scale_fill_manual(values =as.vector(data_c$vcol),guide=FALSE) +
        geom_text(data=data_c_ticks, aes(x=.data$x, y=.data$y, label=.data$value),hjust="inward", size=4) +
        theme_bw() + legendtheme + ylab("Weighted z-score")
      final_plot <-
        gridExtra::grid.arrange(plot_a,plot_b,plot_c, ncol=3, widths=c(0.8,0.1,0.1))
    }
  }

  res<-data.frame(ID=colnames(exp_mtr),MHC=MHC,EC=EC,SC=SC,CP=CP,AZ=AZ,IPS=IPS)

  if(!is.null(project)){
    res$ProjectID <- project
    res <- res[,c(ncol(res),1:ncol(res)-1)]
  }

  res<-tibble::column_to_rownames(res,var = "ID")
  if(plot)
    return(list(res,final_plot))
  return(t(res))
}

#' @title dd
#' @param x ff
#' @export

ipsmap<- function (x) {
  # if(is.na(x)) x<-0
  if (x<=0) {
    ips<-0
  } else {
    if (x>=3) {
      ips<-10
    } else {
      ips<-round(x*10/3, digits=0)
    }
  }
  return(ips)
}


#' @title dd
#' @param x dd
#' @param my_palette description
#' @export

mapcolors<-function (x, my_palette) {
  za<-NULL
  if (x>=3) {
    za=1000
  } else {
    if (x<=-3) {
      za=1
    } else {
      za=round(166.5*x+500.5,digits=0)
    }
  }
  return(my_palette[za])
}

#' @title dd
#' @param x f
#' @param my_palette2 description
#' @export
mapbw <- function (x, my_palette2) {
  za2<-NULL
  if (x>=2) {
    za2=1000
  } else {
    if (x<=-2) {
      za2=1
    } else {
      za2=round(249.75*x+500.5,digits=0)
    }
  }
  return(my_palette2[za2])
}

#' @title epic deconvolution
#' @description use epic to predict TME
#' @param SE an SummarizedExperiment object contains the bulk RNA-seq dataset that you want to use for deconvolution and obtaining its cell fraction.
#' @param ... the arguments
#' @export

epic <- function(SE, ...){
  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)

  t(EPIC::EPIC(bulk = exp_mtr, ...)$cellFractions)
}


#' @title ESTIMATE deconvolution
#' @description use ESTIMATE to predict TME
#' @param SE an SummarizedExperiment object contains the bulk RNA-seq dataset that you want to use for deconvolution and obtaining its cell fraction.
#' @export

ESTIMATE <- function(SE) {
  isList <- is.list(SE)
  gene_expression_matrix <- bind_mtr(SE, isList)

  Common_genes <- readRDS(system.file("extdata", "Common_genes.rds",
                                      package = "tigeR", mustWork = TRUE))
  Signature_genesets <- readRDS(system.file("extdata", "Signature_genesets.rds",
                                            package = "tigeR", mustWork = TRUE))

  merged.df <- merge(Common_genes, gene_expression_matrix, by.x = "GeneSymbol", by.y = "row.names")
  rownames(merged.df) <- merged.df$GeneSymbol
  filtered_matrix <- merged.df[, -1:-(ncol(Common_genes))]
  print(sprintf(
    "Merged dataset includes %d genes (%d mismatched).",
    nrow(filtered_matrix), nrow(Common_genes) - nrow(filtered_matrix)
  ))


  m <- filtered_matrix
  gene.names <- rownames(m)
  sample.names <- colnames(m)

  Ns <- length(m[1, ])
  Ng <- length(m[, 1])


  for (j in 1:Ns) {
    m[, j] <- rank(m[, j], ties.method = "average")
  }
  m <- 10000 * m / Ng
  gs <- as.matrix(Signature_genesets[, -1], dimnames = NULL)
  N.gs <- 2
  gs.names <- row.names(Signature_genesets)
  score.matrix <- matrix(0, nrow = N.gs, ncol = Ns)
  for (gs.i in 1:N.gs) {
    gene.set <- gs[gs.i, ]
    gene.overlap <- intersect(gene.set, gene.names)
    print(paste(
      gs.i, "gene set:", gs.names[gs.i], " overlap=",
      length(gene.overlap)
    ))
    if (length(gene.overlap) == 0) {
      score.matrix[gs.i, ] <- rep(NA, Ns)
      next
    } else {
      ES.vector <- vector(length = Ns)
      for (S.index in 1:Ns) {
        gene.list <- order(m[, S.index], decreasing = TRUE)
        gene.set2 <- match(gene.overlap, gene.names)
        correl.vector <- m[gene.list, S.index]
        TAG <- sign(match(gene.list, gene.set2, nomatch = 0))
        no.TAG <- 1 - TAG
        N <- length(gene.list)
        Nh <- length(gene.set2)
        Nm <- N - Nh
        correl.vector <- abs(correl.vector)^0.25
        sum.correl <- sum(correl.vector[TAG == 1])
        P0 <- no.TAG / Nm
        F0 <- cumsum(P0)
        Pn <- TAG * correl.vector / sum.correl
        Fn <- cumsum(Pn)
        RES <- Fn - F0
        max.ES <- max(RES)
        min.ES <- min(RES)
        if (max.ES > -min.ES) {
          arg.ES <- which.max(RES)
        } else {
          arg.ES <- which.min(RES)
        }
        ES <- sum(RES)
        EnrichmentScore <- list(
          ES = ES, arg.ES = arg.ES,
          RES = RES, indicator = TAG
        )
        ES.vector[S.index] <- EnrichmentScore$ES
      }
      score.matrix[gs.i, ] <- ES.vector
    }
  }
  score.data <- data.frame(score.matrix)
  names(score.data) <- sample.names
  row.names(score.data) <- gs.names
  estimate.score <- apply(score.data, 2, sum)

  convert_row_estimate_score_to_tumor_purity <- function(x) {
    stopifnot(is.numeric(x))
    cos(0.6049872018 + 0.0001467884 * x)
  }
  est.new <- NULL
  for (i in seq_along(estimate.score)) {
    est_i <- convert_row_estimate_score_to_tumor_purity(estimate.score[i])
    est.new <- rbind(est.new, est_i)
    if (est_i >= 0) {
      next
    } else {
      message(paste(sample.names[i], ": out of bounds",
                    sep = ""
      ))
    }
  }

  colnames(est.new) <- c("TumorPurity")
  estimate.t1 <- cbind(estimate.score, est.new)
  x.bad.tumor.purities <- estimate.t1[, "TumorPurity"] < 0
  estimate.t1[x.bad.tumor.purities, "TumorPurity"] <- NA
  score.data <- rbind(score.data, t(estimate.t1))
  rownames(score.data) <- c(
    "StromalScore", "ImmuneScore",
    "ESTIMATEScore", "TumorPurity"
  )

  as.matrix(score.data)
}


#' @title abis deconvolution
#' @description use abis to predict TME
#' @param SE an SummarizedExperiment object contains the bulk RNA-seq dataset that you want to use for deconvolution and obtaining its cell fraction.
#' @param sig_matrix gene expression matrix from isolated cells.
#' @importFrom stats coef
#' @export

ABIS <- function(SE, sig_matrix) {
  if(missing(sig_matrix)){
    LM22 <- NULL
    data(LM22,package = "tigeR", envir = current_env())
    sig_matrix <- LM22
  }

  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)

  genes <- intersect(rownames(exp_mtr), rownames(sig_matrix))
  Dec <- (apply(exp_mtr[genes, , drop = FALSE],2,
                function(x) coef(MASS::rlm(as.matrix(sig_matrix[genes, ]), x, maxit = 100)))) * 100
  Dec <- signif(Dec, 3)
}


#' @title ConsensusTME deconvolution
#' @description use ConsensusTME to predict TME
#' @param SE an SummarizedExperiment object contains the bulk RNA-seq dataset that you want to use for deconvolution and obtaining its cell fraction.
#' @param cancer string passed to indicate which TCGA cancer type samples are most similar to. N.B samples of different cancer types should be run seperately. Available cancer types: "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP","LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM".
#' @param ... the arguments
#' @export

ConsensusTME <- function(SE, cancer="SKCM", ...){
  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)

  ConsensusTME::consensusTMEAnalysis(exp_mtr, cancer, statMethod = "ssgsea", ...)
}


#' @title quanTIseq deconvolution
#' @description use quanTIseq to predict TME
#' @param SE a SummarizedExperiment object contains the bulk RNA-seq dataset that you want to use for deconvolution and obtaining its cell fraction.
#' @param tumor whether the input is tumor sample
#' @param arrays whether you data is array data
#' @param scale_mrna whether perform scaling
#' @param ... other parameter
#' @importFrom magrittr %>%
#' @export

quanTIseq <- function(SE, tumor=TRUE, arrays=FALSE, scale_mrna=FALSE, ...){
  isList <- is.list(SE)
  gene_expression_matrix <- bind_mtr(SE, isList)

  res <- quantiseqr::run_quantiseq(
    expression_data = gene_expression_matrix,
    is_arraydata = arrays,
    is_tumordata = tumor,
    scale_mRNA = scale_mrna,
    ...
  )

  t(as.matrix(res[,-1]))
}
