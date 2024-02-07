#' @title Cibersort functions
#' @description Cibersort functions which perform deconvolution to bulk RNA-seq data. And return the a list which first element is cell fraction and second is a box plot.
#' @param sig_matrix gene expression matrix from isolated cells.
#' @param SE an SummarizedExperiment object contains the bulk RNA-seq dataset that you want to use for deconvolution and obtaining its cell fraction.
#' @param perm the number of permutations.
#' @param QN whether perform quantile normalization or not (TRUE/FALSE).
#' @importFrom magrittr %>%
#' @importFrom stats wilcox.test
#' @export

CIBERSORT <- function(sig_matrix, SE, perm=0, QN=TRUE){
  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)

  result <- Ciber(sig_matrix,exp_mtr,perm,QN)

  TME_data <- as.data.frame(result[,1:22])
  idx <- which(apply(TME_data, 2, mean) > 0.005)
  TME_data <- TME_data[,idx]
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
  for (i in levels(TME_New$Celltype)) {
    m <- TME_New[TME_New$Celltype==i,]
    rs <- c(rs,stats::wilcox.test(m[m$Group=='R',4],m[m$Group=='N',4])$p.value)
  }
  selected_cells <- levels(TME_New$Celltype)[which(rs < 0.05)]
  ciber_theme <- ggplot2::theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                                axis.title = element_text(size = 10,color ="black"),
                                axis.text = element_text(size= 10,color = "black"),
                                axis.text.x = element_text(angle = 45, hjust = 1 ),
                                legend.position = "top",
                                legend.text = element_text(size= 12),
                                legend.title= element_text(size= 12))
  box_TME <-
    ggplot2::ggplot(TME_New[TME_New$Celltype%in%selected_cells,], aes(x = .data$Celltype, y = .data$Composition)) +
    ggplot2::labs(y="Cell composition",x= NULL,title = "TME Cell composition") +
    ggplot2::geom_boxplot(aes(fill = .data$Group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0) +
    ggplot2::scale_fill_manual(values = c("#5f96e8CC", "#ee822fCC")) +
    ggplot2::theme_classic() + ciber_theme
  y_max <- max(ggplot_build(box_TME)$data[[1]]$ymax)
  box_TME <-
    box_TME +
    ggpubr::stat_compare_means(ggplot2::aes(group = .data$Group),
                               label = "p.signif",
                               method = "wilcox.test",
                               hide.ns = T,
                               label.y.npc = y_max*1.4) +
    coord_cartesian(ylim = c(0, y_max*1.1))

  list(result, box_TME)
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


#' @title xCell deconvolution
#' @description use xCell to predict TME
#' @param SE an SummarizedExperiment object contains the bulk RNA-seq dataset that you want to use for deconvolution and obtaining its cell fraction.
#' @param ... the arguments
#' @export

xCell <- function(SE, ...){
  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)

  xCell::xCellAnalysis(exp_mtr,rnaseq = TRUE,
                       scale = TRUE,alpha = 0.5, ...)
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


#' @title ConsensusTME deconvolution
#' @description use ConsensusTME to predict TME
#' @param SE an SummarizedExperiment object contains the bulk RNA-seq dataset that you want to use for deconvolution and obtaining its cell fraction.
#' @param ... the arguments
#' @export

epic <- function(SE, ...){
  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)

  EPIC::EPIC(bulk = exp_mtr, ...)$cellFractions
}


#' @title abis deconvolution
#' @description use abis to predict TME
#' @param sig_matrix gene expression matrix from isolated cells.
#' @param SE an SummarizedExperiment object contains the bulk RNA-seq dataset that you want to use for deconvolution and obtaining its cell fraction.
#' @importFrom stats coef
#' @export

ABIS <- function(sig_matrix,SE) {
  isList <- is.list(SE)
  exp_mtr <- bind_mtr(SE, isList)

  genes <- intersect(rownames(exp_mtr), rownames(sig_matrix))
  Dec <- (apply(exp_mtr[genes, , drop = F],
                2,
                function(x) coef(MASS::rlm(as.matrix(signature[genes, ]), x, maxit = 100)))) * 100
  Dec <- signif(Dec, 3)
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

  score.data
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


#' @title MCPCounter deconvolution
#' @description use MCPCounter to predict TME
#' @param SE an SummarizedExperiment object contains the bulk RNA-seq dataset that you want to use for deconvolution and obtaining its cell fraction.
#' @param featuresType type of identifiers for expression features. Defaults to "affy133P2_probesets" for Affymetrix Human Genome 133 Plus 2.0 probesets. Other options are "HUGO_symbols" (Official gene symbols), "ENTREZ_ID" (Entrez Gene ID) or "ENSEMBL_ID" (ENSEMBL Gene ID)
#' @param ... other parameter
#' @export

MCPCounter <- function(SE, featuresType = "HUGO_symbols", ...) {
  feature_types <- match.arg(feature_types,c("affy133P2_probesets","HUGO_symbols","ENTREZ_ID","ENSEMBL_ID"))
  isList <- is.list(SE)
  gene_expression_matrix <- bind_mtr(SE, isList)

  arguments <- rlang::dots_list(gene_expression_matrix, featuresType = feature_types, ..., .homonyms = "last")
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
                            apply(expression[intersect(row.names(expression),x),,drop=F],
                                  2,mean,na.rm=T)
                            }))))
}


#' @title TIMER deconvolution
#' @description use TIMER to predict TME
#' @param exp_mtr matrix or data.frame with features in rows and samples in columns
#' @param type type of cancer
#' @export

TIMER <- function(exp_mtr,type) {
  TIMER.Immune <- readRDS(system.file("extdata", "TIMER.Immune.rds", package = "tigeR", mustWork = TRUE))

  co_genes <- intersect(rownames(exp_mtr),rownames(TIMER.Immune[[1]]))
  pre_rmBE <- cbind(exp_mtr[co_genes,],TIMER.Immune[[1]][co_genes,])
  batch <- as.factor(c(rep("Tumor",ncol(exp_mtr)),rep("Immune",ncol(TIMER.Immune[[1]]))))
  post_rmBE <- sva::ComBat(pre_rmBE, batch)

  tumor_exp <- post_rmBE[,1:ncol(exp_mtr)]
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
