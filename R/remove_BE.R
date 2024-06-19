#' @title Remove batch effect
#' @description remove batch effect.
#' @param SE a SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param method batch effect correction methods.
#' @export

remove_BE <- function(SE,method){
  batch <- unique(colData(SE)$dataset_id)
  if(length(batch) < 2)
    stop("Only one batch is detected in the data. At least two batches are required in batch effect removal.")

  switch(method,
         Combat = tigeR_Combat(SE),
         Combat_seq = tigeR_Combat_seq(SE),
         limma = tigeR_limma(SE),
         DWD = tigeR_DWD(SE))
}

#' @title Remove batch effect.
#' @description Generate a naive bayes model.
#' @param SE a SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @importFrom stats model.matrix

tigeR_Combat <- function(SE){
  model <- stats::model.matrix(~as.factor(colData(SE)$response))
  mtr <- dataPreprocess(assay(SE), rownames(assay(SE)), turn2HL = FALSE)
  final <- sva::ComBat(dat = mtr,batch = as.factor(colData(SE)$dataset_id),mod = model)
  SummarizedExperiment::SummarizedExperiment(assays = list(counts = final),
                                             colData = colData(SE))
}

#' @title Remove batch effect.
#' @description Generate a naive bayes model.
#' @param SE a SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @importFrom stats model.matrix

tigeR_Combat_seq <- function(SE){
  mtr <- dataPreprocess(assay(SE), rownames(assay(SE)), turn2HL = FALSE)
  final <- sva::ComBat_seq(count = mtr,batch = as.factor(colData(SE)$dataset_id),group = colData(SE)$response)
  SummarizedExperiment::SummarizedExperiment(assays = list(counts = final),
                                             colData = colData(SE))
}

#' @title Remove batch effect.
#' @description Generate a naive bayes model.
#' @param SE a SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @importFrom stats model.matrix

tigeR_DWD <- function(SE){
  train_data <- scale(t(na.omit(SummarizedExperiment::assay(SE))))
  train_labels <- SummarizedExperiment::colData(SE)$dataset_id
  train_labels <- ifelse(train_labels==unique(train_labels)[1],1,-1)
  kern=kernlab::vanilladot()
  dwd_model <- kerndwd::kerndwd(train_data, train_labels,kern=kern)
  nfit = kernlab::kernelMult(kern, train_data, t(na.omit(SummarizedExperiment::assay(SE))),
                             dwd_model$alpha[-1, , drop=FALSE])
  nfit = sweep(nfit, MARGIN=2, dwd_model$alpha[1, , drop=FALSE], "+")

  return(nfit)
}

tigeR_limma <- function(SE){
  meta <- SummarizedExperiment::colData(SE)
  design <- model.matrix(~ meta$response) # 创建设计矩阵

  final <-
    limma::removeBatchEffect(SummarizedExperiment::assay(SE),
                             batch = meta$dataset_id, design = design)
  SummarizedExperiment::SummarizedExperiment(assays = list(counts = final),
                                             colData = colData(SE))
}
