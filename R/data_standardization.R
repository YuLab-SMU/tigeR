#' @title Build machine learning prediction model for immunotherapy response
#' @description Generate immunotherapy prognosis prediction model.
#' @param SE the dataset you wish to use to build your model. A SummarizedExperiment (SE) object, which can be either a single SE object or a list of SE objects. Note that for each SE object, the colData must contain treatment information under the column name Treatment.
#' @param type can be c(1,2,3,4)
#' @export

data_standardization <- function(SE,type){
  batch <- unique(colData(SE)$dataset_id)
  if(length(batch)==1){
    for(i in type){
      SE <-
        switch (i,
                "1" = fpkm2tpm(SE),
                "2" = turNA(SE),
                "3" = log2SE(SE),
                "4" = to1(SE))
    }
  } else {
    SE_list <-
    lapply(seq_along(batch), function(i){
      subset_SE <- SE[, colData(SE)$dataset_id == batch[i]]
      ds_single(subset_SE, type)
    })
    SE <- BiocGenerics::Reduce(function(x, y) SummarizedExperiment::cbind(x, y), SE_list)
  }
  return(SE)
}

ds_single <- function(SE,type){
  for(i in type){
    SE <-
      switch (i,
              "1" = fpkm2tpm(SE),
              "2" = turNA(SE),
              "3" = log2SE(SE),
              "4" = to1(SE))
  }
  return(SE)
}


#' @title count geneset score by different method
#' @description wait to write
#' @param SE an expression matrix.
#' @export

fpkm2tpm <- function(SE)
{
  fpkm <- assay(SE)
  SummarizedExperiment::assay(SE) <- exp(log(fpkm) - log(sum(fpkm,na.rm=TRUE)) + log(1e6))
  SE
}


#' @title Build machine learning prediction model for immunotherapy response
#' @description Generate immunotherapy prognosis prediction model.
#' @param SE the dataset you wish to use to build your model. A SummarizedExperiment (SE) object, which can be either a single SE object or a list of SE objects. Note that for each SE object, the colData must contain treatment information under the column name Treatment.
#' @param by_row if TRUE, scaling by row.
#' @export

to1 <- function(SE,by_row=TRUE){
  exp <- SummarizedExperiment::assay(SE)
  exp <- t(
  apply(exp,ifelse(by_row,1,2),function(x){
    if(!all(is.na(x))){
      return((x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE)))
    }else{
      x
    }
  }))
  if(!by_row)
    exp <- t(exp)
  colnames(exp) <- colnames(assay(SE))
  SummarizedExperiment::assay(SE) <- exp
  SE
}


#' @title Build machine learning prediction model for immunotherapy response
#' @description Generate immunotherapy prognosis prediction model.
#' @param SE the dataset you wish to use to build your model. A SummarizedExperiment (SE) object, which can be either a single SE object or a list of SE objects. Note that for each SE object, the colData must contain treatment information under the column name Treatment.

log2SE <- function(SE){
  exp <- SummarizedExperiment::assay(SE)
  exp <- t(
    apply(exp,1,function(x){
      tmp <- log2(x)
      tmp[tmp==-Inf] <- NA
      tmp
    }))
  colnames(exp) <- colnames(assay(SE))
  SummarizedExperiment::assay(SE) <- exp
  SE
}


#' @title Build machine learning prediction model for immunotherapy response
#' @description Generate immunotherapy prognosis prediction model.
#' @param SE the dataset you wish to use to build your model. A SummarizedExperiment (SE) object, which can be either a single SE object or a list of SE objects. Note that for each SE object, the colData must contain treatment information under the column name Treatment.

turNA <- function(SE){
  exp <- SummarizedExperiment::assay(SE)
  exp <- t(
    apply(exp, 1, function(x){
      if(all(is.na(x))){
        return(x)
      }else if(all(x[!is.na(x)]==0)){
        return(rep(NA,length(x)))
      }else{
        return(x)
      }
    }))
  colnames(exp) <- colnames(assay(SE))
  SummarizedExperiment::assay(SE) <- exp
  SE
}
