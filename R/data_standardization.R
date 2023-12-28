#' @title Build machine learning prediction model for immunotherapy response
#' @description Generate immunotherapy prognosis prediction model.
#' @param SE the dataset you wish to use to build your model. A SummarizedExperiment (SE) object, which can be either a single SE object or a list of SE objects. Note that for each SE object, the colData must contain treatment information under the column name Treatment.
#' @export

data_standardization <- function(SE){
  to1(log2SE(turNA(fpkm2tpm(SE))))
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

to1 <- function(SE){
  exp <- SummarizedExperiment::assay(SE)
  exp <- t(
  apply(exp,1,function(x){
    if(!all(is.na(x))){
      return((x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE)))
    }else{
      x
    }
  }))
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
