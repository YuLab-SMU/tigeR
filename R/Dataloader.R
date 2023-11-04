#' @title Process data before running machine learning algorithm
#' @description Process data before running machine learning algorithm
#' @param pick a number(1-20) or a numeric vector specify the corresponding dataset(s) you wish to load. Alternatively, you can use Dataloader() with pick=NULL to get an overview of all available datasets.
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom AnnotationHub query
#' @export

Dataloader <- function(pick=NULL){
  if(is.null(pick)){
    Dataset_info <- NULL
    data(Dataset_info, package = 'tigeR', envir = current_env())
    return(Dataset_info)
  }
  else if(all(pick >= 1) && all(pick <= 20) && is.numeric(pick)){
    dat = ExperimentHub()
    hub = query(dat,"tigeR.data")
    for (i in pick) {
      assign(Dataset_info[i,1],suppressWarnings(hub[[paste0("EH",8261+i)]]), envir = .GlobalEnv)
    }
    return()
  }
}
