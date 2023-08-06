#' @title Process data before running machine learning algorithm
#' @description Process data before running machine learning algorithm
#' @param pick the number of dataset you want to download. If NULL, return the dataset list.
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom AnnotationHub query
#' @export

Dataloader <- function(pick=NULL){
  if(is.null(pick)){
    Summarized_Statistics <- NULL
    data(Summarized_Statistics, package = 'tigeR', envir = current_env())
    return(Summarized_Statistics)
  }
  else if(all(pick >= 1) && all(pick <= 20) && is.numeric(pick)){
    dat = ExperimentHub()
    hub = query(dat,"tigeR.data")
    for (i in pick) {
      assign(Summarized_Statistics[i,1],suppressWarnings(hub[[paste0("EH",8261+i)]]), envir = .GlobalEnv)
    }
    return()
  }
}
