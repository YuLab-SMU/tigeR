#' @title Process data before running machine learning algorithm
#' @description Process data before running machine learning algorithm
#' @param pick the number of dataset you want to download. If NULL, return the dataset list.
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom AnnotationHub query
#' @export

Dataloader <- function(pick=NULL){
  list <- c("GBM_PRJNA482620",
            "HNSC_GSE93157",
            "LGG_E_MTAB_6270",
            "MEL_GSE78220",
            "MEL_GSE91061",
            "MEL_GSE93157",
            "MEL_GSE96619",
            "MEL_GSE100797",
            "MEL_GSE106128",
            "MEL_GSE115821",
            "MEL_GSE145996",
            "MEL_Nathanson_2017",
            "MEL_phs000452 samples",
            "MEL_PRJEB23709 samples",
            "nonsqNSCLC_GSE93157",
            "NSCLC_GSE126044",
            "NSCLC_GSE135222",
            "RCC_Braun_2020",
            "RCC_GSE67501",
            "STAD_PRJEB25780")
  if(is.null(pick)){
    return(data.frame(Dataset=list))
  }
  else if(all(pick >= 1) && all(pick <= 20) && is.numeric(pick)){
    dat = ExperimentHub()
    hub = query(dat,"tigeR.data")
    for (i in pick) {
      assign(list[i],suppressWarnings(hub[[paste0("EH",8261+i)]]), envir = .GlobalEnv)
    }
    return()
  }
}
