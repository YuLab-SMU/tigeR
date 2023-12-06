#' @title Process data before running machine learning algorithm
#' @description Process data before running machine learning algorithm
#' @param pick a number(1-20) or a numeric vector specify the corresponding dataset(s) you wish to load. Alternatively, you can use Dataloader() with pick=NULL to get an overview of all available datasets.
#' @param use_source specify the source of the data to download ("Web Server" or "ExperimentHub")
#' @export

Dataloader <- function(pick=NULL, use_source="Web Server"){
  if(is.null(pick)){
    Dataset_info <- NULL
    data(Dataset_info, package = 'tigeR', envir = current_env())
    return(Dataset_info)
  }

  if(any(!pick %in% 1:20))
    stop("The parameter pick must be an integer between 1 and 20 or a vector composed of integers within this range.")

  if(use_source == "Web Server"){
    load_from_WebServer(pick)
  }else if(use_source == "ExperimentHub"){
    load_from_ExperimentHub(pick)
  }else{
    stop("The parameter use_source must be 'Web Server' or 'ExperimentHub'.")
  }
}


#' @title Process data before running machine learning algorithm
#' @description Process data before running machine learning algorithm
#' @param pick a number(1-20) or a numeric vector specify the corresponding dataset(s) you wish to load. Alternatively, you can use Dataloader() with pick=NULL to get an overview of all available datasets.
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors SimpleList
#' @importFrom S4Vectors DataFrame
#' @export

load_from_WebServer <- function(pick){
  Dataset_info <- NULL
  data(Dataset_info, package = 'tigeR', envir = current_env())
  Dataset_ID <- c("GBM-PRJNA482620","HNSC-GSE93157","LGG_E-MTAB-6270",
                  "Melanoma-GSE78220","Melanoma-GSE91061","Melanoma-GSE93157",
                  "Melanoma-GSE96619","Melanoma-GSE100797","Melanoma-GSE106128",
                  "Melanoma-GSE115821","Melanoma_GSE145996","Melanoma-Nathanson_2017",
                  "Melanoma-phs000452","Melanoma-PRJEB23709","nonsqNSCLC-GSE93157",
                  "NSCLC_GSE126044","NSCLC_GSE135222","RCC-Braun_2020",
                  "RCC-GSE67501","STAD-PRJEB25780")
  for (i in pick) {
    exp <- read.table(paste0("http://tiger.canceromics.org/tiger/Download/immunotherapy/expression/tsv/",
                                  Dataset_ID[i],".Response.tsv"), sep="\t")
    expr_mtr <- as.matrix(exp[,-1])
    rownames(expr_mtr) <- exp[,1]
    col_data <- read.table(paste0("http://tiger.canceromics.org/tiger/Download/immunotherapy/clinical/tsv/",
                                  Dataset_ID[i],".Response.tsv"), sep="\t")
    rownames(col_data) <- col_data[,1]

    SE_obj <- SummarizedExperiment(assays=SimpleList(expr_mtr),
                                   colData=DataFrame(col_data),
                                   checkDimnames=TRUE)
    assign(Dataset_info[i,1],
           SE_obj,
           envir = .GlobalEnv)
  }
  return()
}


#' @title Process data before running machine learning algorithm
#' @description Process data before running machine learning algorithm
#' @param pick a number(1-20) or a numeric vector specify the corresponding dataset(s) you wish to load. Alternatively, you can use Dataloader() with pick=NULL to get an overview of all available datasets.
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom AnnotationHub query
#' @export

load_from_ExperimentHub <- function(pick){
  Dataset_info <- NULL
  data(Dataset_info, package = 'tigeR', envir = current_env())
  dat = ExperimentHub()
  hub = query(dat,"tigeR.data")
  for (i in pick) {
    assign(Dataset_info[i,1],
           hub[[paste0("EH",8261+i)]],
           envir = .GlobalEnv)
  }
  return()
}
