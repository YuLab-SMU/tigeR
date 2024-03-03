#' @title Process data before running machine learning algorithm
#' @description Process data before running machine learning algorithm
#' @param pick a number(1-20) or a numeric vector specify the corresponding dataset(s) you wish to load. Alternatively, you can use Dataloader() with pick=NULL to get an overview of all available datasets.
#' @param use_source specify the source of the data to download ("Web Server" or "ExperimentHub")
#' @importFrom rlang current_env
#' @export

Dataloader <- function(pick=NULL, use_source="Web Server"){
  use_source <- match.arg(use_source,c("Web Server", "ExperimentHub"))
  if(is.null(pick)){
    Dataset_info <- NULL
    ev <- current_env()
    data(Dataset_info, package = 'tigeR', envir = ev)
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
#' @import utils
#' @importFrom rlang current_env
#' @export

load_from_WebServer <- function(pick){
  Dataset_info <- NULL
  ev <- current_env()
  data(Dataset_info, package = 'tigeR', envir = ev)
  Dataset_ID <- c("GBM-PRJNA482620","HNSC-GSE93157","LGG_E-MTAB-6270",
                  "Melanoma-GSE78220","Melanoma-GSE91061","Melanoma-GSE93157",
                  "Melanoma-GSE96619","Melanoma-GSE100797","Melanoma-GSE106128",
                  "Melanoma-GSE115821","Melanoma_GSE145996","Melanoma-Nathanson_2017",
                  "Melanoma-phs000452","Melanoma-PRJEB23709","nonsqNSCLC-GSE93157",
                  "NSCLC_GSE126044","NSCLC_GSE135222","RCC-Braun_2020",
                  "RCC-GSE67501","STAD-PRJEB25780")
  for (i in pick) {
    data_exp <- RCurl::getURL(paste0("http://tiger.canceromics.org/tiger/Download/immunotherapy/expression/tsv/",
                                     Dataset_ID[i],".Response.tsv"),
                              timeout = 2000)
    exp <- read.delim(textConnection(data_exp), sep="\t")
    expr_mtr <- as.matrix(exp[,-1])
    rownames(expr_mtr) <- exp[,1]

    data_col <- RCurl::getURL(paste0("http://tiger.canceromics.org/tiger/Download/immunotherapy/clinical/tsv/",
                                     Dataset_ID[i],".Response.tsv"),
                              timeout = 2000)
    col_data <- read.delim(textConnection(data_col),
                           sep = "\t",na.strings = c("NA","#N/A"))
    rownames(col_data) <- col_data[,1]

    SE_obj <- SummarizedExperiment::SummarizedExperiment(assays=S4Vectors::SimpleList(expr_mtr),
                                                         colData=S4Vectors::DataFrame(col_data),
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
#' @importFrom rlang current_env
#' @export

load_from_ExperimentHub <- function(pick){
  Dataset_info <- NULL
  ev <- current_env()
  data(Dataset_info, package = 'tigeR', envir = ev)
  dat <- ExperimentHub()
  hub <- query(dat,"tigeR.data")
  for (i in pick) {
    assign(Dataset_info[i,1],
           hub[[paste0("EH",8261+i)]],
           envir = .GlobalEnv)
  }
  return()
}
