#' @title dd
#' @param mtr description
#' @export

cell_name_filter <- function(mtr){
  rn <- rownames(mtr)
  rn <- sub("T\\.cell\\.CD8.","T cells CD8",rn) %>%
    sub("T\\.cells\\.CD8","T cells CD8",.) %>%
    sub("T_cells_CD8","T cells CD8",.,fix = TRUE) %>%
    sub("T cells regulatory \\(Tregs\\)","Tregs",.) %>%
    sub("T\\_regulatory\\_cells","Tregs",.) %>%
    sub("T\\.cell\\.CD4.","T cells CD4",.) %>%
    sub("T\\.cells\\.CD4","T cells CD4",.) %>%
    sub("T\\_cells\\_CD4","T cells CD4",.) %>%
    sub("T cells gamma delta","Yd T cells",.) %>%
    sub("T\\_cells\\_gamma\\_delta","Yd T cells",.) %>%
    sub("^B.cell$","B cells",.) %>%
    sub("B\\-cells","B cells",.) %>%
    sub("B.cells","B cells",.,fixed = TRUE) %>%
    sub("B\\_cells","B cells",.) %>%
    sub("Plasma\\_cells","Plasma cells",.) %>%
    sub("Myeloid\\.dentritic\\.cell","mDCs",.) %>%
    sub("Myeloid dendritic cells","mDCs",.) %>%
    sub("Mast\\_cells","Mast cells",.) %>%
    sub("NK\\_cells","NK cells",.) %>%
    sub("Dendritic cells","DCs",.) %>%
    sub("Dendritic\\.cells","DCs",.) %>%
    sub("Cytotoxic lymphocytes","CTLs",.) %>%
    sub("T cells CD4 memory activated","CD4+ T memory activated",.) %>%
    sub("Immune\\_Score","Immune Score",.)
  rownames(mtr) <- rn
  mtr
}

bar_auc <-
  function(SE, feature,colors=NULL){
    isList <- is.list(SE)
    exp_mtr <- bind_mtr(SE, isList)
    meta <- bind_meta(SE, isList)

    idx <- rownames(meta[meta$response_NR!="UNK",])
    exp_mtr <- exp_mtr[,idx]
    meta <- meta[idx,]
    auc <-
      apply(exp_mtr,1,function(x){
        ROC <- pROC::roc(meta$response_NR, x)
        ROC$auc
      })

    df <- data.frame(Cell_type=names(auc),AUC=auc)
    df1 <- dplyr::arrange(df,AUC)
    if(nrow(df1)>22)
      df1 <- df1[(nrow(df1)-21):nrow(df1),]
    df1$Cell_type <- factor(df1$Cell_type,levels = df1$Cell_type)

    if(is.null(colors))
      colors <- c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c",
                  "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5",
                  "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f",
                  "#c7c7c7", "#bcbd22","#dbdb8d", "#17becf", "#9edae5",
                  "#393b79", "#5254a3")
    ggplot(df1, aes(x = Cell_type, y = AUC,fill=Cell_type)) +
      geom_col(fill = colors[seq_along(df1$Cell_type)]) +
      geom_text(aes(label = round(AUC, 2)), hjust=1.1,vjust = 0.4,color="white") +
      coord_flip() +
      theme(plot.background = element_rect(fill = "transparent",color = "transparent"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "transparent"),
            panel.border = element_rect(fill = "transparent",linewidth = 1.5),
            axis.title.x = element_text(face = "bold", size = "12", color = "black"),
            axis.title.y = element_blank(),
            axis.text = element_text(face = "bold", size = "10", color = "black"),
            aspect.ratio = 2)
  }
dot_Surv <-
  function(SE,PT_drop=FALSE){
    isList <- is.list(SE)
    exp_mtr <- bind_mtr(SE, isList)
    meta <- bind_meta(SE, isList)

    idx_UT <- seq_along(meta[,1])
    if(PT_drop)
      idx_UT <- which(meta$Treatment == 'PRE')

    if(length(idx_UT) == 0)
      stop("Only Untreated patients can be use to perform survival analysis!")

    exp_mtr <- exp_mtr[,idx_UT,drop=FALSE]
    meta <- meta[idx_UT,,drop=FALSE]

    Score <- t(
      apply(exp_mtr, 1, function(x){
        thres <- ifelse(median(x)==0,mean(x),median(x))
        ifelse(x>=thres,1,0)
      }))

    time <- as.numeric(meta$overall.survival..days.)
    status <- sub('Dead','1', meta$vital.status) %>% sub('Alive','0',.)

    dt <-
      apply(Score, 1, function(x){
        df <- data.frame(time,status,Score=x) %>%
          stats::na.omit() %>%
          lapply(as.numeric) %>%
          as.data.frame()
        fit <- survfit(Surv(time, status) ~ Score, data = df)
        cox_md <- coxph(Surv(time, status) ~ Score, data = df)
        summary_cox <- summary(cox_md)
        c(summary_cox$conf.int[,1],summary_cox$coefficients[,5])
      })
    rownames(dt) <- c("HR","P")

    final <- data.frame(
      x = unlist(dt[1,]),
      y = -log10(unlist(dt[2,])),
      cell = colnames(dt)
    )
    final1 <- dplyr::arrange(final,)
    if(nrow(final)>22)
      df1 <- df1[(nrow(df1)-21):nrow(df1),]
    ggplot(final, aes(x = x, y = y, color=cell)) +
      geom_point(size=1.5) +
      scale_y_continuous(position = "right") +
      xlab("Hazard ratio") +
      ylab("-Log10 P value") +
      theme(plot.background = element_rect(fill = "transparent",color = "transparent"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill="transparent"),
            panel.border = element_rect(fill = "transparent",linewidth = 1.5),
            legend.title = element_blank(),
            legend.text = element_text(size = 4),
            legend.position = c(0.79,0.12),
            legend.background = element_rect(fill = "transparent"),
            legend.key.size = unit(0,"cm"),
            legend.spacing.x = unit(0,"cm"),
            legend.spacing.y = unit(0,"cm"),
            axis.text = element_text(face = "bold", size = "10", color = "black"),
            axis.title = element_text(face = "bold", size = "12", color = "black"),
            aspect.ratio = 2)
  }
