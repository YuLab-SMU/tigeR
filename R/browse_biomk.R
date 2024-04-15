#' @title dd
#' @param SE dd
#' @param feature dd
#' @param colors description
#' @param PT_drop description
#' @export

browse_biomk <-
  function(SE, feature,colors=NULL,PT_drop=FALSE){
    isList <- is.list(SE)
    exp_mtr <- bind_mtr(SE, isList)
    meta <- bind_meta(SE, isList)

    idx <- rownames(meta[meta$response_NR!="UNK",])
    exp_mtr_auc <- exp_mtr[,idx]
    meta_auc <- meta[idx,]
    auc <-
      apply(exp_mtr_auc,1,function(x){
        ROC <- pROC::roc(meta_auc$response_NR, x)
        ROC$auc
      })

    df <- data.frame(Cell_type=names(auc),AUC=auc)
    df1 <- dplyr::arrange(df,AUC)
    if(nrow(df1)>10)
      df1 <- df1[(nrow(df1)-9):nrow(df1),]
    df1$Cell_type <- factor(df1$Cell_type,levels = df1$Cell_type)

    if(is.null(colors))
      colors <- c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c",
                  "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5",
                  "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f",
                  "#c7c7c7", "#bcbd22","#dbdb8d", "#17becf", "#9edae5",
                  "#393b79", "#5254a3")
    bar <-
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

    if(nrow(final)>10)
      final <- final[final$cell%in%df1$Cell_type,]

    dot <-
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
            legend.background = element_rect(fill = "transparent",color = "transparent"),
            legend.key = element_rect(color = "transparent"),
            legend.key.size = unit(0,"cm"),
            legend.spacing.x = unit(0,"cm"),
            legend.spacing.y = unit(0,"cm"),
            axis.text = element_text(face = "bold", size = "10", color = "black"),
            axis.title = element_text(face = "bold", size = "12", color = "black"),
            aspect.ratio = 2)
    bar + dot
  }




