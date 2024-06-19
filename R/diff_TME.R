#' @title Visualize differential analysis cell_fraction (Responder vs NonResponder or Pre-Treatment vs Post-Treatment).
#' @description Visualize differential analysis cell_fraction (Responder vs NonResponder or Pre-Treatment vs Post-Treatment).
#' @param SE an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.
#' @param feature the feature or features you are interested in.
#' @param style description
#' @param group_color description
#' @export

diff_TME <- function(SE,feature=NULL,style="elegant",
                     group_color=c("#5f96e8CC", "#ee822fCC")){
  if(is.null(feature)){
    feature <- rownames(SE)
  }

  TME_data <- as.data.frame(t(assay(SE)))

  TME_data$group <- SE$response_NR
  TME_data$sample <- rownames(TME_data)

  TME_New <- reshape2::melt(TME_data)

  colnames(TME_New) <- c("Group","Sample","Celltype","Composition")

  plot_order <- TME_New[TME_New$Group=="R",] %>%
    dplyr::group_by(.data$Celltype) %>%
    dplyr::summarise(m = stats::median(.data$Composition)) %>%
    dplyr::arrange(dplyr::desc(.data$m)) %>%
    dplyr::pull(.data$Celltype)

  TME_New$Celltype <- factor(TME_New$Celltype,levels=plot_order)
  TME_New<-TME_New[!TME_New$Group=='UNK',]
  rs <- c()

  warning_status <- 0
  for (i in levels(TME_New$Celltype)) {
    m <- TME_New[TME_New$Celltype==i,]
    R_S <- m[m$Group == 'R', 4]
    N_S <- m[m$Group == "N", 4]
    if(any(R_S%in%N_S))
      warning_status <- 1
    rs <- c(rs, suppressWarnings(stats::wilcox.test(R_S, N_S)$p.value))
  }
  if(warning_status)
    message("There are identical relative abundance values in groups R and N for the '",i,"'. The p value of the Wilcoxon signed-rank test may not be precise due to ties in the data.")

  #selected_cells <- levels(TME_New$Celltype)[which(rs < 0.05)]
  selected_cells <- levels(TME_New$Celltype)

  if(length(selected_cells) < 5){
    names(rs) <- seq_along(rs)
    selected_cells <- levels(TME_New$Celltype)[as.numeric(names(sort(rs)[1:5]))]
  }

  ciber_theme <- ggplot2::theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                                axis.title = element_text(size = 10,color ="black"),
                                axis.text = element_text(size= 10,color = "black"),
                                axis.text.x = element_text(angle = 45, hjust = 1 ),
                                legend.position = "top",
                                legend.text = element_text(size= 12),
                                legend.title= element_text(size= 12))
  if(style == 'raw'){
    box_TME <-
      ggplot2::ggplot(TME_New[TME_New$Celltype%in%selected_cells,], aes(x = .data$Celltype, y = .data$Composition)) +
      ggplot2::labs(y="Cell composition",x= NULL,title = "TME Cell composition") +
      ggplot2::geom_boxplot(aes(fill = .data$Group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0) +
      ggplot2::scale_fill_manual(values = group_color) +
      ggplot2::theme_classic() + ciber_theme
    y_max <- max(ggplot_build(box_TME)$data[[1]]$ymax)
    box_TME <-
      box_TME +
      ggpubr::stat_compare_means(ggplot2::aes(group = .data$Group),
                                 label = "p.signif",
                                 method = "wilcox.test",
                                 hide.ns = TRUE,
                                 label.y = y_max*1.1) +
      coord_cartesian(ylim = c(0, y_max*1.1))
  }
  if(style == 'elegant'){
    box_TME <- ggplot(TME_New[TME_New$Celltype %in% feature,],
                      aes(x=.data$Celltype,y=.data$Composition,fill=.data$Group)) +
      stat_boxplot(geom = "errorbar",width = 1, color = "black",linetype = "solid",
                   position = position_dodge(0.8),linewidth = 0.7) +
      stat_boxplot(geom = "boxplot", color = "black",linetype = "solid",
                   position = position_dodge(0.8),linewidth = 0.7,
                   width = 0.8, outlier.shape= 19) +
      ggpubr::theme_classic2() + ciber_theme +
      scale_fill_manual(values = group_color)
    y_max <- max(ggplot_build(box_TME)$data[[1]]$ymax)
    box_TME <- box_TME +
      ggpubr::stat_compare_means(ggplot2::aes(group = .data$Group),
                                 label = "p.signif", method = "wilcox.test", hide.ns = FALSE,
                                 label.y = y_max * 1.1,size = 3) + coord_cartesian(ylim = c(0, y_max * 1.1))
  }
  box_TME
}
