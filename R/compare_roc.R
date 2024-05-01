#' @title dd
#' @param ROC dd
#' @export

compare_roc <- function(ROC){
  d <-
    lapply(ROC, function(x){
      dplyr::arrange(
        data.frame(
          TPR = x$sensitivities,
          FPR = x$specificities
        ),.data$TPR)
    })
  ggplot() +
    geom_line(data = d[[1]], aes(x = .data$FPR, y = .data$TPR, color = "Curve 1"), linewidth = 1) +
    geom_line(data = d[[2]], aes(x = .data$FPR, y = .data$TPR, color = "Curve 2"), linewidth = 1) +
    geom_line(data = d[[3]], aes(x = .data$FPR, y = .data$TPR, color = "Curve 3"), linewidth = 1) +
    geom_line(data = d[[4]], aes(x = .data$FPR, y = .data$TPR, color = "Curve 4"), linewidth = 1) +
    geom_line(data = d[[5]], aes(x = .data$FPR, y = .data$TPR, color = "Curve 5"), linewidth = 1) +
    geom_line(data = d[[6]], aes(x = .data$FPR, y = .data$TPR, color = "Curve 6"), linewidth = 1) +
    xlim(1, 0) +
    labs(x = "specificity", y = "sensitivity") +
    scale_color_manual(values = c("black", "#F89B9B", "#4DB867", "#C64D6A", "#FDCEBC", "green"),
                       labels = c(paste0("Naive Bayes ",sprintf("%.3f", ROC[[1]]$auc)),
                                  paste0("SVM ",sprintf("%.3f", ROC[[2]]$auc)),
                                  paste0("Random Forest ",sprintf("%.3f", ROC[[3]]$auc)),
                                  paste0("Cancerclass ",sprintf("%.3f", ROC[[4]]$auc)),
                                  paste0("Adaboost ",sprintf("%.3f", ROC[[5]]$auc)),
                                  paste0("Logitboost ",sprintf("%.3f", ROC[[6]]$auc)))) +
    theme(plot.background = element_rect(fill = "transparent", color = "transparent"),
          panel.background = element_rect(fill = "transparent", color = "black", linewidth = 1.5),
          panel.grid = element_blank(),
          legend.position = c(0.7,0.18),
          legend.background = element_rect(fill = "transparent",color = "transparent"),
          legend.box.background = element_rect(fill = "transparent",color = "transparent"),
          legend.key = element_blank(),
          legend.key.height = unit(0,"mm"),
          legend.key.spacing.y = unit(0,"mm"),
          legend.title = element_blank(),
          legend.text = element_text(size=7),
          legend.byrow = TRUE,
          axis.title = element_text(face = "bold", size = 12, color = "black"),
          axis.text = element_text(face = "bold", size = 9, color = "black"),
          aspect.ratio = 1) +
    guides(color=guide_legend(ncol = 1))
}
