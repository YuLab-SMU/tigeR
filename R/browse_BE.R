#' @title Perform PCA to visualize batch effect
#' @description perform PCA to visualize batch effect.
#' @param SE a SummarizedExperiment(SE) object or a list consists of multiple SE objects. The colData of the SE object(s) must contain treatment information named Treatment.
#' @export

browse_BE <- function(SE){
  scale_mtr <- t(apply(na.omit(SummarizedExperiment::assay(SE)),1,scale))
  pca_result <- stats::prcomp(t(na.omit(scale_mtr)), scale = TRUE)
  summary(pca_result)
  meta <- SummarizedExperiment::colData(SE)
  batch <- as.factor(meta$dataset_id)

  df <- pca_result$x
  ggplot(df, aes(.data$PC1,.data$PC2, col = batch, shape = batch)) +
    stat_ellipse(aes(fill = batch), type = "norm", geom ="polygon", alpha = 0.2) +
    scale_shape_manual(values = c(19, 17)) +
    scale_color_manual(values = c('#E69F00', '#56B4E9')) +
    geom_point(size = 1) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          panel.border = element_blank()) +
    labs(x = "PCA_1", y  ="PCA_2", title = "PCA")
}
