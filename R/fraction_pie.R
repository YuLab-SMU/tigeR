#' @title Draw pie charts for TME deconvolution result
#' @description Generate a pie charts illustrating the cell fraction or relative cell abundance for each sample.
#' @param mtr cell fraction matrix
#' @param feature the color levels of the cell types.
#' @param rows the number of rows in the plotting matrix.
#' @param pie_scale amount to scale the pie size if there is no radius mapping exists.
#' @param color the color of the border between each fraction of the pie charts.
#' @param label_radius numeric the radius of label position (relative the radius of pie), default is NULL, when it is provided, the ratio or value label will be displayed.
#' @param label_show_ratio logical only work when label_radius is not NULL, default is TRUE, meaning the ratio of label will be displayed.
#' @param label_threshold numeric the threshold is to control display the label, the ratio of slice pie smaller than the threshold will not be displayed.
#' @param fontsize the font size of the labels.
#' @param ... other parameters
#' @export

fraction_pie <- function(mtr, feature,
                         rows=1,pie_scale=0.9,
                         color="white",
                         label_radius=NULL,
                         label_show_ratio=FALSE,
                         label_threshold=0.04,fontsize=1.5,...){
  mtr <- mtr[levels(feature),]
  colors <- c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c",
              "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5",
              "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f",
              "#c7c7c7", "#bcbd22","#dbdb8d", "#17becf", "#9edae5",
              "#393b79", "#5254a3")
  colors <- colors[seq_along(feature)]
  names(colors) <- levels(feature)
  mt <- as.data.frame(t(mtr))
  mt$x <- rep(seq(from=1,by=2,length.out=ceiling(nrow(mt)/rows)),rows)[1:nrow(mt)]
  mt$y <- rep(rows:1,each=ceiling(nrow(mt)/rows))[1:nrow(mt)] * 2
  mt$sample <- rownames(mt)
  mt$radius <- rep(0.9,nrow(mt))

  ggplot() +
    scatterpie::geom_scatterpie(data = mt,
                                aes(x,y,r=radius),
                                pie_scale= pie_scale,
                                cols = colnames(mt)[1:nrow(mtr)],
                                color = color,
                                label_radius = label_radius,
                                label_show_ratio = label_show_ratio,
                                label_threshold = label_threshold,
                                fontsize=fontsize,...) +
    scale_fill_manual(values = colors) +
    theme_void() +
    theme(legend.position = c(0.5,1.2),
          legend.title = element_blank(),
          legend.text = element_text(size = 18),

          plot.background = element_rect(fill = "transparent",color = "transparent"),
          panel.background = element_rect(fill="transparent",color = "transparent")) +
    coord_map() +
    guides(fill=guide_legend(nrow=1))
}
