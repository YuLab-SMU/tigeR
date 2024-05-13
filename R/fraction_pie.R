#' @title draw integrate fraction pie
#' @description generate a pie plot illustrating the cell fraction or relative cell abundance for each sample.
#' @param mtr cell fraction matrix
#' @param feature the factor of color levels
#' @param rows row numbers of pie plot.
#' @export

fraction_pie <- function(mtr, feature, rows=1){
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

  ggplot() +
    scatterpie::geom_scatterpie(data = mt,
                                aes(x,y,r=0.90),
                                cols = colnames(mt)[1:nrow(mtr)],
                                color = "white") +
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
