#' @title title
#' @param mtr description
#' @param feature description
#' @export

draw_bar <- function(mtr,feature){
  mtr <- mtr[levels(feature),]
  colors <- c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c",
              "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5",
              "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f",
              "#c7c7c7", "#bcbd22","#dbdb8d", "#17becf", "#9edae5",
              "#393b79", "#5254a3")
  colors <- colors[seq_along(feature)]
  names(colors) <- levels(feature)
  result <-
  lapply(seq_along(mtr[,1]), function(i){
    data <- data.frame(
      sample=colnames(mtr),
      value=unlist(mtr[i,])
    )
    ggplot(data,aes(x=.data$sample, y=.data$value)) +
      geom_bar(stat = "identity",fill=colors[i]) +
      xlab("Samples") +
      ylab(names(colors)[i]) +
      theme(panel.grid = element_blank(),
            panel.background = element_rect(fill = "transparent",
                                            color="transparent"),
            legend.position = "none",
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(color="black",face="bold"),
            axis.title.x = element_text(face="bold",size=16),
            axis.title.y = element_text(face="bold",size=16),
            axis.line = element_line(linewidth = 1))
  })
}
