#' @title Refine a custom reference matrix
#' @description refine a Custom reference matrix for TME deconvolution from single-cell sequencing data.
#' @param Seurat_obj des
#' @param n_ref_genes description
#' @param logfc.threshold description
#' @param min.pct des
#' @param only.pos description
#' @export

refine_Reference <- function(Seurat_obj,n_ref_genes=50,
                            logfc.threshold = 0.15,
                            min.pct=0.1,only.pos=TRUE){
  Markers <- Seurat::FindAllMarkers(object = Seurat_obj,
                                    assay = "RNA",
                                    logfc.threshold = 0.15,
                                    min.pct = min.pct,
                                    only.pos = only.pos)
  Markers <- Markers %>%
    dplyr::filter(p_val_adj < 0.01)
  Markers$d <- Markers$pct.1 - Markers$pct.2

  Markers <- Markers %>%
    dplyr::filter(d > 0.1)
  Markers <- Markers %>%
    dplyr::arrange(cluster,dplyr::desc(avg_log2FC))

  used_gene <- sort(unique(Markers$gene))
  tpm <- Seurat::GetAssayData(Seurat_obj, layer = "count") %>%
    as.matrix() %>%
    count2tpm()

  sig.list <-
    lapply(levels(Seurat_obj$celltype), function(x){
      apply(tpm[,Seurat_obj$celltype==x],1,mean)
    })
  names(sig.list) <- levels(Seurat_obj$celltype)
  sig_matrix <- do.call("cbind",sig.list)
  selected_gene <- intersect(rownames(sig_matrix),used_gene)
  return(sig_matrix[selected_gene,])
}
