#' @title title
#' @description de
#' @param counts d
#' @export

count2tpm <- function(counts){
  symbols <- rownames(counts) %>%
    sub("\\.$","",.)
  ensembl = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  gene_info = biomaRt::getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id'),
                             filters = 'hgnc_symbol',
                             values = symbols,
                             mart = ensembl)
  Homo_sapiens <- AnnotationDbi::loadDb(system.file("extdata",
                                                    "Homo_sapiens.sqlite",
                                                    package = "tigeR",
                                                    mustWork = TRUE))
  exons_gene <- GenomicFeatures::exonsBy(Homo_sapiens, by = "gene")
  ensembl_id <- intersect(unique(gene_info$ensembl_gene_id),names(exons_gene))
  exons_gene <- exons_gene[ensembl_id]
  exons_gene_lens <- lapply(exons_gene,function(x){sum(GenomicRanges::width(GenomicRanges::reduce(x)))})
  length <- unlist(exons_gene_lens)
  gene_info$length <- length[gene_info$ensembl_gene_id]
  final_info <- na.omit(gene_info)

  genes_to_keep <- intersect(final_info$hgnc_symbol, rownames(counts))
  new_counts <- counts[genes_to_keep,]
  gene_length <- final_info$length
  names(gene_length) <- final_info$hgnc_symbol
  gene_length <- gene_length[!duplicated(names(gene_length))]

  length <- gene_length[rownames(new_counts)]
  rpk <- new_counts / length
  tpm <- t(t(rpk) * 1e6 / colSums(rpk))
  return(tpm)
}
