#' @title Standardize the cell names in TME deconvolution
#' @description standarize the cell names in TME deconvolution.
#' @param SE a SummarizedExperiment object contains the cell fraction of deconvolution.
#' @export

cell_name_standardization <- function(SE){
  mtr <- assay(SE)
  rnames <- rownames(mtr)
  ## TIMER
  rnames[rnames == "T.cell.CD8."] <- "CD8T"
  rnames[rnames == "B.cell"] <- "B cell"
  rnames[rnames == "Myeloid.dentritic.cell"] <- "DC"
  rnames[rnames == "Macrophage"] <- "M\u03a6"
  rnames[rnames == "Neutrophil"] <- "NEU"
  rnames[rnames == "T.cell.CD4."] <- "CD4T"

  ## CIBERSORT Neutrophils
  rnames[rnames == "B cells naive"] <- "Naive B"
  rnames[rnames == "T cells regulatory (Tregs)"] <- "Treg"
  rnames[rnames == "Neutrophils"] <- "NEU"
  rnames[rnames == "Mast cells resting"] <- "RMC"
  rnames[rnames == "T cells gamma delta"] <- "\u03b3\u03b4 T"
  rnames[rnames == "Macrophages M2"] <- "M2 M\u03a6"
  rnames[rnames == "T cells CD4 naive"] <- "naive CD4T"
  rnames[rnames == "T cells CD8"] <- "CD8T"
  rnames[rnames == "Dendritic cells resting"] <- "RDC"

  ## MCPCounter
  rnames[rnames == "CD8 T cells"] <- "CD8T"
  rnames[rnames == "Endothelial cells"] <- "Endo"
  rnames[rnames == "NK cells"] <- "NK"
  rnames[rnames == "Monocytic lineage"] <- "Mo lineage"
  rnames[rnames == "Myeloid dendritic cells"] <- "DC"
  rnames[rnames == "T cells"] <- "T cell"
  rnames[rnames == "Fibroblasts"] <- "Fibro"
  rnames[rnames == "Cytotoxic lymphocytes"] <- "CTL"

  ## xCell
  rnames[rnames == "Smooth muscle"] <- "SM"
  rnames[rnames == "mv Endothelial cells"] <- "mvEndo"
  rnames[rnames == "Eosinophils"] <- "EOS"
  rnames[rnames == "ly Endothelial cells"] <- "lyEndo"
  rnames[rnames == "CD8+ Tcm"] <- "CD8Tcm"
  rnames[rnames == "Megakaryocytes"] <- "MK"
  rnames[rnames == "Osteoblast"] <- "OB"
  rnames[rnames == "CD8+ Tem"] <- "CD8Tem"
  rnames[rnames == "Tgd cells"] <- "\u03b3\u03b4 T"

  # EPIC
  rnames[rnames == "otherCells"] <- "other"
  rnames[rnames == "Endothelial"] <- "Endo"
  rnames[rnames == "Macrophages"] <- "M\u03a6"
  rnames[rnames == "Bcells"] <- "B cell"
  rnames[rnames == "CD8_Tcells"] <- "CD8T"
  rnames[rnames == "CD4_Tcells"] <- "CD4T"
  rnames[rnames == "NKcells"] <- "NK"

  # ESTIMATE
  rnames[rnames == "ImmuneScore"] <- "Immu-Sc"
  rnames[rnames == "TumorPurity"] <- "T-Purity"
  rnames[rnames == "ESTIMATEScore"] <- "ES-Sc"
  rnames[rnames == "StromalScore"] <- "Str-Sc"

  # ABIS替换
  rnames[rnames == "T cell CD8"] <- "CD8T"
  rnames[rnames == "NK cells activated"] <- "aNK"
  rnames[rnames == "Dendritic cells activated"] <- "aDC"
  rnames[rnames == "T cells follicular helper"] <- "Tfh"
  rnames[rnames == "T cells CD4 memory activated"] <- "aCD4Tm"
  rnames[rnames == "Plasma cells"] <- "PC"
  rnames[rnames == "Macrophages M1"] <- "M1 M\u03a6"

  # ConsensusTME
  rnames[rnames == "B_cells"] <- "B cell"
  rnames[rnames == "Plasma_cells"] <- "PC"
  rnames[rnames == "T_regulatory_cells"] <- "Treg"
  rnames[rnames == "Immune_Score"] <- "Immu-Sc"
  rnames[rnames == "T_cells_CD8"] <- "CD8T"
  rnames[rnames == "Mast_cells"] <- "MC"
  rnames[rnames == "T_cells_gamma_delta"] <- "\u03b3\u03b4 T"
  rnames[rnames == "NK_cells"] <- "NK"
  rnames[rnames == "T_cells_CD4"] <- "CD4T"

  # quanTIseq
  rnames[rnames == "T.cells.CD8"] <- "CD8T"
  rnames[rnames == "Macrophages.M1"] <- "M1 M\u03a6"
  rnames[rnames == "Macrophages.M2"] <- "M2 M\u03a6"
  rnames[rnames == "B.cells"] <- "B cell"
  rnames[rnames == "Tregs"] <- "Treg"
  rnames[rnames == "Dendritic.cells"] <- "DC"
  rnames[rnames == "T.cells.CD4"] <- "CD4T"

  rownames(mtr) <- rnames
  SE_new <- SummarizedExperiment::SummarizedExperiment(
    assays = S4Vectors::SimpleList(mtr = mtr),
    colData = SummarizedExperiment::colData(SE))
  SE_new
}
