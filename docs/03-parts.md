# 🌱 Tumor Microenviroment Analysis

## Availiable TME Analysis Method in tigeR
<div style="width:780px; height:200px; overflow-y: scroll; overflow-x: hidden;">
|Algorithm |license |PMID|
|:----------:|:--------------------:|:--------------:|
|[TIMER](http://cistrome.org/TIMER/)|free (GPL 2.0)|[27549193](https://pubmed.ncbi.nlm.nih.gov/27549193/)|
|[CIBERSORT](https://cibersort.stanford.edu/)|free for non-commercial use only|[25822800](https://pubmed.ncbi.nlm.nih.gov/25822800/)|
|[MCPCounter](https://github.com/ebecht/MCPcounter) |free ([GPL 3.0](https://github.com/ebecht/MCPcounter/blob/master/Source/License))|[27765066](https://pubmed.ncbi.nlm.nih.gov/27765066/)|
|[xCell](http://xcell.ucsf.edu/)|free ([GPL 3.0](https://github.com/dviraran/xCell/blob/master/DESCRIPTION))|[29141660](https://pubmed.ncbi.nlm.nih.gov/29141660/)|
|[IPS](https://github.com/icbi-lab/Immunophenogram)| free (BSD)|[28052254](https://pubmed.ncbi.nlm.nih.gov/28052254/)|
|[EPIC](https://gfellerlab.shinyapps.io/EPIC_1-1/)|free for non-commercial use only ([Academic License](https://github.com/GfellerLab/EPIC/blob/master/LICENSE))|[29130882](https://pubmed.ncbi.nlm.nih.gov/29130882/)|
|[ESTIMATE](https://gfellerlab.shinyapps.io/EPIC_1-1/)| free ([GPL 2.0](https://bioinformatics.mdanderson.org/public-software/estimate/))|[24113773](https://pubmed.ncbi.nlm.nih.gov/24113773/)|
|[ABIS](https://giannimonaco.shinyapps.io/ABIS/)|free ([GPL 2.0](https://github.com/giannimonaco/ABIS))|[30726743](https://pubmed.ncbi.nlm.nih.gov/30726743/)|
|[ConsensusTME](https://olliecast.shinyapps.io/Deconvolution_Benchmarking/)|free ([GPL 3.0](https://github.com/cansysbio/ConsensusTME/blob/master/LICENSE.md))|[31641033](https://pubmed.ncbi.nlm.nih.gov/31641033/)|
|[quanTIseq](http://icbi.at/software/quantiseq/doc/index.html)|free (BSD)|[31126321](https://pubmed.ncbi.nlm.nih.gov/31126321/)|
</div>

## Derive the Proportions of Different TME Cell Types
 You can use the function **deconv_TME()** to derive the proportions of different tumor microenvironment cell types from gene expression data with these tools.

<div style="width:780px; height:200px; overflow-y: scroll; overflow-x: hidden;">
```
packages <- c("xCell", "EPIC", "ConsensusTME", "quantiseqr")
for (package in packages) {
  if (!require(package, character.only = TRUE)) {
    BiocManager::install(package)
  }
}

## TIMER
frac1 <- deconv_TME(MEL_GSE78220,method="TIMER")

## CIBERSORT
frac2 <- deconv_TME(MEL_GSE78220,method="CIBERSORT")

## MCPCounter
frac3 <- deconv_TME(MEL_GSE78220,method="MCPCounter")

## xCell
frac4 <- deconv_TME(MEL_GSE78220,method="xCell")

## IPS
frac5 <- deconv_TME(MEL_GSE78220,method="IPS")

## EPIC
frac6 <- deconv_TME(MEL_GSE78220,method="epic")

## ESTIMATE
frac7 <- deconv_TME(MEL_GSE78220,method="ESTIMATE")

## ABIS
frac8 <- deconv_TME(MEL_GSE78220,method="ABIS")

## ConsensusTME
frac9 <- deconv_TME(MEL_GSE78220,method="ConsensusTME")

## quanTIseq
frac10 <- deconv_TME(MEL_GSE78220,method="quanTIseq")
```
</div>

## Visualization and Comparing the Cell Proportions
```
cell1 <- c("T cells CD4","Neutrophil", "Macrophage","mDCs","B cells", "T cells CD8")
pie1 <- fraction_pie(cell_name_filter(frac1),feature=factor(cell1, levels = cell1))

cell2 <- c("DCs resting", "T cells CD8", "T cells CD4 naive", "Macrophages M2", "Yd T cells", "Monocytes","Mast cells resting", "Neutrophils", "Tregs","B cells naive")
pie2 <- fraction_pie(cell_name_filter(frac2[[1]][1:22,]),feature=factor(cell2, levels = cell2))
```

## Searching for Key Cell Types Influencing Immune Therapy Response
```
## TIMER
TM <- deconv_TME(MEL_GSE91061,method = "TIMER")
TM_SE <- SummarizedExperiment(assays=SimpleList(TM),
                               colData=colData(MEL_GSE91061))
browse_biomk(SE=TM_SE)
```

## 📝 More Details about TME Analysis{-}
 **TIMER** is a comprehensive tool for systematical analysis of immune infiltrates across diverse cancer types.

 **CIBERSORT** is an analytical tool from the Alizadeh Lab and Newman Lab to impute gene expression profiles and provide an estimation of the abundances of member cell types in a mixed cell population, using gene expression data.

 **xCell** is a gene signatures-based method learned from thousands of pure cell types from various sources. xCell applies a novel technique for reducing associations between closely related cell types. 

 **Consensus<sup>TME</sup>** a consensus approach to generating cancer specific signatures for multiple cell types found within the tumour microenvironment.

