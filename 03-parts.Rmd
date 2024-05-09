# Tumor Microenviroment Analysis

## TIMER
&emsp;TIMER is a comprehensive tool for systematical analysis of immune infiltrates across diverse cancer types.
&emsp;**TIMER()** function will return a cell type relative abundance matrix.

```
TIMER(assay(MEL_GSE91061),type = "SKCM")
```

## CIBERSORT
 CIBERSORT is an analytical tool from the Alizadeh Lab and Newman Lab to impute gene expression profiles and provide an estimation of the abundances of member cell types in a mixed cell population, using gene expression data.
 tigeR offer you an built-in function **CIBERSORT()** for estimation of the abundances of member cell types.
**CIBERSORT()** function will return a list which first element is the cell type proportion matrix and second element is a boxplot.

```         
result <- CIBERSORT(sig_matrix = LM22, SE = MEL_GSE78220, perm = 0, QN = TRUE)
```

## xCell
 xCell is a gene signatures-based method learned from thousands of pure cell types from various sources. xCell applies a novel technique for reducing associations between closely related cell types. 
 tigeR offer you an built-in function **xCell()** for estimation of the abundances of member cell types.
**CIBERSORT()** function will return a 

```         
xCell(MEL_GSE91061)
```

## ConsensusTME
 Consensus<sup>TME</sup> a consensus approach to generating cancer specific signatures for multiple cell types found within the tumour microenvironment.

```
ConsensusTME(MEL_GSE91061,cancer = "SKCM")
```
