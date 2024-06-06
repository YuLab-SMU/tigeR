# 🚩 Biomarker Evaluation

## Integrate analysis
 The **integrate_analysis()** function returns the results of both the differential analysis and survival analysis for a gene or gene set within a dataset (or datasets).
```         
integrate_analysis(SE=MEL_GSE91061, geneSet="CD274")
```

## Differential analysis
 You can use **diff_biomk()** to visualize differential analysis result between Pre-Treatment and Post-Treatment patients or Responders and Non-Responders in specified gene.

***Pre-Treatment vs Post-Treatment***

```         
diff_biomk(SE=MEL_GSE91061,gene='CD274',type='Treatment') +
  ggtitle("Pre-Treatment vs Post-Treatment) +
  theme(plot.title = element_text(hjust = 0.5)) 
```

***Responder vs Non-Responder***

```         
diff_biomk(SE=MEL_GSE91061,gene='CD274',type='Response') +
  ggtitle("Responder vs Non-Responder") +
  theme(plot.title = element_text(hjust = 0.5))
```

## Suvival analysis
 You can use **diff_biomk()** to visualize survival analysis result in specified gene.

```         
P <- surv_biomk(SE=MEL_GSE91061,gene='CD274')
P$plot <- P$plot +
  ggtitle("Survival analysis") +
  theme(plot.title = element_text(hjust = 0.5))
P
```

## Calculate comprehensive signature score

 By employing the **score_biomk()** function, you can obtain a comprehensive signature score matrix for the 23 signatures in TigeR.
In this matrix, the columns represent the signature scores, and the rows denote the sample names.

<div style="max-width: 200%; width:780px; height:200px; overflow: scroll">
|Signature|Full Name|Method|PMID|Cancer Type|Target|
|:----:|:-----:|:-------------------------------------:|:------:|:-------:|:----------:|:----------:|
|IRS|immunosenescence-related score|Multivariate cox analysis|35280438|Urothelial Cancer|PD-1/PD-L1|
|tGE8|predefined eight-gene cytotoxic T cell transcriptional signature|Median of Z-score|31686036|Muscle-invasive Urothelial Cancer|PD-L1|
|MEMTS|metastasis-related EMT signature|Average mean|35769483|Gastric Cancer|PD-1|
|PRGScore|pyroptosis-related gene score|Average mean|35479097|Urothelial Cancer; Melanoma|PD-1/PD-L1|
|Angiogenesis|Angiogenesis|Average mean|29867230|Metastatic Renal Cell Carcinoma|PD-L1;VEGF|
|Teffector||Average mean|29867230|Metastatic Renal Cell Carcinoma|PD-L1;VEGF|
|Myeloid_inflammatory||Average mean|29867230|Metastatic Renal Cell Carcinoma|PD-L1;VEGF|
|IFNG_Sig|IFNG-response gene expression signature|Average mean|29150430|Melanoma|CTLA-4|
|TLS|gene signature associated with tertiary lymphoid structures|Weighted mean|31942071|Melanoma|PD-1;CTLA4|
|MSKCC|signature derived from Memorial Sloan Kettering Cancer Center's data|Weighted mean|34421886|Bladder Cancer|PD-1/PD-L1;CTLA4|
|LMRGPI|lipid metabolism-related gene prognostic index|Weighted mean|35582412|Urothelial Cancer|PD-L1|
|PRS|pyroptosis risk score|Weighted mean|35085103|Breast Carcinoma|PD-1;CTLA4|
|Stemnesssignatures|Stemness signatures|Weighted mean|35681225|Colorectal Cancer;Urothelial Cancer;Melanoma|PD-1/PD-L1|
|GRIP|genes related to both inflammation and pyroptosis|Weighted mean|35492358|Melanoma|PD-1; CTLA4|
|IPS|immune prognostic signature|Weighted mean|32572951|Glioblastoma|PD-1|
|Tcell_inflamed_GEP|T cell–inflamedgene expression profile|Weighted mean|30309915|Pan-tumor|PD-1|
|DDR|DNA replication and DNA damage response|Z-score; PCA|29443960|Urothelial Cancer|PD-L1|
|CD8Teffector|CD8+ T-effector|Z-score; PCA|29443960|Non-small Cell Lung Carcinoma|PD-L1|
|CellCycleReg|cell cycle regulator gene sets|Z-score; PCA|29443960|Urothelial Cancer|PD-L1|
|PanFTBRs|pan-fibroblast TGFβ response signature|Z-score; PCA|29443960|Urothelial Cancer|PD-L1|
|EMT1|tumour cell epithelial-to-mesenchymal transition1|Z-score; PCA|29443960|Urothelial Cancer|PD-L1|
|EMT2|tumour cell epithelial-to-mesenchymal transition2|Z-score; PCA|29443960|Urothelial Cancer|PD-L1|
|EMT3|tumour cell epithelial-to-mesenchymal transition3|Z-score; PCA|29443960|Urothelial Cancer|PD-L1|
</div>
   
 In this matrix, the columns represent the signature scores, and the rows denote the sample names.

```
sig_res <- score_biomk(MEL_GSE78220)
```

 Columns represent signatures and rows represent sample.

## Assess the Performance of Signature

 By employing the **roc_biomk()** function, you can assess the performance of built-in and custom signatures in different datasets.
 The function will generate a roc object and a curve to assess the predictive performance.

```
sig_roc <- 
roc_biomk(MEL_PRJEB23709,
          Weighted_mean_Sigs$Tcell_inflamed_GEP,
          method = "Weighted_mean",
          rmBE=TRUE,
          response_NR=TRUE)
```
