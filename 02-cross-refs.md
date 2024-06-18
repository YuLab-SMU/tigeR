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

<div style="width:780px; height:200px; overflow-y: scroll; overflow-x: hidden;">
|Signature|Method|PMID|Cancer Type|
|:------:|:------:|:-------:|:----------:|
|IRS|Multivariate cox analysis|35280438|Urothelial Cancer|
|tGE8|Median of Z-score|31686036|Muscle-invasive Urothelial Cancer|
|MEMTS|Average mean|35769483|Gastric Cancer|
|PRGScore|Average mean|35479097|Urothelial Cancer; Melanoma|
|Angiogenesis|Average mean|29867230|Metastatic Renal Cell Carcinoma|
|Teffector|Average mean|29867230|Metastatic Renal Cell Carcinoma|
|Myeloid_inflammatory|Average mean|29867230|Metastatic Renal Cell Carcinoma|
|IFNG_Sig|Average mean|29150430|Melanoma|
|TLS|Weighted mean|31942071|Melanoma|
|MSKCC|Weighted mean|34421886|Bladder Cancer|
|LMRGPI|Weighted mean|35582412|Urothelial Cancer|
|PRS|Weighted mean|35085103|Breast Carcinoma|
|Stemnesssignatures|Weighted mean|35681225|Colorectal Cancer; Urothelial Cancer; Melanoma|
|GRIP|Weighted mean|35492358|PD-1; CTLA4|
|IPS|Weighted mean|32572951|Glioblastoma|
|Tcell_inflamed_GEP|Weighted mean|30309915|Pan-tumor|
|DDR|Z-score; PCA|29443960|Urothelial Cancer|
|CD8Teffector|Z-score; PCA|29443960|Non-small Cell Lung Carcinoma|
|CellCycleReg|Z-score; PCA|29443960|Urothelial Cancer|
|PanFTBRs|Z-score; PCA|29443960|Urothelial Cancer|
|EMT1|Z-score; PCA|29443960|Urothelial Cancer|
|EMT2|Z-score; PCA|29443960|Urothelial Cancer|
|EMT3|Z-score; PCA|29443960|Urothelial Cancer|
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
