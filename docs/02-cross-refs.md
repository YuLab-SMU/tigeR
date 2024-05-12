# 🚩 Biomarker Evaluation

## Integrate analysis
 The **integrate_analysis()** function returns the results of both the differential analysis and survival analysis for a gene or gene set within a dataset (or datasets).
```         
integrate_analysis(SE=MEL_GSE91061, geneSet="CD274")
```

## Differential analysis
 You can use **diff_biomk()** to visualize differential analysis result between Pre-Treatment and Post-Treatment patients or Responders and Non-Responders in specified gene.

***Pre-Treament vs Post-Treatment***

```         
diff_biomk(SE=MEL_GSE91061,gene='CD274',type='Treatment') +
  ggtitle("Pre-Treament vs Post-Treatment) +
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

|Signature|Method|PMID|
|------|------|-------|
|IRS|multivariate Cox analysis|35280438|
|tGE8|median of Z-score|31686036|
|MEMTS|Average mean|35769483|
|PRGScore|Average mean|35479097|
|Angiogenesis|Average mean|29867230|
|Teffector|Average mean|29867230|
|Myeloid_inflammatory|Average mean|29867230|
|IFNG_Sig|Average mean|29150430|
|TLS|Weighted mean|31942071|
|MSKCC|Weighted mean|34421886|
|LMRGPI|Weighted mean|35582412|
|PRS|Weighted mean|35085103|
|Stemnesssignatures|Weighted mean|35681225|
|GRIP|Weighted mean|35492358|
|IPS|Weighted mean|32572951|
|Tcell_inflamed_GEP|Weighted mean|30309915|
|DDR|Z-score;PCA|29443960|
|CD8Teffector|Z-score;PCA|29443960|
|CellCycleReg|Z-score;PCA|29443960|
|PanFTBRs|Z-score;PCA|29443960|
|EMT1|Z-score;PCA|29443960|
|EMT2|Z-score;PCA|29443960|
|EMT3|Z-score;PCA|29443960|

```
score_biomk(MEL_GSE78220)
```

 Columns represent signatures and rows represent sample.

## Assess the Performance of Signature

 By employing the **roc_biomk()** function, you can assess the performance of built-in and custom signatures in different datasets.
 The function will generate a roc object and a curve to assess the predictive performance.

```
roc_biomk(MEL_PRJEB23709,
          Weighted_mean_Sigs$Tcell_inflamed_GEP,
          method = "Weighted_mean",
          rmBE=TRUE,
          response_NR=TRUE)
```
