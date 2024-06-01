# 🚩 Biomarker Evaluation

## Integrate analysis

 The **integrate_analysis()** function returns the results of both the differential analysis and survival analysis for a gene or gene set within a dataset (or datasets).


```r
tigeR::integrate_analysis(SE=MEL_GSE91061, geneSet="CD274")
```

```
## $`Response vs Non-Response`
##    log2(FC)         P     Score
## 1 0.3162897 0.5257107 0.2792532
## 
## $`Pre-Therapy vs Post-Therapy`
##     log2(FC)          P     Score
## 1 -0.8121796 0.01780415 -1.749479
## 
## $Survival
##         HR          P      Score 
##  0.9203840  0.8175588 -0.0874810
```
 The score of differential expression analysis is calculated by following formula:

  $-SIGN(log2(FC)) \times log10(p)$

 where $FC$ represents the fold change and $p$ represents the P value derived from the Wilcoxon rank-sum test

 The score of survival analysis is calculated by the following formula:

  $-SIGN(log2(HR)) \times log10(p)$

 where $HR$ represents the hazard ratio and $p$ represents the P value derived from univariate Cox regression analysis.

## Differential analysis
 You can use **diff_biomk()** to visualize differential analysis result between Pre-Treatment and Post-Treatment patients or Responders and Non-Responders in specified gene.

***Pre-Treatment vs Post-Treatment***


```r
diff_biomk(SE=MEL_GSE91061,gene='CD274',type='Treatment') +
  ggtitle("Pre-Treatment vs Post-Treatment") +
  theme(plot.title = element_text(hjust = 0.5)) 
```

<img src="02-cross-refs_files/figure-html/unnamed-chunk-3-1.png" width="672" />

***Responder vs Non-Responder***


```r
diff_biomk(SE=MEL_GSE91061,gene='CD274',type='Response') +
  ggtitle("Responder vs Non-Responder") +
  theme(plot.title = element_text(hjust = 0.5))
```

<img src="02-cross-refs_files/figure-html/unnamed-chunk-4-1.png" width="672" />

## Suvival analysis
 You can use **diff_biomk()** to visualize survival analysis result in specified gene.


```r
P <- surv_biomk(SE=MEL_GSE91061,gene='CD274')
P$plot <- P$plot +
  ggtitle("Survival analysis") +
  theme(plot.title = element_text(hjust = 0.5))
P
```

<img src="02-cross-refs_files/figure-html/unnamed-chunk-5-1.png" width="672" />

## Calculate comprehensive signature score

 By employing the **score_biomk()** function, you can obtain a comprehensive signature score matrix for the 23 signatures in TigeR.

<div style="width:780px; height:200px; overflow-y: scroll; overflow-x: hidden;">
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
</div>
   
 In this matrix, the columns represent the signature scores, and the rows denote the sample names.

```r
score_biomk(MEL_GSE78220)
```

```
## 2 Signature genes are not found in expression matrix. The function can execute properly, but the performance of the model may be compromised.
## 2 Signature genes are not found in expression matrix. The function can execute properly, but the performance of the model may be compromised.
```

```
## 1 Signature genes are not found in expression matrix. The function can execute properly, but the performance of the model may be compromised.
```

```
## 3 Signature genes are not found in expression matrix. The function can execute properly, but the performance of the model may be compromised.
```

```
## 1 Signature genes are not found in expression matrix. The function can execute properly, but the performance of the model may be compromised.
```

```
##               IRS          tGE8     MEMTS  PRGScore Angiogenesis  Teffector
## SRR3184279  1.142 -0.1274466254 11.940531 11.771754    13.484165  3.1620737
## SRR3184280  0.563 -0.2455777711  2.942161 15.670362     4.718858  2.4850581
## SRR3184281  0.002 -0.2939281368  7.624309  6.816627     8.210066  2.3183895
## SRR3184282  0.673 -0.6286860863 12.379746 12.612995     3.466574  1.3979274
## SRR3184283  1.142 -0.5815497283  7.371334  6.633541     5.090331  0.6548324
## SRR3184284  0.178 -0.5199876731  7.142007  5.804861     4.044871  2.1538591
## SRR3184285  1.142 -0.4435026383  4.701351 13.242594     6.753375  1.4941546
## SRR3184286  1.142  0.0247948123 11.944323 16.733508     5.829289  3.2385169
## SRR3184287  0.483 -0.6093308999 13.386835  4.038779    12.556777  0.2898800
## SRR3184288  0.978 -0.5721122027 12.341645  7.804835    12.991586  0.4104385
## SRR3184289 -0.153  0.3179526130  9.868764 16.361679    16.475330  5.8068096
## SRR3184290  1.142 -0.2112894772 17.120574  7.901416    18.544065  5.9687710
## SRR3184291  0.647 -0.2246103754 11.037060 12.285773     4.101229  2.3418387
## SRR3184293  0.178 -0.4445736563  3.502348  7.806506     1.116110  1.3780202
## SRR3184294 -0.927  4.2997061585  8.799417 37.257834     4.160336 24.3030619
## SRR3184295  0.766  0.2893208324 11.202873 14.826802    17.870372  4.2952791
## SRR3184296  0.766  0.0922226574  8.417918 12.755647    13.464519  4.5893147
## SRR3184297  0.068  0.0006039458 12.279312 17.480623     3.324207  5.2497422
## SRR3184298 -0.282  0.1014041880  5.494192 19.969046     4.096930  3.9572014
## SRR3184299  0.094 -0.3306502913  5.900445 13.145733     1.934650  1.5245628
## SRR3184300  1.365 -0.5966846648  3.025450 10.208271    14.647850  1.5278598
## SRR3184301  0.673  0.1562774235  5.668299 13.787251     3.984429  4.4077161
## SRR3184302  0.952 -0.4878345506  5.036720  8.087117     6.113358  2.4108979
## SRR3184303  1.175 -0.5120525308 12.718229  5.550579     8.257533  0.6931763
## SRR3184304  1.142 -0.2146564598 13.498657 12.107371     6.090359  2.6757308
## SRR3184305  0.647 -0.3030426880  4.506788 15.700597     1.185254  2.4191943
## SRR3184306  0.187  0.9173511615  4.878808 26.143796     4.544727  8.8775902
##            Myeloid_inflammatory   IFNG_Sig        TLS      MSKCC    LMRGPI
## SRR3184279            1.0387236  53.557553  12.915505  -7.359213 2228.2565
## SRR3184280           12.5713176  58.944631  10.167293  -5.301346 1995.6762
## SRR3184281            0.2850853  13.038679   4.112142  -6.058954  389.2787
## SRR3184282            0.2398302   7.214410   1.189567  -9.994039  711.2312
## SRR3184283            1.4416500  13.734622  25.598985  -8.526128 1227.0611
## SRR3184284            0.9009998   7.478617   6.002477  -7.427001 1108.5554
## SRR3184285           14.4184585  29.185365   4.655421  -6.597305 2060.4087
## SRR3184286            1.3936457  22.239602   7.665335 -12.556601 1570.2545
## SRR3184287           10.1350623  15.243892  26.373111  -9.463943 3115.7273
## SRR3184288            1.1661599  12.086776   2.749075  -7.165436 6622.5043
## SRR3184289            3.0921195 140.828609  19.711734  -9.099641 4762.2341
## SRR3184290           10.2093485  17.081297   6.390681 -10.891821 4260.5843
## SRR3184291           11.4705920  23.884921   5.595248  -6.737246 2346.6850
## SRR3184293            1.0430889  25.300356   7.496732  -7.901521  732.6875
## SRR3184294            1.5226401  71.407357  23.872684  -6.593064  448.9391
## SRR3184295            3.6507347  82.897006  15.510554  -8.787397 1289.2336
## SRR3184296            9.4160121  30.453679   8.590768  -6.104721 1495.6680
## SRR3184297            0.1214514  51.329264   6.230481  -9.162395  456.2589
## SRR3184298            0.6266078  55.996400   6.512297 -11.382099 1195.3969
## SRR3184299            0.2881308  33.255712   3.458758 -12.201856  626.8247
## SRR3184300           10.5246139  10.526540  11.048209  -5.439693 6965.1493
## SRR3184301            1.0215001  33.060018  16.997479 -10.657366  317.0471
## SRR3184302            0.9821260   9.197015   3.462501  -4.161005 1602.2385
## SRR3184303            2.0218686   9.692717   2.326838 -11.425945 1955.7707
## SRR3184304            3.0919466  21.606968   4.344373  -8.418186 1146.6345
## SRR3184305            0.4213449  13.765030   2.669343  -5.821114  506.0899
## SRR3184306            1.5901618  65.185269 352.070341 -17.961746  776.7463
##                  PRS Stemnesssignatures      GRIP       IPS Tcell_inflamed_GEP
## SRR3184279  65.64614          6.0551307 0.7086824 11.968088          2.5122324
## SRR3184280  48.09594          0.6347839 0.3318125  7.610523          1.8561422
## SRR3184281  51.95532         10.8826249 0.1079182 40.114020          1.1043452
## SRR3184282  53.21504          0.5174279 0.7894171  1.216898          1.4095208
## SRR3184283  35.56187          3.1279351 0.2258769  5.557282          0.9858985
## SRR3184284  43.48133          1.2896477 0.5710363  3.358411          0.6684792
## SRR3184285  48.26292          3.3684499 0.3071354 28.515524          0.9053413
## SRR3184286  56.22087          9.1048007 0.5117116 12.905714          1.9154800
## SRR3184287  59.74546          4.9152321 0.4100044 67.712721          0.6420798
## SRR3184288  56.08383          2.8460659 0.4089724 24.906047          0.9800455
## SRR3184289  78.14104          6.5724987 0.4878133 21.561503          2.6782935
## SRR3184290  91.18922          6.0947345 0.4023069 34.593054          1.5020027
## SRR3184291  47.85892          1.5540122 0.3056658 18.890261          1.3345470
## SRR3184293  29.13646          2.0998963 1.0319277 11.190063          0.9089890
## SRR3184294  84.70110          2.0782220 0.5243760 14.328041          9.7224509
## SRR3184295  78.40198          6.9836082 0.4056731 17.217759          3.0515033
## SRR3184296  58.83162          2.5230547 0.3465551 15.716402          1.9470141
## SRR3184297  60.97359          3.9215207 0.3233180 30.681658          2.3443459
## SRR3184298  70.20240          1.2193907 0.3635129  5.833457          3.5044841
## SRR3184299  52.61078          0.6038819 0.2964016  2.655690          1.9560381
## SRR3184300 238.46248          0.7878857 0.1440482  6.586125          0.4794043
## SRR3184301  57.80541          2.2348079 0.1844206 10.864492          2.2706593
## SRR3184302  36.28359          3.3680218 0.6856813  9.693200          1.2748896
## SRR3184303  57.34977         17.4874321 0.3597441 27.285943          0.9915240
## SRR3184304  63.41676          2.6713151 0.4259503  7.565359          1.7037176
## SRR3184305  30.89488          1.0163282 0.5628645  3.387794          1.2979521
## SRR3184306  94.58496          2.8473755 0.2749565 48.883193          4.6277454
##                   DDR CD8Teffector CellCycleReg    PanFTBRs       EMT1
## SRR3184279 -0.1975844   0.18568449    0.1893995 -0.21084908 0.19653441
## SRR3184280 -0.1915641   0.22239930    0.2128977 -0.21192319 0.19491950
## SRR3184281 -0.1949043   0.22718870    0.2093351 -0.21394746 0.19079958
## SRR3184282 -0.1979084   0.15267379    0.2102159 -0.20981667 0.19649051
## SRR3184283 -0.1828094   0.01990807    0.2182126 -0.20434785 0.19273395
## SRR3184284 -0.1899734   0.22455493    0.1867713 -0.19124966 0.19548481
## SRR3184285 -0.1928915   0.22023056    0.2148923 -0.21580511 0.19318426
## SRR3184286 -0.1976339   0.21349610    0.1983098 -0.18235451 0.19645178
## SRR3184287 -0.1895723   0.12582645    0.1797330 -0.07789503 0.19536573
## SRR3184288 -0.1854981   0.08131488    0.1899819 -0.20374203 0.19557528
## SRR3184289 -0.2011675   0.22262781    0.1449800 -0.04490327 0.19486434
## SRR3184290 -0.1908081   0.17288617    0.1616639 -0.11617947 0.19558674
## SRR3184291 -0.1932547   0.22474075    0.2144560 -0.20991059 0.19559778
## SRR3184293 -0.1656256   0.17609773    0.2163997 -0.21410889 0.19608374
## SRR3184294 -0.1978316   0.21934150    0.1706078 -0.21464499 0.19575590
## SRR3184295 -0.2012071   0.22251676    0.1949373 -0.12477122 0.19602983
## SRR3184296 -0.1894981   0.18912942    0.1524055 -0.19277433 0.19567387
## SRR3184297 -0.1930148   0.20284639    0.1964675 -0.21293053 0.19567933
## SRR3184298 -0.1960135   0.22265583    0.1935689 -0.20042117 0.19641906
## SRR3184299 -0.1947792   0.21965300    0.1918343 -0.19794594 0.19645731
## SRR3184300 -0.1965829   0.20929758    0.1176218 -0.21420376 0.09200499
## SRR3184301 -0.1964928   0.22206368    0.2105635 -0.21445934 0.19569386
## SRR3184302 -0.1848634   0.15431025    0.1793152 -0.21428410 0.19609495
## SRR3184303 -0.1969067   0.09703136    0.2147835 -0.13999062 0.19409875
## SRR3184304 -0.1955291   0.22397810    0.1701928 -0.21493569 0.19397346
## SRR3184305 -0.1866733   0.18566219    0.2174461 -0.21262036 0.19571413
## SRR3184306 -0.1921057   0.17167210    0.1962221 -0.19336206 0.19601393
##                   EMT2       EMT3
## SRR3184279 -0.21999489 0.22368777
## SRR3184280 -0.19427188 0.25123034
## SRR3184281 -0.22460166 0.04031939
## SRR3184282 -0.17621775 0.12836018
## SRR3184283 -0.14874344 0.23743344
## SRR3184284 -0.15836417 0.19684673
## SRR3184285 -0.21875967 0.25953850
## SRR3184286 -0.22425060 0.25674730
## SRR3184287 -0.13749826 0.19659773
## SRR3184288 -0.15804581 0.08880141
## SRR3184289 -0.23235635 0.24175944
## SRR3184290 -0.14660386 0.06552324
## SRR3184291 -0.21928083 0.22340912
## SRR3184293 -0.22048507 0.21131955
## SRR3184294 -0.16567045 0.06541941
## SRR3184295 -0.21647219 0.17629116
## SRR3184296 -0.23015242 0.19765431
## SRR3184297 -0.05907203 0.18270602
## SRR3184298 -0.14745456 0.24405227
## SRR3184299 -0.08584916 0.19203458
## SRR3184300 -0.20418504 0.03081705
## SRR3184301 -0.22348333 0.23091879
## SRR3184302 -0.22320867 0.25571116
## SRR3184303 -0.14165774 0.11139717
## SRR3184304 -0.23293684 0.21767458
## SRR3184305 -0.22391943 0.22951677
## SRR3184306 -0.21247507 0.01099541
```

## Assess the Performance of Signature

 By employing the **roc_biomk()** function, you can assess the performance of built-in and custom signatures in different datasets.
 The function will generate a roc object and a curve to assess the predictive performance.


```r
roc_biomk(MEL_PRJEB23709,
          Weighted_mean_Sigs$Tcell_inflamed_GEP,
          method = "Weighted_mean",
          rmBE=TRUE,
          response_NR=TRUE)
```

```
## 3 Signature genes are not found in expression matrix. The function can execute properly, but the performance of the model may be compromised.
```

```
##  HLA.DQA1 HLA.DRB1 HLA.E does not exist in expression matrix.
```

```
## Setting levels: control = N, case = R
```

```
## Setting direction: controls < cases
```

```
## [[1]]
## 
## Call:
## roc.default(response = data[[2]]$response, predictor = value)
## 
## Data: value in 33 controls (data[[2]]$response N) < 40 cases (data[[2]]$response R).
## Area under the curve: 0.8364
## 
## [[2]]
```

<img src="02-cross-refs_files/figure-html/unnamed-chunk-7-1.png" width="672" />
