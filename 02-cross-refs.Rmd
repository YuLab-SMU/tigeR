# BioMarker Evaluation

## Integrate analysis
```         
integrate_analysis(SE=MEL_GSE91061, geneSet="CD274")
```

## Differential analysis
 You can use diff_biomk() to visualize differential analysis result between Pre-Treatment and Post-Treatment patients or Responders and Non-Responders in specified gene.

***Pre-Treament vs Post-Treatment***

```         
diff_biomk(SE=MEL_GSE91061,gene='CD274',type='Treatment') +
  ggtitle("Treatment vs UnTreatment") +
  theme(plot.title = element_text(hjust = 0.5)) 
```

***Responder vs Non-Responder***

```         
diff_biomk(SE=MEL_GSE91061,gene='CD274',type='Response') +
  ggtitle("Responder vs Non-Responder") +
  theme(plot.title = element_text(hjust = 0.5))
```

## Suvival analysis
 You can use diff_biomk() to visualize survival analysis result in specified gene.

```         
P <- surv_biomk(SE=MEL_GSE91061,gene='CD274')
P$plot <- P$plot +
  ggtitle("Survival analysis") +
  theme(plot.title = element_text(hjust = 0.5))
P
```

## Calculate comprehensive signature score

 By employing the score_biomk() function, you can obtain a comprehensive signature score matrix for the 23 signatures in TigeR.
In this matrix, the columns represent the signature scores, and the rows denote the sample names.

|Signature|Method|Citation|PMID|
|------|------|------|-------|
|IRS|multivariate Cox analysis|Zhou R, Liang J, Tian H, Chen Q, Yang C, Liu C. An Immunosenescence-Related Gene Signature to Evaluate the Prognosis, Immunotherapeutic Response, and Cisplatin Sensitivity of Bladder Cancer. Dis Markers. 2022 Mar 2;2022:2143892. doi: 10.1155/2022/2143892. PMID: 35280438; PMCID: PMC8915927.|35280438|
|tGE8|median of Z-score|Powles T, Kockx M, Rodriguez-Vida A, Duran I, Crabb SJ, Van Der Heijden MS, Szabados B, Pous AF, Gravis G, Herranz UA, Protheroe A, Ravaud A, Maillet D, Mendez MJ, Suarez C, Linch M, Prendergast A, van Dam PJ, Stanoeva D, Daelemans S, Mariathasan S, Tea JS, Mousa K, Banchereau R, Castellano D. Clinical efficacy and biomarker analysis of neoadjuvant atezolizumab in operable urothelial carcinoma in the ABACUS trial. Nat Med. 2019 Nov;25(11):1706-1714. doi: 10.1038/s41591-019-0628-7. Epub 2019 Nov 4. Erratum in: Nat Med. 2020 Jun;26(6):983. Erratum in: Nat Med. 2023 Dec;29(12):3271. PMID: 31686036.|31686036|
|MEMTS|Average mean|Song J, Wei R, Huo S, Gao J, Liu X. Metastasis Related Epithelial-Mesenchymal Transition Signature Predicts Prognosis and Response to Immunotherapy in Gastric Cancer. Front Immunol. 2022 Jun 13;13:920512. doi: 10.3389/fimmu.2022.920512. PMID: 35769483; PMCID: PMC9234207.|35769483|
|PRGScore|Average mean|Zhang Q, Tan Y, Zhang J, Shi Y, Qi J, Zou D, Ci W. Pyroptosis-Related Signature Predicts Prognosis and Immunotherapy Efficacy in Muscle-Invasive Bladder Cancer. Front Immunol. 2022 Apr 11;13:782982. doi: 10.3389/fimmu.2022.782982. PMID: 35479097; PMCID: PMC9035667.|35479097|
|Angiogenesis|Average mean|McDermott DF, Huseni MA, Atkins MB, Motzer RJ, Rini BI, Escudier B, Fong L, Joseph RW, Pal SK, Reeves JA, Sznol M, Hainsworth J, Rathmell WK, Stadler WM, Hutson T, Gore ME, Ravaud A, Bracarda S, Suárez C, Danielli R, Gruenwald V, Choueiri TK, Nickles D, Jhunjhunwala S, Piault-Louis E, Thobhani A, Qiu J, Chen DS, Hegde PS, Schiff C, Fine GD, Powles T. Clinical activity and molecular correlates of response to atezolizumab alone or in combination with bevacizumab versus sunitinib in renal cell carcinoma. Nat Med. 2018 Jun;24(6):749-757. doi: 10.1038/s41591-018-0053-3. Epub 2018 Jun 4. Erratum in: Nat Med. 2018 Dec;24(12):1941. PMID: 29867230; PMCID: PMC6721896.|29867230|
|Teffector|Average mean|McDermott DF, Huseni MA, Atkins MB, Motzer RJ, Rini BI, Escudier B, Fong L, Joseph RW, Pal SK, Reeves JA, Sznol M, Hainsworth J, Rathmell WK, Stadler WM, Hutson T, Gore ME, Ravaud A, Bracarda S, Suárez C, Danielli R, Gruenwald V, Choueiri TK, Nickles D, Jhunjhunwala S, Piault-Louis E, Thobhani A, Qiu J, Chen DS, Hegde PS, Schiff C, Fine GD, Powles T. Clinical activity and molecular correlates of response to atezolizumab alone or in combination with bevacizumab versus sunitinib in renal cell carcinoma. Nat Med. 2018 Jun;24(6):749-757. doi: 10.1038/s41591-018-0053-3. Epub 2018 Jun 4. Erratum in: Nat Med. 2018 Dec;24(12):1941. PMID: 29867230; PMCID: PMC6721896.|29867230|
|Myeloid_inflammatory|Average mean|McDermott DF, Huseni MA, Atkins MB, Motzer RJ, Rini BI, Escudier B, Fong L, Joseph RW, Pal SK, Reeves JA, Sznol M, Hainsworth J, Rathmell WK, Stadler WM, Hutson T, Gore ME, Ravaud A, Bracarda S, Suárez C, Danielli R, Gruenwald V, Choueiri TK, Nickles D, Jhunjhunwala S, Piault-Louis E, Thobhani A, Qiu J, Chen DS, Hegde PS, Schiff C, Fine GD, Powles T. Clinical activity and molecular correlates of response to atezolizumab alone or in combination with bevacizumab versus sunitinib in renal cell carcinoma. Nat Med. 2018 Jun;24(6):749-757. doi: 10.1038/s41591-018-0053-3. Epub 2018 Jun 4. Erratum in: Nat Med. 2018 Dec;24(12):1941. PMID: 29867230; PMCID: PMC6721896.|29867230|
|IFNG_Sig|Average mean|Mo X, Zhang H, Preston S, Martin K, Zhou B, Vadalia N, Gamero AM, Soboloff J, Tempera I, Zaidi MR. Interferon-γ Signaling in Melanocytes and Melanoma Cells Regulates Expression of CTLA-4. Cancer Res. 2018 Jan 15;78(2):436-450. doi: 10.1158/0008-5472.CAN-17-1615. Epub 2017 Nov 17. PMID: 29150430; PMCID: PMC5771950.|29150430|
|TLS|Weighted mean|Cabrita R, Lauss M, Sanna A, Donia M, Skaarup Larsen M, Mitra S, Johansson I, Phung B, Harbst K, Vallon-Christersson J, van Schoiack A, Lövgren K, Warren S, Jirström K, Olsson H, Pietras K, Ingvar C, Isaksson K, Schadendorf D, Schmidt H, Bastholt L, Carneiro A, Wargo JA, Svane IM, Jönsson G. Tertiary lymphoid structures improve immunotherapy and survival in melanoma. Nature. 2020 Jan;577(7791):561-565. doi: 10.1038/s41586-019-1914-8. Epub 2020 Jan 15. Erratum in: Nature. 2020 Apr;580(7801):E1. PMID: 31942071.|31942071|
|MSKCC|Weighted mean|Pan YH, Zhang JX, Chen X, Liu F, Cao JZ, Chen Y, Chen W, Luo JH. Predictive Value of the TP53/PIK3CA/ATM Mutation Classifier for Patients With Bladder Cancer Responding to Immune Checkpoint Inhibitor Therapy. Front Immunol. 2021 Aug 4;12:643282. doi: 10.3389/fimmu.2021.643282. PMID: 34421886; PMCID: PMC8371040.|34421886|
|LMRGPI|Weighted mean|Jiang A, Chen X, Zheng H, Liu N, Ding Q, Li Y, Fan C, Fu X, Liang X, Tian T, Ruan Z, Yao Y. Lipid metabolism-related gene prognostic index (LMRGPI) reveals distinct prognosis and treatment patterns for patients with early-stage pulmonary adenocarcinoma. Int J Med Sci. 2022 Mar 28;19(4):711-728. doi: 10.7150/ijms.71267. PMID: 35582412; PMCID: PMC9108406.|35582412|
|PRS|Weighted mean|Yu H, Fu Y, Tang Z, Jiang L, Qu C, Li H, Tan Z, Shu D, Peng Y, Liu S. A novel pyroptosis-related signature predicts prognosis and response to treatment in breast carcinoma. Aging (Albany NY). 2022 Jan 27;14(2):989-1013. doi: 10.18632/aging.203855. Epub 2022 Jan 27. PMID: 35085103; PMCID: PMC8833126.|35085103|
|Stemnesssignatures|Weighted mean|Zheng H, Liu H, Li H, Dou W, Wang J, Zhang J, Liu T, Wu Y, Liu Y, Wang X. Characterization of stem cell landscape and identification of stemness-relevant prognostic gene signature to aid immunotherapy in colorectal cancer. Stem Cell Res Ther. 2022 Jun 9;13(1):244. doi: 10.1186/s13287-022-02913-0. PMID: 35681225; PMCID: PMC9185878.|35681225|
|GRIP|Weighted mean|Xu Y, Chen Y, Niu Z, Xing J, Yang Z, Yin X, Guo L, Zhang Q, Qiu H, Han Y. A Novel Pyroptotic and Inflammatory Gene Signature Predicts the Prognosis of Cutaneous Melanoma and the Effect of Anticancer Therapies. Front Med (Lausanne). 2022 Apr 15;9:841568. doi: 10.3389/fmed.2022.841568. PMID: 35492358; PMCID: PMC9053829.|35492358|
|IPS|Weighted mean|Zhao B, Wang Y, Wang Y, Chen W, Liu PH, Kong Z, Dai C, Wang Y, Ma W. Systematic identification, development, and validation of prognostic biomarkers involving the tumor-immune microenvironment for glioblastoma. J Cell Physiol. 2021 Jan;236(1):507-522. doi: 10.1002/jcp.29878. Epub 2020 Jun 22. PMID: 32572951.|32572951|
|Tcell_inflamed_GEP|Weighted mean|Cristescu R, Mogg R, Ayers M, Albright A, Murphy E, Yearley J, Sher X, Liu XQ, Lu H, Nebozhyn M, Zhang C, Lunceford JK, Joe A, Cheng J, Webber AL, Ibrahim N, Plimack ER, Ott PA, Seiwert TY, Ribas A, McClanahan TK, Tomassini JE, Loboda A, Kaufman D. Pan-tumor genomic biomarkers for PD-1 checkpoint blockade-based immunotherapy. Science. 2018 Oct 12;362(6411):eaar3593. doi: 10.1126/science.aar3593. Erratum in: Science. 2019 Mar 1;363(6430): PMID: 30309915; PMCID: PMC6718162.|30309915|
|DDR|Z-score;PCA|Mariathasan S, Turley SJ, Nickles D, Castiglioni A, Yuen K, Wang Y, Kadel EE III, Koeppen H, Astarita JL, Cubas R, Jhunjhunwala S, Banchereau R, Yang Y, Guan Y, Chalouni C, Ziai J, Şenbabaoğlu Y, Santoro S, Sheinson D, Hung J, Giltnane JM, Pierce AA, Mesh K, Lianoglou S, Riegler J, Carano RAD, Eriksson P, Höglund M, Somarriba L, Halligan DL, van der Heijden MS, Loriot Y, Rosenberg JE, Fong L, Mellman I, Chen DS, Green M, Derleth C, Fine GD, Hegde PS, Bourgon R, Powles T. TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. Nature. 2018 Feb 22;554(7693):544-548. doi: 10.1038/nature25501. Epub 2018 Feb 14. PMID: 29443960; PMCID: PMC6028240.|29443960|
|CD8Teffector|Z-score;PCA|Mariathasan S, Turley SJ, Nickles D, Castiglioni A, Yuen K, Wang Y, Kadel EE III, Koeppen H, Astarita JL, Cubas R, Jhunjhunwala S, Banchereau R, Yang Y, Guan Y, Chalouni C, Ziai J, Şenbabaoğlu Y, Santoro S, Sheinson D, Hung J, Giltnane JM, Pierce AA, Mesh K, Lianoglou S, Riegler J, Carano RAD, Eriksson P, Höglund M, Somarriba L, Halligan DL, van der Heijden MS, Loriot Y, Rosenberg JE, Fong L, Mellman I, Chen DS, Green M, Derleth C, Fine GD, Hegde PS, Bourgon R, Powles T. TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. Nature. 2018 Feb 22;554(7693):544-548. doi: 10.1038/nature25501. Epub 2018 Feb 14. PMID: 29443960; PMCID: PMC6028240.|29443960|
|CellCycleReg|Z-score;PCA|Mariathasan S, Turley SJ, Nickles D, Castiglioni A, Yuen K, Wang Y, Kadel EE III, Koeppen H, Astarita JL, Cubas R, Jhunjhunwala S, Banchereau R, Yang Y, Guan Y, Chalouni C, Ziai J, Şenbabaoğlu Y, Santoro S, Sheinson D, Hung J, Giltnane JM, Pierce AA, Mesh K, Lianoglou S, Riegler J, Carano RAD, Eriksson P, Höglund M, Somarriba L, Halligan DL, van der Heijden MS, Loriot Y, Rosenberg JE, Fong L, Mellman I, Chen DS, Green M, Derleth C, Fine GD, Hegde PS, Bourgon R, Powles T. TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. Nature. 2018 Feb 22;554(7693):544-548. doi: 10.1038/nature25501. Epub 2018 Feb 14. PMID: 29443960; PMCID: PMC6028240.|29443960|
|PanFTBRs|Z-score;PCA|Mariathasan S, Turley SJ, Nickles D, Castiglioni A, Yuen K, Wang Y, Kadel EE III, Koeppen H, Astarita JL, Cubas R, Jhunjhunwala S, Banchereau R, Yang Y, Guan Y, Chalouni C, Ziai J, Şenbabaoğlu Y, Santoro S, Sheinson D, Hung J, Giltnane JM, Pierce AA, Mesh K, Lianoglou S, Riegler J, Carano RAD, Eriksson P, Höglund M, Somarriba L, Halligan DL, van der Heijden MS, Loriot Y, Rosenberg JE, Fong L, Mellman I, Chen DS, Green M, Derleth C, Fine GD, Hegde PS, Bourgon R, Powles T. TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. Nature. 2018 Feb 22;554(7693):544-548. doi: 10.1038/nature25501. Epub 2018 Feb 14. PMID: 29443960; PMCID: PMC6028240.|29443960|
|EMT1|Z-score;PCA|Mariathasan S, Turley SJ, Nickles D, Castiglioni A, Yuen K, Wang Y, Kadel EE III, Koeppen H, Astarita JL, Cubas R, Jhunjhunwala S, Banchereau R, Yang Y, Guan Y, Chalouni C, Ziai J, Şenbabaoğlu Y, Santoro S, Sheinson D, Hung J, Giltnane JM, Pierce AA, Mesh K, Lianoglou S, Riegler J, Carano RAD, Eriksson P, Höglund M, Somarriba L, Halligan DL, van der Heijden MS, Loriot Y, Rosenberg JE, Fong L, Mellman I, Chen DS, Green M, Derleth C, Fine GD, Hegde PS, Bourgon R, Powles T. TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. Nature. 2018 Feb 22;554(7693):544-548. doi: 10.1038/nature25501. Epub 2018 Feb 14. PMID: 29443960; PMCID: PMC6028240.|29443960|
|EMT2|Z-score;PCA|Mariathasan S, Turley SJ, Nickles D, Castiglioni A, Yuen K, Wang Y, Kadel EE III, Koeppen H, Astarita JL, Cubas R, Jhunjhunwala S, Banchereau R, Yang Y, Guan Y, Chalouni C, Ziai J, Şenbabaoğlu Y, Santoro S, Sheinson D, Hung J, Giltnane JM, Pierce AA, Mesh K, Lianoglou S, Riegler J, Carano RAD, Eriksson P, Höglund M, Somarriba L, Halligan DL, van der Heijden MS, Loriot Y, Rosenberg JE, Fong L, Mellman I, Chen DS, Green M, Derleth C, Fine GD, Hegde PS, Bourgon R, Powles T. TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. Nature. 2018 Feb 22;554(7693):544-548. doi: 10.1038/nature25501. Epub 2018 Feb 14. PMID: 29443960; PMCID: PMC6028240.|29443960|
|EMT3|Z-score;PCA|Mariathasan S, Turley SJ, Nickles D, Castiglioni A, Yuen K, Wang Y, Kadel EE III, Koeppen H, Astarita JL, Cubas R, Jhunjhunwala S, Banchereau R, Yang Y, Guan Y, Chalouni C, Ziai J, Şenbabaoğlu Y, Santoro S, Sheinson D, Hung J, Giltnane JM, Pierce AA, Mesh K, Lianoglou S, Riegler J, Carano RAD, Eriksson P, Höglund M, Somarriba L, Halligan DL, van der Heijden MS, Loriot Y, Rosenberg JE, Fong L, Mellman I, Chen DS, Green M, Derleth C, Fine GD, Hegde PS, Bourgon R, Powles T. TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. Nature. 2018 Feb 22;554(7693):544-548. doi: 10.1038/nature25501. Epub 2018 Feb 14. PMID: 29443960; PMCID: PMC6028240.|29443960|
### algorithm

 The Signature scores are mostly represented by the average mean or weighted mean or ZScore of the expression value of specified genes.\
 For example, MSKCC Signature are calculated by the following formula:\
 \
  $MSKCC=-0.492\times exp(TP53)+0.562\times exp(PIK3CA)+1.454\times exp(ATM)$\
 \
 where exp() represents the expression value of the gene.

```
score_biomk(MEL_GSE78220)
```

 Columns represent Signatures and rows represent sample.

## Assess the performance of Signature

 By employing the auc_biomk() function, you can assess the performance of Signature(including user-built Signature) in different datasets.
The function will return a "roc"" object, a list of class "roc".

```
roc_biomk(MEL_PRJEB23709,
          Weighted_mean_Sigs$Tcell_inflamed_GEP,
          method = "Weighted_mean",
          rmBE=TRUE,
          response_NR=TRUE)
```
