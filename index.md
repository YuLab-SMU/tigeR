--- 
title: "tigeR: Tumor Immunotherapy Gene Expression Data Analysis R package"
author: "Yihao Chen, Li-Na He, Yuanzhe Zhang, Yuelong Shu, Di Zhang, Guangchuang Yu, Zhixiang Zuo"
date: "2024-06-13"
site: bookdown::bookdown_site
output: bookdown::gitbook
documentclass: book
bibliography: [book.bib, packages.bib]
# url: your book url like https://bookdown.org/yihui/bookdown
# cover-image: path to the social sharing image like images/cover.jpg
description:
link-citations: yes
github-repo: rstudio/bookdown-demo
---

# 📖 **Getting start with tigeR**{-}

## 🔬 Introduction
<p align="center">
<img src="./figs/Figure 1.svg" alt="SVG Image">
</p>

 tigeR encompasses **4** distinct yet closely interconnected modules:

- The **Biomarker Evaluation module** enables researchers to evaluate whether the biomarkers of interest are associated with immunotherapy response via built-in or custom immunotherapy gene expression data. 
- The **Tumor Microenvironment Deconvolution module** integrates 10 open-source algorithms to obtain the proportions of different cell types within the tumor microenvironment, facilitating the investigation of the association between immune cell populations and immunotherapy response. 
- The **Prediction Model Construction module** equips users with the ability to construct sophisticated prediction models using a range of built-in machine learning algorithms.
- The **Response Prediction module** predict the immunotherapy response for the patients from gene expression data using our pre-trained machine learning models or public gene expression signatures.

## 🏞 The workflow of tigeR
<p align="center">
<img src="./figs/Figure 2.svg" alt="SVG Image">
</p>
 As depicted in this figure, users have the flexibility to load built-in gene expression data with immunotherapy outcome information or to utilize their own data for subsequent analysis. 

 The Biomarker Evaluation module serves to assess the correlation between biomarkers and immunotherapy outcomes. 

 Furthermore, the Tumor Microenvironment Deconvolution module enables the derivation of cell type proportions within the tumor microenvironment using 10 open-source algorithms. This module also provides functionality for evaluating the interplay between fractions of tumor microenvironment cells and immunotherapy outcomes. 

 Subsequently, based on the features selected from these two modules, users can leverage the Prediction Model Construction module, which incorporates a range of machine learning models, to train a model for predicting immunotherapy response using transcriptome gene expression data. 

 Users can use the Response Prediction module to predict the immunotherapy response for the patients from gene expression data using our pre-trained machine learning models or public gene expression signatures.

## 🛠️ Installation

```
packages <- c("BiocManager", "devtools", "ggplot2", "pROC", "RobustRankAggreg")
for (package in packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
}
devtools::install_github("YuLab-SMU/tigeR")
```

## 🗺️ Function Overview

**Data Loading**

Dataloader(): load data online (ExperimentHub or TIGER Web Server).

**Biomarkder Evaluation**

score_biomk(): generate a comprehensive signature score matrix for the 23 signatures in tigeR. Columns represent the signature scores, while rows denote the sample names.

score_biomk_SE(): generate a comprehensive signature score matrix for the 23 signatures in tigeR. Columns represent the signature scores, while rows denote the sample names. (return a SummarizedExperiment object)

roc_biomk(): generate a Receiver Operating Characteristic (ROC) object and a curve to assess the predictive performance.

diff_biomk(): plot differential result (Responder vs Non-Responder or Pre-Treatment vs Post-Treatment).

surv_biomk(): calculates hazard ratios, confidence intervals and P value of cox-ph analysis as well as draw KM curve.

compare_biomk(): generate a heatmap of signature auc of datasets.

browser_biomk(): generate an integration diagram comprising a bar plot representing AUC and a dot plot denoting Hazard Ratio and P-value.

integrate_analysis(): perform differential expression analysis and survival analysis.

**TME Deconvolution**

TME_deconvolution(): perform Tumor Microenvironment deconvolution through 10 open-source algorithms.

fraction_pie(): generate a pie plot illustrating the cell fraction or relative cell abundance for each sample.

**Prediction Model Construction**

data_standardization(): perform data standardization, including converting FPKM to TPM, removing NA values, applying log transformation, and scaling the data.

compare_roc(): plot all the roc curves on a single plot.

gini_gene(): calculating the Gini index and get an overview of the classification efficiency of genes.

diff_gene(): return differential expression gene between Responder and Non-Responder.

build_Model(): generate immunotherapy response prediction machine-learning model.

test_Model(): test model generated by build_Model function.

**Response Prediction**

pred_response(): predict immunotherapy response and generate heatmap of signatures.

## 💥 Troubleshooting
Please report bugs in [Github issue](https://github.com/YuLab-SMU/tigeR/issues).



