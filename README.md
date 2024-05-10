# tigeR
tigeR is an R package designed for the analysis of gene expression in tumor immunotherapy.

### 1.Introduction
- 1060 samples with immunotherapy clinical information from 11 melanoma datasets, including 3 lung cancer datasets, 2 kidney cancer datasets, 1 gastric cancer dataset, 1 low-grade glioma dataset, 1 glioblastoma dataset and 1 Head and Neck Squamous data set (all organized into R language ‘SummarizedExperiment’ objects).

- 23 immunotherapy response related biomarkers from literature, multiple methods for analysis and visualization.

- 10 open source tumor microenvironment deconvolution methods including CIBERSORT, TIMER, ESTIMATE, IPS, xCell, EPIC, ConsensusTME, ABIS, quanTIseq, and MCPCounter. Several downstream method for analysis and visualization.

- 6 machine learning method for multi-modal prediction model construction and testing.

<center>
    <img src="https://raw.githubusercontent.com/Chengxugorilla/tigeR.extra/main/logo.png">
    <p class="caption">**Overall design of tigeR**</p>
</center>

### 2.Installation
```
packages <- c("BiocManager", "devtools", "ggplot2", "pROC")
for (package in packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
}
devtools::install_github("YuLab-SMU/tigeR")
```

### 3.Quick Start
The workflow of tigeR is below, see more details in [tigeR documentation](https://chengxugorilla.github.io/tigeR-book/).

<center>
    <img src="https://raw.githubusercontent.com/Chengxugorilla/tigeR.extra/main/Figure 2.png">
</center>

<p class="caption">Workflow of tigeR</p>

## 4.TIGER web server
http://tiger.canceromics.org/#/
