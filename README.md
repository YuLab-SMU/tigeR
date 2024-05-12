# tigeR
tigeR is an R package designed for the analysis of gene expression in tumor immunotherapy.

### 1. Introduction
- Built-in datasets: 1060 samples with immunotherapy clinical information from 11 melanoma datasets, 3 lung cancer datasets, 2 kidney cancer datasets, 1 gastric cancer dataset, 1 low-grade glioma dataset, 1 glioblastoma dataset and 1 head and neck squamous cell cancer dataset (all organized into R language ‘SummarizedExperiment’ objects).

- 23 immunotherapy response-related biomarkers from literature, multiple methods for analysis and visualization.

- 10 open source tumor microenvironment deconvolution methods including CIBERSORT, TIMER, ESTIMATE, IPS, xCell, EPIC, ConsensusTME, ABIS, quanTIseq, and MCPCounter. Several downstream method for analysis and visualization.

- 7 machine learning method for multi-modal prediction model construction and testing.

<p align="center">
<img src="https://raw.githubusercontent.com/Chengxugorilla/tigeR.extra/main/Figure 1.png" width="80%">
</p>
<p align="center"><b>Overall design of tigeR</b></p>

### 2. Installation
```
packages <- c("BiocManager", "devtools", "ggplot2", "pROC")
for (package in packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
}
devtools::install_github("YuLab-SMU/tigeR")
```

### 3. Quick Start
The workflow of tigeR is below, see more details in [tigeR documentation](https://chengxugorilla.github.io/tigeR-book/).

<p align="center">
<img src="https://raw.githubusercontent.com/Chengxugorilla/tigeR.extra/main/Figure 2.svg" alt="SVG Image" width="80%">
</p>
<p align="center"><b>Workflow of tigeR</b></p>

### 4. TIGER web server
[TIGER Web Server](http://tiger.canceromics.org/#/)
