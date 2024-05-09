# 1 Load built-in data into R environment

All chapters start with a first-level heading followed by your chapter title, like the line above. There should be only one first-level heading (`#`) per .Rmd file.

## 1.1 Obtain data from data package(recommandate)
```
devtools::install.github("YuLab-SMU/tigeR")
```
## 1.2 Obtain data from tigeR web server
```
Dataloader(c(1,2,3), use_source="Web Server")
```
## 1.3 Obtain data from ExperimentHub
```
Dataloader(c(1,2,3), use_source="ExperimentHub)
```
