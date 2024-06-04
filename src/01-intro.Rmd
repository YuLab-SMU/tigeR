# 📥 Load data into R environment 

## Organization of built-in data
The built-in data follows several specifications for constructing: 
  
① They are constructed into a SummarizedExperiment object using the **SummarizedExperiment()** function in SummarizedExperiment package.  
  
② The assays of the SummarizedExperiment object are stored as a SimpleList containing an expression matrix, where rows represent genes and columns represent samples. 
  
③ The colData of the SummarizedExperiment object is a table(DFrame object of S4), where the row names must be the same as the column names of the expression matrix. The column names of the table represent the clinical information of patients. The detailed information about the column names is presented below.

<div style="width:780px; height:400px; overflow-y: scroll; overflow-x: hidden;">

|          column name          | Recommend value |                                             Description                                              | Necessity in tigeR analysis |
|:-----------:|:-----------:|:---------------------------------:|:-----------:|
|        sample id        |    character    |                                      The sample identity.                                      |              Required              |
|      patient_name       |    character    |                                     The patient identity.                                      |              Alternative           |
|       dataset_id        |    character    |                                     The dataset identity.                                      |              Alternative           |
|        Treatment        |    PRE/POST     |                 The patient sample is collected after Treatment or before Treatment                  |              Required              |
|        response         | CR/MR/PR/PD/SD/NR/R/N |                                     The immunotherapy response.                                |              Alternative           |
|       response_NR       |       R/N       |                        The immunotherapy response which only contains R or N.                        |              Required              |
|         M.stage         |    character    |                                       The M stage of patients.                                       |              Alternative           |
| overall.survival..days. |      value      | The number of days, months or years of survival time of patients.(all samples must be the same unit) |              Required (for survival analysis)              |
|      vital.status       |   Alive/Dead    |                                   The survival status of patient.                                    |              Required (for survival analysis)            |
|     Total.Mutation      |     numeric     |                                   The total mutation gene numbers.                                   |              Alternative           |
|         Gender          |       M/F       |                                        The gender of patient.                                        |              Alternative           |
|         Therapy         |    character    |                                 The anti-tumor therapy on patients.                                  |              Alternative           |
|        age_start        |     numeric     |                                     The age of the patient at diagnosis.                                    |              Alternative           |
|       tumor_type        |    character    |                                          The type of tumor.                                          |              Alternative           |
|        seq_type         |    character    |                                         The sequencing type.                                         |              Alternative           |
|           id            |    character    |                                    The identity of dataset.                                    |              Alternative           |
    
</div>



## Obtain data from tigeR web server
```
Dataloader(c(1,2,3), use_source="Web Server")
```
## Obtain data from ExperimentHub
```
Dataloader(c(1,2,3), use_source="ExperimentHub)
```
## Pre-processing of Custom Data
 When conducting analysis using custom data, you need to pre-process your data and construct a SummarizedExperiment object. Prepare at least one gene expression matrix (rows for genes, columns for samples) and a data frame including corresponding clinical information of the samples. The data frame should include at least the following information:

- sample id
- Treatment (PRE/POST)
- response_NR (R/N)
- overall.survival..days. (for survival analysis)
- vital.status (for survival analysis)

 Here is a brief example:

<div style="width:780px; height:300px; overflow-y: scroll; overflow-x: hidden;">
```
## Construt FPKM Matrix
genes <- c("5S_rRNA", "A1BG", "A1BG-AS1")
samples <- c("SRR3184279", "SRR3184280", "SRR3184281")

mtr <- matrix(
  c(0.04115157, 0.01311394, 0.02275499,
    0.07774626, 0.29494916, 0.11743900,
    2.01782650, 0.83763027, 0.77378394),
  nrow = 3, ncol = 3, byrow = TRUE, 
  dimnames = list(genes, samples))

## Construct Clinical Information Table
data <- data.frame(
  sample_id = c("SRR3184279", "SRR3184280", "SRR3184281"),
  Treatment = c("PRE", "PRE", "PRE"),
  response_NR = c("N", "R", "R"),
  overall_survival_days = c(607, 927, 948),
  vital_status = c("Dead", "Alive", "Alive"),
  stringsAsFactors = FALSE)

## Construct SummarizedExperiment Object
SE <- SummarizedExperiment::SummarizedExperiment(
  assays = list(fpkm = mtr),
  colData = data)

SE
# class: SummarizedExperiment
# dim: 3 3
# metadata(0):
#   assays(1): fpkm
# rownames(3): 5S_rRNA A1BG A1BG-AS1
# rowData names(0):
#   colnames(3): SRR3184279 SRR3184280 SRR3184281
# colData names(5): sample_id Treatment response_NR overall_survival_days
# vital_status
```
</div>

