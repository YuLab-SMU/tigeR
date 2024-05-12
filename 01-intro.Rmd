# 📥 Load data into R environment 

## Organization of built-in data
The built-in data follows several specifications for constructing: 
  
① They are constructed into a SummarizedExperiment object using the **SummarizedExperiment()** function in SummarizedExperiment package.  
  
② The assays of the SummarizedExperiment object are stored as a SimpleList containing an expression matrix, where rows represent genes and columns represent patients. 
  
③ The colData of the SummarizedExperiment object is a table(DFrame object of S4), where the row names must be the same as the column names of the expression matrix. The column names of the table represent the clinical information of patients. The detailed information about the column names is presented below.

<div style="width:780px; height:400px; overflow-y: scroll; overflow-x: hidden;">

|          column name          | Recommend value |                                             Description                                              | Necessity in tigeR analysis |
|:-----------:|:-----------:|:---------------------------------:|:-----------:|
|        sample id        |    character    |                                      The sample identification.                                      |              ✕              |
|      patient_name       |    character    |                                     The patient identification.                                      |              ✕              |
|       dataset_id        |    character    |                                     The dataset identification.                                      |              ✕              |
|      dataset_group      |    character    |                                                                                                      |              ✕              |
|        Treatment        |    PRE/POST     |                 The patient sample is collected after Treatment or before Treatment                  |              ✓              |
|        response         | CR/MR/PR/PD/SD/NR/R/N |                                     The Immunotherapy response.                                      |              ✕              |
|       response_NR       |       R/N       |                        The Immunotherapy response which only contains R or N.                        |              ✓              |
|         M.stage         |    character    |                                       The M stage of patients.                                       |              ✕              |
| overall.survival..days. |      value      | The number of days, months or years of survival time of patients.(all samples must be the same unit) |              ✓              |
|      vital.status       |   Alive/Dead    |                                   The survival status of patient.                                    |              ✓              |
|     Total.Mutation      |     numeric     |                                   The total mutation gene numbers.                                   |              ✕              |
|         Gender          |       M/F       |                                        The gender of patient.                                        |              ✕              |
|         Therapy         |    character    |                                 The anti-tumor therapy on patients.                                  |              ✕              |
|        age_start        |     numeric     |                                     When the tumor is diagnosed.                                     |              ✕              |
|       tumor_type        |    character    |                                          The type of tumor.                                          |              ✕              |
|        seq_type         |    character    |                                         The sequencing type.                                         |              ✕              |
|           id            |    character    |                                    The identification of dataset.                                    |              ✕              |
    
</div>



## Obtain data from tigeR web server
```
Dataloader(c(1,2,3), use_source="Web Server")
```
## Obtain data from ExperimentHub
```
Dataloader(c(1,2,3), use_source="ExperimentHub)
```
