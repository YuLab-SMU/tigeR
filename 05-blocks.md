# 🔮 Response Prediction
**pred_response** function predict the immunotherapy response for the patients from gene expression data using our pre-trained machine learning models or public gene expression signatures.

```
pred_response(SE=MEL_GSE93157,Signature = ipt,
method = "Weighted_mean",threshold = 0.8,
PT_drop = FALSE,sort_by = "Customed.Signature",
group_by = "Customed.Signature",show.real = TRUE,
rankscore = FALSE)
```
