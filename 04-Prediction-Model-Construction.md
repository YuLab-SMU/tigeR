# 🤖 Prediction Model Construction

## Build Model
 You can use function **build_Model** to build machine learning predictive model including naive Bayes models, Support vector machine model, random forest model, Cancerclass model, Adaboost model, Logitboost model, Logistics regression model.
 The function returns a trained model.

```         
train_set <- list(MEL_GSE91061, MEL_phs000452, RCC_Braun_2020)
mymodel <- build_Model(Model='NB', SE=train_set, feature_genes=Stem.Sig, response_NR = TRUE)
```

## Test Model
 **test_Model** are designed for testing model built by **build_Model** function. This function returns a ROC object.

```
test_Model(mymodel,MEL_78220)
```

