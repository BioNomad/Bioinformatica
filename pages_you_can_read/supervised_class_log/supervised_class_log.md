Supervised Learning - Classification/Logistic Regression
================

So in the last topic note, we explored predicting a continuous variable
with another continuous variable. But what if we wanted to predict a
categorical outcome from a continuous variable? As you may recall from
the alternative models topic notes, logistic regression is quite suited
to such a task. So let’s make a model\!

``` r
#load the data
load("./lgg.rda")
library(caret)
library(ggplot2)
library(hrbrthemes)
library(viridis)
#grab age and IDH1 mutant status
age <- lgg$PatientData$age_at_index
idh1_mutant <- as.numeric(gsub("Mutant",1,gsub("WT",0,lgg$PatientData$paper_IDH.status)))
status <- lgg$PatientData$paper_IDH.status
df <- data.frame(
  idh1_mutant,
  age,
  status
)
df <- na.omit(df)
#split data into training and test data
set.seed(42)
train_ind <- sample(1:nrow(df), 0.8*nrow(df))
train <- df[train_ind, ]
test  <- df[-train_ind, ]
#make the model
model <- glm(idh1_mutant~age,data = train,family = "binomial")
#our are terms significant?
summary(model)$coef
```

    ##               Estimate Std. Error   z value     Pr(>|z|)
    ## (Intercept)  4.3770834 0.53195066  8.228363 1.897885e-16
    ## age         -0.0631723 0.01055479 -5.985180 2.161512e-09

``` r
#maybe a quick visualization our our data?
ggplot(df, aes(x=status,y=age,fill=status)) + 
  geom_violin(alpha=.4)+
  scale_fill_manual(values = c(viridis(10)[3],viridis(10)[7]))
```

![](supervised_class_log_files/figure-gfm/data-1.svg)<!-- -->

So our model is trying to predict IDH1 mutant status using age. We
construct our logistic regression model and can see that the `P(>|z|)`
value is indicative that this term is indeed significant. A visual
examination also confirms there is a difference in the ages of patients
with IDH1 mutations and those without. So now let’s get some model
stats\!

``` r
#awesome! they are!
#so let's see how our model classifies
model_pred <- ifelse(predict(model, type = "link") > 0, "1", "0")
stats <- confusionMatrix(table(predicted = model_pred, actual = train$idh1_mutant), positive = "1")
#let's look at some stats
c(stats$overall["Accuracy"],
stats$byClass["Sensitivity"],
stats$byClass["Specificity"])
```

    ##    Accuracy Sensitivity Specificity 
    ##  0.81173594  0.98198198  0.06578947

``` r
#what about the classification error?
ClassError = function(actual, predicted) {
  mean(actual != predicted)
}
ClassError(train$idh1_mutant,model_pred)
```

    ## [1] 0.1882641

Here we see that the accuracy is relatively high at `0.81173594` and
same with the sensitivity at `0.98198198`. It should also be noted that
the classification error is `1-Accuracy`. The specificity is extremely
low at `0.06578947`. While the model accuracy and error are a little
more straightforward, what are sensistivity and specificity? Going back
to the topic note on prediction
([Wikipedia](https://en.wikipedia.org/wiki/Sensitivity_and_specificity)):

  - **Sensitivity/True Positive Rate** = True Positives/Actual Positives
    = True Positives/(True Positives + False Negatives)

  - **Specificity/ True Negative Rate** = True Negatives/Actual
    Negatives = True Negatives/(True Negatives + False Positives) = 1 -
    False Positive Rate

So here we can see that sensitivity is the ability to classify positives
and specificity is the ability to classify negatives. So this model is
great at classifying those with IDH1 mutant from age but isn’t great at
catching those without the IDH1 mutation from age. Which makes sense,
your age doesn’t necessarily mean you have an IDH1 mutation, this metric
isn’t **specific** to IDH1 mutations.

``` r
#let's plot our model information
ggplot(df, aes(x=age, y=idh1_mutant)) + geom_point(color="violetred4") + 
  stat_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)+
  geom_vline(xintercept = -coef(model)[1] / coef(model)[2])+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = 1)+
  geom_hline(yintercept = 0.5)
```

![](supervised_class_log_files/figure-gfm/plot-1.svg)<!-- -->

``` r
#what is our decision boundry?
-coef(model)[1] / coef(model)[2]
```

    ## (Intercept) 
    ##    69.28802

Here we plot our model information and calculate the **decision
boundary**. In classification is the boundary that separates the
classes; in this case it is `69.28802`
([Wikipedia](https://en.wikipedia.org/wiki/Decision_boundary)). Enough
of evaluating the training model, let’s test it\!

``` r
#let's test our model on test data
test_pred <- ifelse(predict(model,newdata = test, type = "response") > 0.5, "1", "0")
#let's get the stats
test_stats <- confusionMatrix(table(predicted = test_pred, actual = test$idh1_mutant), positive = "1")
c(test_stats$overall["Accuracy"],
test_stats$byClass["Sensitivity"],
test_stats$byClass["Specificity"])
```

    ##    Accuracy Sensitivity Specificity 
    ##  0.79611650  0.95294118  0.05555556

``` r
#now something new
#let's plot specificty against sensitivity
library(pROC)
test_prob <- predict(model, newdata = test, type = "response")
test_roc <- roc(test$idh1_mutant ~ test_prob, plot = TRUE, print.auc = TRUE)
```

![](supervised_class_log_files/figure-gfm/test-1.svg)<!-- -->

So we tested our model and got similar stats to our training model\! We
go further here and plot specificity against sensitivity. This is a
special plot that allows us to calculate the area under the curve (AUC)
or the **Accuracy** of our model. We look for a value above 0.5 as we
would like our model to be be accurate more than 50% of the time. If it
isn’t, the model is essentially worse than randomly choosing a class.

## References

1.  <https://en.wikipedia.org/wiki/Sensitivity_and_specificity>

2.  <https://en.wikipedia.org/wiki/Decision_boundary>
