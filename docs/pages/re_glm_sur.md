Random Effects/General Linear Models/Survival Analysis
================

Thus far we have only talked about linear models. And as we saw in the
topic note on the 2-way ANOVA, that might not be the best fit for the
data. So let’s talk about some other models. We will start with the
random effects models, but before we jump into random effects let’s talk
about fixed effects. Fixed effects are well, fixed. They are part of the
study design and non-random
([Wikipedia](https://en.wikipedia.org/wiki/Fixed_effects_model)). Groups
like gender and diagnosis would be fixed. Random effects are those that
you don’t control for, they are random
([Wikipedia](https://en.wikipedia.org/wiki/Random_effects_model)). So
how do you go about modeling something random? We have an equation for
that ([Wikipedia](https://en.wikipedia.org/wiki/Random_effects_model))\!

![re](images/re.PNG)

Where **u** is the grand mean, **i** is the group, **j** is the
individual observation, **a<sub>i</sub>** is the random effect for the
group, and **e<sub>ij</sub>** is the individual observation random
effect
([Wikipedia](https://en.wikipedia.org/wiki/Random_effects_model)). In a
random effects model, there is the assumption of random effects between
groups and that observations within a group are correlated
([Wikipedia](https://en.wikipedia.org/wiki/Intraclass_correlation)). We
can measure the degree of this correlation with the intraclass
correlation coefficient
([Wikipedia](https://en.wikipedia.org/wiki/Intraclass_correlation)):

![icc](images/icc.PNG)

Where **o<sub>a</sub>** is the variance of the group, and
**o<sub>e</sub>** is the variance of the individual observation. We can
estimate these using the `lme4` library. So let’s try it out with our
LGG data\!

``` r
load("./lgg.rda")
library(lme4)
library(limma)
library(ggplot2)
#remember to normalize!
log.trans <- log2(lgg$ExpressionData + 1)
norm.data <- normalizeQuantiles(log.trans)

#now what if we ask the question what effect does the hospital have on
#IDH1 expression, are there random effects between hospitals?
hospitals <- as.factor(
  as.character(lgg$PatientData$paper_Tissue.source.site)
  )
idh1 <- idh1 <- as.numeric(norm.data["IDH1|3417",])
df <- data.frame(
  hospitals = hospitals,
  idh1=idh1
)
#let's take a quick peek at the data
ggplot(df, aes(x=hospitals, y=idh1,fill=hospitals)) + 
  geom_violin()+
  coord_flip()+ 
  theme(legend.position="bottom",
        legend.key.size = unit(.1, 'cm'))
```

![](re_glm_sur_files/figure-gfm/re-1.svg)<!-- -->

``` r
random.fit <- lmer(idh1 ~ (1 | hospitals), data = df)
summary(random.fit)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: idh1 ~ (1 | hospitals)
    ##    Data: df
    ## 
    ## REML criterion at convergence: 955.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.6444 -0.5114  0.0380  0.6570  2.7728 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  hospitals (Intercept) 0.01021  0.1011  
    ##  Residual              0.36350  0.6029  
    ## Number of obs: 516, groups:  hospitals, 26
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  11.9792     0.0384     312

Visualizing by hospital does reveal quite a bit of variation between
them. But how do we interpret the results of the `lmer()` function? For
that we turn to the `Random effects` section. The variance of the
hospitals is our **o<sub>a</sub>** term and the variance of the
residuals is our **o<sub>e</sub>**. So plugging that into the intraclass
correlation coefficient we get `0.01021 /(0.36350 + 0.01021)
= 0.02732065`. So while there is some visual variation between
hospitals, only \~3% of the variance can be explained by different
hospitals.

## Generalized Linear Model

Going back to linear models, we have noted several times that the models
constructed fail to have their residuals normally distributed. This is
indicative that the relationship isn’t linear
([Wikipedia](https://en.wikipedia.org/wiki/Generalized_linear_model)).
So what do we do in these cases? Well considering the title of section
is the generalized linear model, you may have guessed this is the
answer. Generalized linear models are special in that they use a *link*
function. A link function relates the response variable to the model and
that function varies linearly with the predictors
([Wikipedia](https://en.wikipedia.org/wiki/Generalized_linear_model)).
By doing this you get away with assuming the response variable and
predictor vary linearly. While we can stick with a typical gaussian
function, we have the opportunity to explore logistic regression.

So what is so special about logistic regression? Well logistic
regression differs from linear regression in that logistic regression
models a binary outcome
([STHDA](http://www.sthda.com/english/articles/36-classification-methods-essentials/151-logistic-regression-essentials-in-r/)).
So predicting a yes or no/ 1 or 0. Logistic regression models are linked
by a logit function
([STHDA](http://www.sthda.com/english/articles/36-classification-methods-essentials/151-logistic-regression-essentials-in-r/)):

![logreg](images/logreg.PNG)

Where **p** is the probability of y=1 and the **B** terms are the
parameters of the model. Let’s try it in R\!

``` r
#normalized IDH1 gene expression
idh1 <- idh1 <- as.numeric(norm.data["IDH1|3417",])
#now grab prior malignancy
prior.mal <- as.factor(
  lgg$PatientData$prior_malignancy
)
df1 <- data.frame(
  idh1,
  prior.mal
)
glm.prior <- glm(prior.mal~idh1, data = df1, family = "binomial")
summary(glm.prior)$coef
```

    ##               Estimate Std. Error    z value  Pr(>|z|)
    ## (Intercept) -7.7536696   5.571221 -1.3917360 0.1640023
    ## idh1         0.3439527   0.460009  0.7477087 0.4546359

To better understand our model we turn to the coefficients
([STHDA](http://www.sthda.com/english/articles/36-classification-methods-essentials/151-logistic-regression-essentials-in-r/)):

  - Our `Estimate` gives us the Beta value for the coefficients in our
    model
  - The `Std. Error` is the standard error of these model coefficients
    (larger `Std. Error` = less confidence in the estimate )
  - The `z value` is the z-statistic (`z value` = `Estimate` / `Std.
    Error`)
  - And finally, `Pr(>|z|)` is the p-value associated with the
    z-statistic (smaller `Pr(>|z|)` = more significant the estimate)

Given these we can see right away that the p-value indicates that our
model estimation isn’t terribly significant. Beta coefficients in a
logistic regression model can also offer us an odds ratio - a measure of
the association between a predictor and an outcome
([STHDA](http://www.sthda.com/english/articles/36-classification-methods-essentials/151-logistic-regression-essentials-in-r/)).
So in our model we have a coefficient of 0.344 for IDH1 expression. This
is indicative that for every one unit of IDH1 gene expression the odds
that the subject has a prior history of malignancy increases by 1.344.

## Survival Analysis

When we talk about Survival analysis, we often hear about Kaplan-Meier
survival curves. These are a good way of visualizing the survival
probability over time. The following is the probability of survival at
time i, **S(t<sub>i</sub>)**, using the Kaplan-Meier method
([STHDA](http://www.sthda.com/english/wiki/survival-analysis-basics)):

![surv](images/surv.PNG)

Where, **S(t<sub>i - 1</sub>)** is the probability of being alive at
time **t<sub>i - 1</sub>**, **n<sub>i</sub>** is the number of patients
alive just before **t<sub>i</sub>**, **d<sub>i</sub>** is the number of
events at **t<sub>i</sub>** and **S<sub>0</sub>** = 1/ **t<sub>0</sub>**
= 0. Now not all patients in a study will be deceased at the end of a
study for one reason or another. These patients are noted as *censored*
. Now let’s get to coding\!

``` r
library(survival)
library(survminer)
#grab the censored status, gender, and time to death
censored <- as.numeric(
  gsub("(Not Reported)|(Alive)","1",
  gsub("Dead","2",
    lgg$PatientData$vital_status
    )
  )
)
gender <- lgg$PatientData$gender
time <- lgg$PatientData$paper_Survival..months.
#create the data frame
df2 <- data.frame(
  censored,
  gender,
  time
)
#construct the survival model and print some summary statistics
surv.model <- survfit(Surv(time, censored) ~ gender, data = df2)
print(surv.model)
```

    ## Call: survfit(formula = Surv(time, censored) ~ gender, data = df2)
    ## 
    ##    77 observations deleted due to missingness 
    ##                 n events median 0.95LCL 0.95UCL
    ## gender=female 201     46   75.1    62.0     114
    ## gender=male   256     61   75.0    50.8     124

``` r
ggsurvplot(surv.model,
          pval = TRUE, conf.int = TRUE,
          risk.table = TRUE, 
          risk.table.col = "strata", 
          linetype = "strata", 
          surv.median.line = "hv", 
          ggtheme = theme_bw() 
          )
```

![](re_glm_sur_files/figure-gfm/surv-1.svg)<!-- -->

Here we can see there really isn’t too much difference between the
survival curves between males and females. The summary of our model
demonstrates that the median survival time for both genders is about 75
months and shows that their confidence intervals overlap considerably.
We also see the number of males/females and the number of events per
gender. Using the `ggsurvplot()` function we can visualize these curves
and by setting `risk.table` to true we can see the number of individuals
at risk at a certain time. The p-value generated is from the **Log-Rank
Test**, which tests for differences between survival curves between two
or more groups
([STHDA](http://www.sthda.com/english/wiki/survival-analysis-basics)).
We could also accomplish this using the `survdiff()` function.

However, these tests (Kaplan-Meier and Log-Rank) are examples of
*univariate* analysis. Meaning, that they work when only *one* factor is
being investigated
([STHDA](http://www.sthda.com/english/wiki/cox-proportional-hazards-model)).
Additionally, that one factor needs to be categorical. In our last
example that variable was gender. Now what if you wanted to test
something quantitative, like IDH1 gene expression? We would leverage Cox
proportional hazards regression analysis - which can take more than one
predictor variable and both categorical and quantitative variables
([STHDA](http://www.sthda.com/english/wiki/cox-proportional-hazards-model)).
The Cox model is the following hazard function - put another way the
risk of dying at time t
([STHDA](http://www.sthda.com/english/wiki/cox-proportional-hazards-model)):

![cox](images/cox.PNG)

Where, **t** is time, **h(t)** is the hazard function for p covariates,
the **b** coefficients are the measure of impact each one has, and
**h<sub>0</sub>** is the baseline hazard and is the hazard if all the
covariates were equal to 0. Enough math though, let’s turn to R\!

``` r
library(survival)
library(survminer)
#univariate analysis
cox.model.uni <- coxph(Surv(time, censored) ~ gender, data = df2)
cox.model.uni
```

    ## Call:
    ## coxph(formula = Surv(time, censored) ~ gender, data = df2)
    ## 
    ##               coef exp(coef) se(coef)     z     p
    ## gendermale 0.08782   1.09179  0.19649 0.447 0.655
    ## 
    ## Likelihood ratio test=0.2  on 1 df, p=0.6543
    ## n= 457, number of events= 107 
    ##    (77 observations deleted due to missingness)

``` r
ggforest(cox.model.uni,data = df2)
```

![](re_glm_sur_files/figure-gfm/cox-1.svg)<!-- -->

``` r
#let's try multivariate analysis
age <- lgg$PatientData$age_at_diagnosis
df3 <- data.frame(
  idh1,
  gender,
  time,
  censored,
  age
)
cox.model.multi <- res.cox <- coxph(Surv(time, censored) ~ age + gender + idh1, data =  df3)
cox.model.multi
```

    ## Call:
    ## coxph(formula = Surv(time, censored) ~ age + gender + idh1, data = df3)
    ## 
    ##                 coef exp(coef)  se(coef)     z        p
    ## age        1.554e-04 1.000e+00 2.185e-05 7.114 1.13e-12
    ## gendermale 7.747e-02 1.081e+00 2.005e-01 0.386   0.6992
    ## idh1       4.828e-01 1.621e+00 1.818e-01 2.656   0.0079
    ## 
    ## Likelihood ratio test=58.25  on 3 df, p=1.389e-12
    ## n= 456, number of events= 107 
    ##    (78 observations deleted due to missingness)

``` r
ggforest(cox.model.multi,data = df3)
```

![](re_glm_sur_files/figure-gfm/cox-2.svg)<!-- -->

Well, there is quite a bit to unpack here. The `coef` column contains
the regression coefficients. So in our univariate analysis we see that
`gendermale` has a coef of 0.08782, indicating that males have a
slightly higher risk of death. However, less stock should be put into
gender as the p-value doesn’t suggest this term is significant. The
`exp(coef)` column is the hazard ratio or effect of the covariate. `z`
is the z-score and `p` is the corresponding p-value. These models can be
better visualized by the `ggforest()` function, where we see that IDH1
expression increase the risk of death.

## References

1.  <https://en.wikipedia.org/wiki/Fixed_effects_model>

2.  <https://en.wikipedia.org/wiki/Random_effects_model>

3.  <https://en.wikipedia.org/wiki/Intraclass_correlation>

4.  <https://en.wikipedia.org/wiki/Generalized_linear_model>

5.  <http://www.sthda.com/english/articles/36-classification-methods-essentials/151-logistic-regression-essentials-in-r/>

6.  <http://www.sthda.com/english/wiki/survival-analysis-basics>

7.  <http://www.sthda.com/english/wiki/cox-proportional-hazards-model>
