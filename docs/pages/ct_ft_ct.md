Contingency Table/Fisher Test/Chi-Square Test
================

So far we have been dealing with continuous variables - where the values
well are continuous and not discrete. So counts would be a discrete
value, for example having 10 patients. Gene expression is an example of
a continuous variable as you can have values in between, say like having
a value of 10, 10.3, 10.5, 11, 11.1, etc.. So what if we want to compare
two groups of categorical variables? The contingency table can help us
here\! A contingency table shows us the frequency distribution between
variables
([Wikipedia](https://en.wikipedia.org/wiki/Contingency_table)). Let’s
try it out with our glioma data\!

``` r
load("./lgg.rda")

#let's look at gender versus history of prior malignancy
table(
  lgg$PatientData$gender,
  lgg$PatientData$prior_malignancy
  )
```

    ##         
    ##           no yes
    ##   female 233   8
    ##   male   286   6

So here we see that out of all the female patients 233 did not have a
history of prior malignancy and 8 did. Conversely, 286 male patients did
not have a history of prior malignancy and 6 did. Now looking at this we
can see that there really isn’t a big difference between these
categorical variables. But what if we wanted a p-value out of this
table? Enter the Chi-square test\! But how does it work? Well for that
we need to look at the Chi-square test statistic:

![chi\_test](images/chi_test.PNG)

  - Where O is the observed value and E is the expected value
  - degrees of freedom = (number of rows - 1)(number of columns - 1)

Let’s apply this in R:

``` r
load("./lgg.rda")

#now let's use the chi-square test!
chisq.test(
table(
  lgg$PatientData$gender,
  lgg$PatientData$prior_malignancy
  )
)
```

    ## 
    ##  Pearson's Chi-squared test with Yates' continuity correction
    ## 
    ## data:  table(lgg$PatientData$gender, lgg$PatientData$prior_malignancy)
    ## X-squared = 0.40523, df = 1, p-value = 0.5244

``` r
#what about getting the expected values?
chisq.test(
table(
  lgg$PatientData$gender,
  lgg$PatientData$prior_malignancy
  )
)$expected
```

    ##         
    ##                no      yes
    ##   female 234.6698 6.330206
    ##   male   284.3302 7.669794

So, as we expected the p-value tells us there isn’t a significant
difference between these categorical variables. But, the chi-square test
is really suited to small sample sizes. For that we need the Fisher
Exact test. What is special about the Fisher Exact test is that it is
non-parametric, meaning it makes no assumption about the underlying
distribution
([Wikipedia](https://en.wikipedia.org/wiki/Fisher%27s_exact_test)).
Which is good for small sample sizes given that the underlying
distribution is better approximated by larger sample sizes. Let’s take a
peak under the hood.

Using our example with gender v. primary diagnosis:

![ct](images/ct.PNG)

![fisher\_p](images/fisher_p.PNG)

While it would be fun to calculate these factorials, we can leverage the
`fisher.test()` function in R:

``` r
load("./lgg.rda")

#now let's use the fisher test!
fisher.test(
table(
  lgg$PatientData$gender,
  lgg$PatientData$prior_malignancy
  )
)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  table(lgg$PatientData$gender, lgg$PatientData$prior_malignancy)
    ## p-value = 0.4212
    ## alternative hypothesis: true odds ratio is not equal to 1
    ## 95 percent confidence interval:
    ##  0.1723462 2.0422983
    ## sample estimates:
    ## odds ratio 
    ##  0.6115911

## References

1.  <https://en.wikipedia.org/wiki/Contingency_table>

2.  <https://en.wikipedia.org/wiki/Fisher%27s_exact_test>
