Error of the Mean/ Central Limit Theorem/ Confidence Intervals
================

## Error of the Mean

Here we will explore error of the mean, the central limit theorem and
confidence intervals. Before we dive into error of the mean it is
helpful to mention variance or the expectation that the squared
deviation of a random variable from the mean
([Wikipedia](https://en.wikipedia.org/wiki/Variance)).

<center>

**Variance = E\[(X - u)^2\]**

</center>

Where X is the random variable and u is the mean. Variance is also
related to standard deviation
([Wikipedia](https://en.wikipedia.org/wiki/Variance)):

<center>

**Standard Deviation = square root (Variance)**

</center>

So both variance and standard deviation are measures of spread of the
data away from the average value. When dealing with statistics like say
the mean, the mean calculated from a sample will have a standard error
of the mean value of
([Wikipedia](https://en.wikipedia.org/wiki/Standard_error)):

<center>

**Standard Error of the Mean = Standard Deviation / square root (number
of observations)**

</center>

Basically, this error term is the standard deviation of the sampling
distribution and is the estimate of the true standard deviation
([Wikipedia](https://en.wikipedia.org/wiki/Standard_error)). So we can
see, to decrease the error we need to increase the number of
observations. In this way, we get closer to approximating the true
underlying distribution. So to illustrate:

``` r
library(ggplot2)
library(dplyr)
#lets grab some sample sizes, set the vectors to fill, and the
#number of simulations
sample.sizes <- 2^(0:8)
means <- numeric()
standard.error.mean <- numeric()
theoretical <- 1/sqrt(sample.sizes)
number.of.simulations <- 500

#now we will loop through make a matrix where each sample has 500
#observations and test the error of the mean with different sample
#sizes.
for (i in sample.sizes){
  matrix.i <- matrix(rnorm(i*number.of.simulations),
                     nrow=number.of.simulations,ncol=i)
  means.tmp <- apply(matrix.i,1,mean)
  means <- c(means,mean(means.tmp))
  standard.error.mean <- c(standard.error.mean,sd(means.tmp))
}

means.df <- cbind.data.frame(
  sample.sizes = c(rep(sample.sizes,length(means)),
                   rep(sample.sizes,length(standard.error.mean)),
                   rep(sample.sizes,length(theoretical))),
  values = c(means,
             standard.error.mean,
             theoretical),
  group = c(rep("means",length(means)),
            rep("standard error of mean",length(standard.error.mean)),
            rep("theoretical",length(theoretical)))
)
p1 <- means.df %>%
  ggplot(aes(x=sample.sizes, y=values,
             group=group, color=group)) +
  geom_line()+
  ggtitle("")
p1
```

![](error_clt_ci_files/figure-gfm/sem-1.svg)<!-- -->

So we can see that increasing sample size will decrease the standard
error of the mean.

## Central limit theorem

Increasing sample size is also integral to the central limit theorem.
The central limit theorem states that the sum of random variables added
together will approach a normal distribution as sample size increases
([Wikipedia](https://en.wikipedia.org/wiki/Central_limit_theorem)). The
easiest way to visualize this is with histograms:

``` r
library(ggplot2)
library(gridExtra)

#generate histograms of different sample sizes
MakeHistogram <- function(sample.size,color){
  data <- data.frame(value=rnorm(sample.size))
  p <- ggplot(data, aes(x=value)) + 
    geom_histogram(fill=color)+
    ggtitle(paste("sample size = ",as.character(sample.size)))
  return(p)
}
grid.arrange(MakeHistogram(4,"thistle"),
             MakeHistogram(8,"thistle"),
             MakeHistogram(16,"thistle"),
             MakeHistogram(32,"seagreen1"),
             MakeHistogram(64,"seagreen1"),
             MakeHistogram(128,"seagreen1"),
             MakeHistogram(256,"paleturquoise4"),
             MakeHistogram(512,"paleturquoise4"),
             MakeHistogram(1024,"paleturquoise4"),
             ncol=3)
```

![](error_clt_ci_files/figure-gfm/clt-1.svg)<!-- -->

So we see that as sample size increases the distribution better
approximates a normal distribution or “bell-curve”.

## Confidence Intervals

A confidence interval is a range for an unknown parameter, with an
associated confidence level that the *true* value is within that range.

![ci](images/ci.PNG)

Where x is the sample mean, t is the critical value of the
t-distribution, o is the sample standard deviation, and n is the number
of observations
([Wikipedia](https://en.wikipedia.org/wiki/Confidence_interval)). In R
we can use the `confint()` function to accomplish this.

``` r
#let's create a simple model based on random data
data <- rnorm(200)
lm.data <- lm(data~1)

#now let's get that confidence interval for the intercept
confint(lm.data,level = 0.95)
```

    ##                  2.5 %   97.5 %
    ## (Intercept) -0.1583668 0.110836
