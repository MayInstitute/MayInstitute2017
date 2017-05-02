---
title: 'Day 3: Beginner''s statistics in R'
author: "Laurent Gatto and Meena Choi"
output: 
  html_document:
    self_contained: true
    toc: true
    toc_float: true
    fig_caption: no
---




# Objectives

- Hypothesis testing: t-test
- Categorical data
- Power calculation
- Linear models and correlation

---

# Part 4: Statistical hypothesis test

First, we are going to prepare the session for further analyses.


```r
load("./data/summaryresults.rda")
load("./data/iprg.rda")
```

## Two sample t-test for one protein with one feature

Now, we'll perform a t-test whether protein `sp|P44015|VAC2_YEAST` has
a change in abundance between Condition 1 and Condition 2.

### Hypothesis

* $H_0$: no change in abundance, mean(Condition1) - mean(Condition2) = 0
* $H_a$: change in abundance, mean(Condition1) - mean(Condition 2) $\neq$ 0

### Statistics

* Observed $t = \frac{\mbox{difference of group means}}{\mbox{estimate of variation}} = \frac{(mean_{1} - mean_{2})}{SE} \sim t_{\alpha/2, df}$
* Standard error, $SE=\sqrt{\frac{s_{1}^2}{n_{1}} + \frac{s_{2}^2}{n_{2}}}$

with 

* $n_{1}$: number of replicates
* $s_{1}^2 = \frac{1}{n_{1}-1} \sum (Y_{1i} - \bar{Y_{1 \cdot}})^2$: sample variance

### Data preparation


```r
## Let's start with one protein, named "sp|P44015|VAC2_YEAST"
oneproteindata <- iprg[iprg$Protein == "sp|P44015|VAC2_YEAST", ]

## Then, get two conditions only, because t.test only works for two
## groups (conditions).
oneproteindata.condition12 <- oneproteindata[oneproteindata$Condition %in% 
                                             c('Condition1', 'Condition2'), ]
oneproteindata.condition12
```

```
##                    Protein Log2Intensity                       Run
## 21096 sp|P44015|VAC2_YEAST      26.30163 JD_06232014_sample1_B.raw
## 21097 sp|P44015|VAC2_YEAST      26.11643 JD_06232014_sample1_C.raw
## 21098 sp|P44015|VAC2_YEAST      26.29089 JD_06232014_sample1-A.raw
## 21099 sp|P44015|VAC2_YEAST      25.81957 JD_06232014_sample2_A.raw
## 21100 sp|P44015|VAC2_YEAST      26.11527 JD_06232014_sample2_B.raw
## 21101 sp|P44015|VAC2_YEAST      26.08498 JD_06232014_sample2_C.raw
##        Condition BioReplicate Intensity
## 21096 Condition1            1  82714388
## 21097 Condition1            1  72749239
## 21098 Condition1            1  82100518
## 21099 Condition2            2  59219741
## 21100 Condition2            2  72690802
## 21101 Condition2            2  71180513
```

```r
table(oneproteindata.condition12[, c("Condition", "BioReplicate")])
```

```
##             BioReplicate
## Condition    1 2
##   Condition1 3 0
##   Condition2 0 3
##   Condition3 0 0
##   Condition4 0 0
```

If we want to remove the levels that are not relevant anymore, we can
use `droplevels`:


```r
oneproteindata.condition12 <- droplevels(oneproteindata.condition12)
table(oneproteindata.condition12[, c("Condition", "BioReplicate")])
```

```
##             BioReplicate
## Condition    1 2
##   Condition1 3 0
##   Condition2 0 3
```

To perform the t-test, we use the `t.test` function. Let's first
familiarise ourselves with it by looking that the manual 


```r
?t.test
```

And now apply to to our data


```r
# t test for different abundance (log2Int) between Groups (Condition)
result <- t.test(Log2Intensity ~ Condition,
                 data = oneproteindata.condition12,
                 var.equal = FALSE)

result
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  Log2Intensity by Condition
## t = 2.0608, df = 3.4001, p-value = 0.1206
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.1025408  0.5619598
## sample estimates:
## mean in group Condition1 mean in group Condition2 
##                 26.23632                 26.00661
```

> **Challenge**
>
> Repeat the t-test above but with calculating a 90% confidence interval
> for the log2 fold change.



### The `htest` class

The `t.test` function, like other hypothesis testing function, return
a result of a type we haven't encountered yet, the `htest` class:


```r
class(result)
```

```
## [1] "htest"
```

which stores typical results from such tests. Let's have a more
detailed look at what information we can learn from the results our
t-test. When we type the name of our `result` object, we get a short
textual summary, but the object contains more details:


```r
names(result)
```

```
## [1] "statistic"   "parameter"   "p.value"     "conf.int"    "estimate"   
## [6] "null.value"  "alternative" "method"      "data.name"
```

and we can access each of these by using the `$` operator, like we
used to access a single column from a `data.frame`, but the `htest`
class is not a `data.frame` (it's actually a `list`). For example, to
access the group means, we would use


```r
result$estimate
```

```
## mean in group Condition1 mean in group Condition2 
##                 26.23632                 26.00661
```

> **Challenge**
> 
> * Calculate the (log2-transformed) fold change between groups
> * Extract the value of the t-statistics
> * Calculate the standard error (fold-change/t-statistics)
> * Extract the degrees of freedom (parameter)
> * Extract the p values
> * Extract the 95% confidence intervals




We can also manually compute our t-test statistic using the formulas
we descibed above and compare it with the `summaryresult`.

Recall the `summaryresult` we generated last section


```r
summaryresult
```

```
##        Group     mean         sd         se length ciw.lower.95
## 1 Condition1 26.23632 0.10396539 0.06002444      3     25.97805
## 2 Condition2 26.00661 0.16268179 0.09392438      3     25.60248
## 3 Condition3 23.25609 0.09467798 0.05466236      3     23.02090
## 4 Condition4 20.97056 0.73140174 0.42227499      3     19.15366
##   ciw.upper.95
## 1     26.49458
## 2     26.41073
## 3     23.49128
## 4     22.78746
```

```r
summaryresult12 <- summaryresult[1:2, ]

## test statistic, It is the same as 'result$statistic' above.
diff(summaryresult12$mean) ## same as result$estimate[1]-result$estimate[2]
```

```
## [1] -0.2297095
```

```r
sqrt(sum(summaryresult12$sd^2/summaryresult12$length)) ## same as stand error
```

```
## [1] 0.1114662
```

```r
## the t-statistic
diff(summaryresult12$mean)/sqrt(sum(summaryresult12$sd^2/summaryresult12$length))
```

```
## [1] -2.060799
```

## Re-calculating the p values

See the previous
[*Working with statistical distributions*](https://htmlpreview.github.io/?https://github.com/MayInstitute/MayInstitute2017/blob/master/Program3_Intro%20stat%20in%20R/02-rstats.html#working_with_statistical_distributions)
for a reminder.

Referring back to our t-test results above, we can manually calculate
the one- and two-side tests p-values using the t-statistics and the
test parameter (using the `pt` function).


Our result t statistic was 2.0607988 (accessible
with `result$statistic`). Let's start by visualising it along a t
distribution. Let's create data from such a distribution, making sure
we set to appropriate parameter.


```r
xt <- rt(1e5, result$parameter)
plot(density(xt), xlim = c(-10, 10))
abline(v = result$statistic, col = "red")
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png)

The area on the left of that point is given by `pt(result$statistic,
result$parameter)`, which is 0.939693. The p-value for a one-sided test is this given by


```r
1 - pt(result$statistic, result$parameter)
```

```
##          t 
## 0.06030697
```

And the p-value for a two-sided test is 


```r
2 * (1 - pt(result$statistic, result$parameter))
```

```
##         t 
## 0.1206139
```

which is the same as the one calculated by the t-test.



# Part 5a: Sample size calculation

To calculate the required sample size, you’ll need to know four
things:

* $\alpha$: confidence level
* $power$: 1 - $\beta$, where $\beta$ is probability of a true positive discovery
* $\Delta$: anticipated fold change
* $\sigma$: anticipated variance

## R code

Assuming equal varaince and number of samples across groups, the
following formula is used for sample size estimation:

$$\frac{2{\sigma}^2}{n}\leq(\frac{\Delta}{z_{1-\beta}+z_{1-\alpha/2}})^2$$



```r
library("pwr")

## ?pwr.t.test

# Significance level alpha
alpha <- 0.05

# Power = 1 - beta
power <- 0.95

# anticipated log2 fold change 
delta <- 1

# anticipated variability
sigma <- 1.5

# Effect size
# It quantifies the size of the difference between two groups
d <- delta/sigma

#Sample size estimation
pwr.t.test(d = d, sig.level = alpha, power = power, type = 'two.sample')
```

```
## 
##      Two-sample t test power calculation 
## 
##               n = 59.45415
##               d = 0.6666667
##       sig.level = 0.05
##           power = 0.95
##     alternative = two.sided
## 
## NOTE: n is number in *each* group
```

Then, we investigate the effect of required fold change and variance on the sample size estimation.


```r
# anticipated log2 fold change 
delta <- seq(0.1, 0.7, .1)
nd <- length(delta)

# anticipated variability
sigma <- seq(0.1,0.5,.1)
ns <- length(sigma)

# obtain sample sizes
samsize <- matrix(0, nrow=ns*nd, ncol = 3)
counter <- 0
for (i in 1:nd){
  for (j in 1:ns){
    result <- pwr.t.test(d = delta[i]/sigma[j],
                         sig.level = alpha, power = power,
                         type = "two.sample")
    counter <- counter + 1
    samsize[counter,1] <- delta[i]
    samsize[counter,2] <- sigma[j]
    samsize[counter,3] <- ceiling(result$n)
  }
}
colnames(samsize) <- c("fd","var","value")


library("ggplot2")
samsize <- as.data.frame(samsize)
samsize$var <- as.factor(samsize$var)
ggplot(data=samsize, aes(x=fd, y=value, group = var, colour = var)) +
  geom_line() +
  geom_point(size=2, shape=21, fill="white") +
  labs(title="Sig=0.05 Power=0.05", x="Anticipated log2 fold change", y='Sample Size (n)') +
  theme(plot.title = element_text(size=20, colour="darkblue"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=13)) 
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-1.png)

# Part 5b: Comparison of two proportions

For this part, we are going to use a new dataset, which contains the
patient information from TCGA colorectal cohort. This data is from
*Proteogenomic characterization of human colon and rectal cancer*
(Zhang et al. 2014).

Rows in the data array are patients and columns are patient
information. The column definition is shown as following:

| Variable            |
|---------------------|
| TCGA participant ID |
| Gender              |
| Cancer type         |
| BTAF mutation status|
| History of colon polyps |

## Generate 2-way contingency tables

We first need to calculate 2-way contingency tables for the
statistical tests.


Let's read in and explore the TCGA colorectal cancer data:


```r
TCGA.CRC <- read.csv("./data/TCGA_sample_information.csv")
head(TCGA.CRC)
```

```
##   TCGA.participant.ID Gender Cancer BRAF.mutation history_of_colon_polyps
## 1        TCGA-A6-3807 Female  Colon             0                      NO
## 2        TCGA-A6-3808   Male  Colon             0                     YES
## 3        TCGA-A6-3810   Male  Colon             0                     YES
## 4        TCGA-AA-3518 Female  Colon             0                      NO
## 5        TCGA-AA-3525   Male  Colon             1                      NO
## 6        TCGA-AA-3526   Male  Colon             0                     YES
```

We are interested in the cancer type and history of colon polyps;
let's select columns from TCGA dataset:


```r
TCGA.CRC.gc <- TCGA.CRC[, c('Cancer', 'history_of_colon_polyps')]
nrow(TCGA.CRC.gc)
```

```
## [1] 78
```

and generate a 2-way contingency table using the `table` function and
visualise the count data


```r
ov <- table(TCGA.CRC.gc)
ov
```

```
##         history_of_colon_polyps
## Cancer   NO YES
##   Colon  31  22
##   Rectum 20   5
```

```r
dotchart(t(ov), xlab = "Observed counts")
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19-1.png)


## Chi-squared test

**Hypothesis** : 

$H_0$ : each population has the same proportion of observations, $\pi_{j=1|i=1} = \pi_{j=1|i=2}$

$H_a$ : different population have different proportion of observations


$$\chi^2 =\sum_{i=1}^2 \sum_{j=1}^2 \frac{(O_{ij}-E_{ij})^2}{E_{ij}} \sim \chi^2_{(2-1)(2-1)}$$

$O_{ij}$ : $n_{ij}$, which is the count within the cells

$E_{ij}$ : $n_{i+}n_{+j}/n$, where $n_{i+}$ is the row count sum, $n_{+j}$ is the column count sum and n is the total count.


**Hypothesis**: whether the proportion of patients who have history of
colon polyps in the patients with colon cancer is different from that
in the patients with rectal cancer



```r
pt <- prop.test(ov)
pt
```

```
## 
## 	2-sample test for equality of proportions with continuity
## 	correction
## 
## data:  ov
## X-squared = 2.5871, df = 1, p-value = 0.1077
## alternative hypothesis: two.sided
## 95 percent confidence interval:
##  -0.44991310  0.01972442
## sample estimates:
##    prop 1    prop 2 
## 0.5849057 0.8000000
```

```r
# name of output
names(pt)
```

```
## [1] "statistic"   "parameter"   "p.value"     "estimate"    "null.value" 
## [6] "conf.int"    "alternative" "method"      "data.name"
```

```r
# proportion in each group
pt$estimate 
```

```
##    prop 1    prop 2 
## 0.5849057 0.8000000
```

```r
# test statistic value
pt$statistic 
```

```
## X-squared 
##  2.587111
```

```r
# degree of freedom
pt$parameter
```

```
## df 
##  1
```

## Fisher's exact test

The Fisher's exact test can be used with small sample sizes. It
compares distributions of counts within the 4 cells.

> **Challenge**
>
> Apply the Fisher's exact test using the `fisher.test` function and
> extract the odds ratio.



## Z-test

Finally, we could also apply a z-test, defined as

$$Z=\frac{\widehat{\pi}_1-\widehat{\pi}_2}{\sqrt{\frac{\widehat{\pi}_1(1 - \widehat{\pi}_1)}{n_1}+\frac{\widehat{\pi_2}(1 - \widehat{\pi_2})}{n_2}}}$$

where $\widehat{\pi}_1 = \frac{y_{1}}{n_1}$ and $\widehat{\pi}_2 = \frac{y_{2}}{n_2}$.

We are going to use this test to illustrate how to write functions in
R. 

An R function is created with the function constructor, named
`function`, and is composed of:

1. a name: we will call our function to calculate the p-value
   `z.prop.p`;
2. some inputs: our inputs are the number of observations for the
   outcome of interest, `x1` and `x2`, and the total number of
   observation in each category, `n1` and `n2`;
3. a returned value (output): the computed p-value, named `pvalue`;
4. a body, i.e. the code that, given the inputs, calculates the output.



```r
z.prop.p <- function(x1, x2, n1, n2) {
    pi_1 <- x1/n1
    pi_2 <- x2/n2
    numerator <- pi_1 - pi_2
    denominator <- sqrt(((pi_1*(1-pi_1))/n1 + (pi_2*(1-pi_2))/n2))
    stat <- numerator/denominator
    pvalue <- 2 * (1 - pnorm(abs(stat)))
    return(pvalue)
}

z.prop.p(ov[1,1], ov[2,1], sum(ov[1,]), sum(ov[2,]))
```

```
## [1] 0.04010935
```

Same as above, for the confidence interval at a given alpha level.


```r
z.prop.ci <- function(x1, x2, n1, n2, alpha = 0.05){
  pi_1 <- x1/n1
  pi_2 <- x2/n2
  numerator <- pi_1 - pi_2
  denominator <- sqrt(((pi_1*(1-pi_1))/n1 + (pi_2*(1-pi_2))/n2))
  cri_value <- qnorm(1-alpha/2)
  prop.ci <- c(numerator + cri_value * denominator,
               numerator - cri_value * denominator)
  return(prop.ci)
}

z.prop.ci(ov[1,1], ov[2,1], sum(ov[1,]), sum(ov[2,]))
```

```
## [1] -0.009709539 -0.420479141
```

> **Challenge**
>
> Write a function named `f2c` (`c2f`) that converts a temperature
> from Fahrenheit to Celsium (Celsium to Fahrenheit) using the
> following formula $F = C \times 1 + 32$ ($C = \frac{F - 32}{1.8}$).


# Part 6: Linear models and correlation


When considering correlations and modelling data, visualisation is key. 

Let's use the famous
[*Anscombe's quartet*](https://en.wikipedia.org/wiki/Anscombe%27s_quartet)
data as a motivating example. This data is composed of 4 pairs of
values, $(x_1, y_1)$ to $(x_4, y_4)$:


| x1| x2| x3| x4|    y1|   y2|    y3|    y4|
|--:|--:|--:|--:|-----:|----:|-----:|-----:|
| 10| 10| 10|  8|  8.04| 9.14|  7.46|  6.58|
|  8|  8|  8|  8|  6.95| 8.14|  6.77|  5.76|
| 13| 13| 13|  8|  7.58| 8.74| 12.74|  7.71|
|  9|  9|  9|  8|  8.81| 8.77|  7.11|  8.84|
| 11| 11| 11|  8|  8.33| 9.26|  7.81|  8.47|
| 14| 14| 14|  8|  9.96| 8.10|  8.84|  7.04|
|  6|  6|  6|  8|  7.24| 6.13|  6.08|  5.25|
|  4|  4|  4| 19|  4.26| 3.10|  5.39| 12.50|
| 12| 12| 12|  8| 10.84| 9.13|  8.15|  5.56|
|  7|  7|  7|  8|  4.82| 7.26|  6.42|  7.91|
|  5|  5|  5|  8|  5.68| 4.74|  5.73|  6.89|

Each of these $x$ and $y$ sets have the same variance, mean and
correlation:




|         |          1|          2|          3|          4|
|:--------|----------:|----------:|----------:|----------:|
|var(x)   | 11.0000000| 11.0000000| 11.0000000| 11.0000000|
|mean(x)  |  9.0000000|  9.0000000|  9.0000000|  9.0000000|
|var(y)   |  4.1272691|  4.1276291|  4.1226200|  4.1232491|
|mean(y)  |  7.5009091|  7.5009091|  7.5000000|  7.5009091|
|cor(x,y) |  0.8164205|  0.8162365|  0.8162867|  0.8165214|

But...

While the *residuals* of the linear regression clearly indicate
fundamental differences in these data, the most simple and
straightforward approach is *visualisation* to highlight the
fundamental differences in the datasets.

![plot of chunk anscombefig](figure/anscombefig-1.png)

See also another, more recent example:
[The Datasaurus Dozen dataset](https://www.autodeskresearch.com/publications/samestats).
<details>

![The Datasaurus Dozen dataset](https://d2f99xq7vri1nk.cloudfront.net/DinoSequentialSmaller.gif)
</details>


## Correlation

Here is an example where a wide format comes very handy. We are going
to convert our iPRG data using the `spread` function from the `tidyr`
package:



```r
library("tidyr")
iprg2 <- spread(iprg[, 1:3], Run, Log2Intensity)
rownames(iprg2) <- iprg2$Protein
iprg2 <- iprg2[, -1]
```



And lets focus on the 3 runs, i.e. 2 replicates from condition
1 and 


```r
x <- iprg2[, c(1, 2, 10)]
head(x)
```

```
##                       JD_06232014_sample1-A.raw JD_06232014_sample1_B.raw
## sp|D6VTK4|STE2_YEAST                   26.58301                  26.81232
## sp|O13297|CET1_YEAST                   24.71809                  24.71912
## sp|O13329|FOB1_YEAST                   23.47075                  23.37678
## sp|O13539|THP2_YEAST                   24.29661                  27.52021
## sp|O13547|CCW14_YEAST                  27.11638                  27.22234
## sp|O13563|RPN13_YEAST                  26.17056                  26.09476
##                       JD_06232014_sample4-A.raw
## sp|D6VTK4|STE2_YEAST                   26.65573
## sp|O13297|CET1_YEAST                   24.50814
## sp|O13329|FOB1_YEAST                   23.03473
## sp|O13539|THP2_YEAST                   25.07576
## sp|O13547|CCW14_YEAST                  27.07526
## sp|O13563|RPN13_YEAST                  25.77958
```

```r
pairs(x)
```

![plot of chunk unnamed-chunk-26](figure/unnamed-chunk-26-1.png)

We can use the `cor` function to calculate the Pearson correlation
between two vectors of the same length (making sure the order is
correct), or a dataframe. But, we have missing values in the data,
which will stop us from calculating the correlation:


```r
cor(x)
```

```
##                           JD_06232014_sample1-A.raw
## JD_06232014_sample1-A.raw                         1
## JD_06232014_sample1_B.raw                        NA
## JD_06232014_sample4-A.raw                        NA
##                           JD_06232014_sample1_B.raw
## JD_06232014_sample1-A.raw                        NA
## JD_06232014_sample1_B.raw                         1
## JD_06232014_sample4-A.raw                        NA
##                           JD_06232014_sample4-A.raw
## JD_06232014_sample1-A.raw                        NA
## JD_06232014_sample1_B.raw                        NA
## JD_06232014_sample4-A.raw                         1
```

We first need to omit the proteins/rows that contain missing values


```r
x2 <- na.omit(x)
cor(x2)
```

```
##                           JD_06232014_sample1-A.raw
## JD_06232014_sample1-A.raw                 1.0000000
## JD_06232014_sample1_B.raw                 0.9794954
## JD_06232014_sample4-A.raw                 0.9502142
##                           JD_06232014_sample1_B.raw
## JD_06232014_sample1-A.raw                 0.9794954
## JD_06232014_sample1_B.raw                 1.0000000
## JD_06232014_sample4-A.raw                 0.9502517
##                           JD_06232014_sample4-A.raw
## JD_06232014_sample1-A.raw                 0.9502142
## JD_06232014_sample1_B.raw                 0.9502517
## JD_06232014_sample4-A.raw                 1.0000000
```

### A note on correlation and replication

It is often assumed that high correlation is a halmark of good
replication. Rather than focus on the correlation of the data, a
better measurement would be to look a the log2 fold-changes, i.e. the
distance between or repeated measurements. The ideal way to visualise
this is on an MA-plot:



```r
par(mfrow = c(1, 2))
r1 <- x2[, 1]
r2 <- x2[, 2]
M <- r1 - r2
A <- (r1 + r2)/2
plot(A, M); grid()
suppressPackageStartupMessages(library("affy"))
affy::ma.plot(A, M)
```

![plot of chunk unnamed-chunk-29](figure/unnamed-chunk-29-1.png)

See also this
[post](http://simplystatistics.org/2015/08/12/correlation-is-not-a-measure-of-reproducibility/)
on the *Simply Statistics* blog.

## Linear modelling

On the plots above, `abline(0, 1)` was used to add a line with
intercept 0 and slop 1. It we want to add the line that models the
data linearly, we can calculate the parameters using the `lm`
function:


```r
lmod <- lm(r2 ~ r1)
lmod
```

```
## 
## Call:
## lm(formula = r2 ~ r1)
## 
## Coefficients:
## (Intercept)           r1  
##      0.3482       0.9859
```

which can be used to add the adequate line that reflects the (linear)
relationship between the two data


```r
plot(r1, r2)
abline(lmod, col = "red")
```

![plot of chunk unnamed-chunk-31](figure/unnamed-chunk-31-1.png)

As we have seen in the beginning of this section, it is essential not
to rely solely on the correlation value, but look at the data. This
also holds true for linear (or any) modelling, which can be done by
plotting the model:


```r
par(mfrow = c(2, 2))
plot(lmod)
```

![plot of chunk unnamed-chunk-32](figure/unnamed-chunk-32-1.png)

* *Cook's distance* is a commonly used estimate of the influence of a
  data point when performing a least-squares regression analysis and
  can be used to highlight points that particularly influence the
  regression.
  
* *Leverage* quantifies the influence of a given observation on the
  regression due to its location in the space of the inputs.

See also `?influence.measures`.


> **Challenge**
> 
> 1. Take any of the `iprg2` replicates, model and plot their linear
>    relationship. The `iprg2` data is available as an `rda` file, or
>    regenerate it as shown above.
> 2. The Anscombe quartet is available as `anscombe`. Load it, create
>    a linear model for one $(x_i, y_i)$ pair of your choice and
>    visualise/check the model.


--- 
Back to course [home page](https://github.com/MayInstitute/MayInstitute2017/blob/master/Program3_Intro%20stat%20in%20R/README.md)
