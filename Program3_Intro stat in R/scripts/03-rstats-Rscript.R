# Part 4: Statistical hypothesis test

load("./data/summaryresults.rda")
load("./data/iprg.rda")


## Two sample t-test for one protein with one feature

### Data preparation
## Let's start with one protein, named "sp|P44015|VAC2_YEAST"
oneproteindata <- iprg[iprg$Protein == "sp|P44015|VAC2_YEAST", ]

## Then, get two conditions only, because t.test only works for two groups (conditions).
oneproteindata.condition12 <- oneproteindata[oneproteindata$Condition %in% 
                                                 c('Condition1', 'Condition2'), ]
oneproteindata.condition12
table(oneproteindata.condition12[, c("Condition", "BioReplicate")])

## remove duplicated
oneproteindata.condition12 <- droplevels(oneproteindata.condition12)
table(oneproteindata.condition12[, c("Condition", "BioReplicate")])


?t.test

# t test for different abundance (log2Int) between Groups (Condition)
result <- t.test(Log2Intensity ~ Condition,
                 data = oneproteindata.condition12,
                 var.equal = FALSE)

result

## **Challenge**
## Repeat the t-test above but with calculating a 90% confidence interval for the log2 fold change.


### The htest class
class(result)

names(result)

result$estimate

## **Challenge**
## * Calculate the (log2-transformed) fold change between groups
## * Extract the value of the t-statistics
## * Calculate the standard error (fold-change/t-statistics)
## * Extract the degrees of freedom (parameter)
## * Extract the p values
## * Extract the 95% confidence intervals


summaryresult
summaryresult12 <- summaryresult[1:2, ]

## test statistic, It is the same as 'result$statistic' above.
diff(summaryresult12$mean) ## different sign, but absolute values are same as result$estimate[1]-result$estimate[2]
sqrt(sum(summaryresult12$sd^2/summaryresult12$length)) ## same as stand error

## the t-statistic : sign is different
diff(summaryresult12$mean)/sqrt(sum(summaryresult12$sd^2/summaryresult12$length))

# Part 5a: Sample size calculation

library("pwr")
## ?pwr.t.test
# Significance level alpha
alpha <- 0.05

# Power = 1 - beta
power <- 0.95

# anticipated log2 fold change 
delta <- 1

# anticipated variability
sigma <- 0.9

# Effect size
# It quantifies the size of the difference between two groups
d <- delta/sigma

#Sample size estimation
pwr.t.test(d = d, sig.level = alpha, power = power, type = 'two.sample')

## **Challenge**
## * Calculate power with 10 samples and the same parameters as above.

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
                             sig.level = alpha, 
                             power = power,
                             type = "two.sample")
        counter <- counter + 1
        samsize[counter,1] <- delta[i]
        samsize[counter,2] <- sigma[j]
        samsize[counter,3] <- ceiling(result$n)
    }
}

colnames(samsize) <- c("desiredlog2FC","variability","samplesize")

library("ggplot2")
samsize <- as.data.frame(samsize)
samsize$variability <- as.factor(samsize$variability)
ggplot(data=samsize, aes(x=desiredlog2FC, y=samplesize, group = variability, colour = variability)) +
    geom_line() +
    geom_point(size=2, shape=21, fill="white") +
    labs(title="Significance level=0.05, Power=0.95", x="Anticipated log2 fold change", y='Sample Size (n)') +
    theme(plot.title = element_text(size=20, colour="darkblue"),
          axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15),
          axis.text.x = element_text(size=13)) 

# Part 5c: Analysis of categorical data

TCGA.CRC <- read.csv("./data/TCGA_sample_information.csv")
head(TCGA.CRC)

## check the dimension of data
dim(TCGA.CRC)
dim(TCGA.CRC$Gender) # dim function is for matrix, array or data frame.

## check unique information of Gender column.
unique(TCGA.CRC$Gender)
class(TCGA.CRC$Gender)

## **Challenge**
#* Get unique information and class for Cancer information
#* Get unique information and class for BRAF.mutation information
#* Get unique information and class for history_of_colon_polyps information

## check how many female and male are in the dataset
table(TCGA.CRC$Gender)

## check unique information if paticipant ID
unique(TCGA.CRC$TCGA.participant.ID)
length(unique(TCGA.CRC$TCGA.participant.ID))

countID <- table(TCGA.CRC$TCGA.participant.ID)
countID

unique(countID)
countID[countID > 1]

TCGA.CRC[TCGA.CRC$TCGA.participant.ID == 'TCGA-AA-A00A', ]
TCGA.CRC <- TCGA.CRC[!duplicated(TCGA.CRC), ]

## **Challenge**
## ** Check whether dimension and number of participants ID are changed after removing duplicated rows.

## check the dimension of data
dim(TCGA.CRC)

## check unique information if paticipant ID
unique(TCGA.CRC$TCGA.participant.ID)
length(unique(TCGA.CRC$TCGA.participant.ID))

countID <- table(TCGA.CRC$TCGA.participant.ID)
countID

unique(countID)

## Generate 2-way contingency tables
cancer.polyps <- table(TCGA.CRC$history_of_colon_polyps, TCGA.CRC$Cancer)
cancer.polyps
dotchart(cancer.polyps, xlab = "Observed counts")

## Comparison of two proportions
polyps <- cancer.polyps[2, ]
cancer <- apply(cancer.polyps, 2, sum)
pt <- prop.test(polyps, cancer)
pt

# name of output
names(pt)

# proportion in each group
pt$estimate 

# test statistic value
pt$statistic 

# degree of freedom
pt$parameter

## Test of independence
### Pearson Chi-squared test
chisq.test(cancer.polyps)

### Fisher's exact test
ft <- fisher.test(cancer.polyps)
ft
ft$estimate 

## **Challenge**
## * Compare the proportion of male patients in the patients with colon cancer is different from that in the patients with rectal cancer.


## **Challenge**
# Write a function named f2c (c2f) that converts a temperature
# from Fahrenheit to Celsium (Celsium to Fahrenheit) using the
# following formula $F = C \times 1.8 + 32$ ($C = \frac{F - 32}{1.8}$).

# Part 6: Linear models and correlation

x <- iprg2[, c(1, 2, 10)]
head(x)
pairs(x)

cor(x)

x2 <- na.omit(x)
cor(x2)

## Linear modelling

lmod <- lm(r2 ~ r1)
summary(lmod)

plot(r1, r2)
abline(lmod, col = "red")

par(mfrow = c(2, 2))
plot(lmod)

# **Challenge**
# 1. Take any of the iprg2 replicates, model and plot their linear
#    relationship. The iprg2 data is available as an rda file, or
#    regenerate it as shown above.
# 2. The Anscombe quartet is available as anscombe. Load it, create
#    a linear model for one $(x_i, y_i)$ pair of your choice and
#    visualise/check the model.

library("ggplot2")
dfr <- data.frame(r1, r2, M, A)
p <- ggplot(aes(x = r1, y = r2), data = dfr) + geom_point()
p + geom_smooth(method = "lm") +
    geom_quantile(colour = "red")

# **Challenge**
# Replicate the MA plot above using ggplot2. Then add a
# non-parametric lowess regression using geom_smooth().

