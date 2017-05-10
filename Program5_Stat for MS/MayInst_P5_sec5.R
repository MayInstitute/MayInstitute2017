
# 1. Protein quantification from the output of `dataProcess` in MSstats
## 1.1 Read the output of `dataProcess` 
load('./data_day3/quant.skyline.rda')

head(quant.skyline$RunlevelData)
runquant <- quant.skyline$RunlevelData
head(runquant)

## 1.2 Reformat to wide format
library(reshape2)

input.wide <- dcast(Protein ~ originalRUN, data=runquant, value.var = 'LogIntensities')
head(input.wide)
rownames(input.wide) <- input.wide$Protein
input.wide <- input.wide[, -1]
head(input.wide)
dim(input.wide) # there are 3027 rows (one row per protein) and 12 columns (each colomn for run)

sum(is.na(input.wide))

samplemissing <- apply(input.wide, 1, function(x) sum(is.na(x)))
unique(samplemissing)
selectedSamples <- which(samplemissing == 0) 
# Only keep the samples with no missing values
input.wide <- input.wide[selectedSamples, ]
dim(input.wide)

input.wide.spike <- input.wide[rownames(input.wide) %in% c("sp|P44015|VAC2_YEAST",
                                                           "sp|P55752|ISCB_YEAST",
                                                           "sp|P44374|SFG2_YEAST",
                                                           "sp|P44983|UTR6_YEAST",
                                                           "sp|P44683|PGA4_YEAST",
                                                           "sp|P55249|ZRT4_YEAST"),]

save(input.wide.spike, file = 'input.wide.spike.rda')
save(input.wide.spike, file = 'input.wide.spike.rda')


# 2. Perform principal components analysis (PCA) on the given data matrix

## 2.1 Only 6 spike-in proteins

### 2.1.1 PCA with `prcomp` function
head(input.wide.spike)
pc <- prcomp(t(input.wide.spike))

# Inspect PCA object
summary(pc)
names(pc)

### 2.1.2 Check the proportion of explained variance

percent_var <- pc$sdev^2/sum(pc$sdev^2)
barplot(percent_var, xlab="Principle component", ylab="% of variance")

cum_var <- cumsum(pc$sdev^2/sum(pc$sdev^2))
barplot(cum_var, xlab="Principle component", ylab="Cumulative % of variance" )

### 2.1.3 Visualization for PC1 vs PC2

head(pc$x)

library(ggplot2)

ggplot(aes(x=PC1, y=PC2), data=data.frame(pc$x))+
    geom_point(size=4, alpha=0.5)+
    theme_bw()

colnames(input.wide.spike)
Condition <- c(rep("Condition1", 3), rep("Condition2", 3),
               rep("Condition3", 3), rep("Condition4", 3))
Condition

# Create PC1 vs PC2 scatterplot with Condition colors
ggplot(aes(x=PC1, y=PC2, color=Group), data=data.frame(pc$x, Group=Condition))+
    geom_point(size=4, alpha=0.5)+
    theme_bw()

ggplot(aes(x=PC1, y=PC2, color=Group, shape=Group), data=data.frame(pc$x, Group=Condition))+
    geom_point(size=4, alpha=0.5)+
    theme_bw()

ggplot(aes(x=PC1, y=PC3, color=Group, shape=Group), data=data.frame(pc$x, Group=Condition))+
    geom_point(size=4, alpha=0.5)+
    theme_bw()

ggplot(aes(x=PC2, y=PC3, color=Group, shape=Group), data=data.frame(pc$x, Group=Condition))+
    geom_point(size=4, alpha=0.5)+
    theme_bw()

### 2.1.4 Check the PC1 loadings

head(pc$rotation)

ggplot(aes(x=Protein, xend=Protein, y=0, yend=PC1), 
data=data.frame(pc$rotation, Protein=rownames(pc$rotation)))+
geom_segment()+
labs(y="PC1 loading")+
theme_bw()+
theme(axis.text.x = element_text(angle=90))

### 2.1.5 loadings vs p-values

load('./data_day3/Skyline.result.rda')
head(Skyline.result)

Skyline.result.spikein <- Skyline.result[Skyline.result$Protein %in% c("sp|P44015|VAC2_YEAST",
                                                                       "sp|P55752|ISCB_YEAST",
                                                                       "sp|P44374|SFG2_YEAST",
                                                                       "sp|P44983|UTR6_YEAST",
                                                                       "sp|P44683|PGA4_YEAST",
                                                                       "sp|P55249|ZRT4_YEAST"), ]
Skyline.result.spikein

smoothScatter(pc$rotation[,1], Skyline.result.spikein[Skyline.result.spikein$Label == 'C2-C1', 'adj.pvalue'])
smoothScatter(pc$rotation[,1], Skyline.result.spikein[Skyline.result.spikein$Label == 'C3-C1', 'adj.pvalue'])


## 2.2 All 3025 proteins

### 2.2.1 PCA with `prcomp` function

head(input.wide)
dim(input.wide)

pc <- prcomp(t(input.wide))

# Inspect PCA object
summary(pc)
names(pc)

### 2.2.2 Check the proportion of explained variance

percent_var <- pc$sdev^2/sum(pc$sdev^2)
barplot(percent_var, xlab="Principle component", ylab="% of variance")

cum_var <- cumsum(pc$sdev^2/sum(pc$sdev^2))
barplot(cum_var, xlab="Principle component", ylab="Cumulative % of variance" )

### 2.1.3 Visualization for PC1 vs PC2
# 'x' include PC components for each subject.
head(pc$x)

ggplot(aes(x=PC1, y=PC2), data=data.frame(pc$x))+
geom_point(size=4, alpha=0.5)+
theme_bw()

colnames(input.wide)
Condition <- c(rep("Condition1", 3), rep("Condition2", 3),
rep("Condition3", 3), rep("Condition4", 3))
Condition

# Create PC1 vs PC2 scatterplot with Condition colors
ggplot(aes(x=PC1, y=PC2, color=Group), data=data.frame(pc$x, Group=Condition))+
geom_point(size=4, alpha=0.5)+
theme_bw()

# 3. Heatmap

## 3.1 Only 6 spike-in proteins

### 3.1.1 matrix format
# check the class
class(input.wide.spike)

# It is data.frame. Convert to numeric matrix
input.wide.spike <- as.matrix(input.wide.spike)
class(input.wide.spike)

# Visually no difference
head(input.wide.spike)

### 3.1.2 `heatmap` function in base `stats` package

heatmap(input.wide.spike, scale='none')
heatmap(input.wide.spike, scale='row')
heatmap(input.wide.spike, scale='column')

# Change the font of row and column label
heatmap(input.wide.spike, cexRow = 0.6, cexCol = 0.6, margins = c(2, 2))
# Don't do cluster on rows
heatmap(input.wide.spike, cexRow = 0.6, cexCol = 0.6, margins = c(2, 2), Rowv = NA)
# Don't do cluster on columns
heatmap(input.wide.spike, cexRow = 0.6, cexCol = 0.6, margins = c(2, 2), Colv = NA)

### 3.1.3 Color bar for group information

group.color <- c(rep("red", 3), rep("orange", 3),
                 rep("yellow", 3), rep("blue", 3))
heatmap(input.wide.spike, ColSideColors=group.color)


### 3.1.4 Color scale

library(marray)
my.colors <- c(maPalette(low = "darkblue", high = "white", k = 7)[-7],
               "white",maPalette(low = "white", high = "darkred", k = 7)[-1])
heatmap(input.wide.spike, 
        col=my.colors, 
        ColSideColors=group.color,
        scale='none')

### 3.1.5 Different distance and clustering

#* Distance options: euclidean (default), maximum, canberra, binary, minkowski, manhattan
#* Cluster options: complete (default), single, average, mcquitty, median, centroid, ward

# can change method for distance calculation
col_distance <- dist(t(input.wide.spike), method = "euclidean")
# can change clustering method
col_cluster <- hclust(col_distance, method = "ward.D")

heatmap(input.wide.spike, 
        col=my.colors, 
        ColSideColors=group.color,
        Colv = as.dendrogram(col_cluster),
        scale='none')


## 3.2 All proteins (3025 proteins)

### 3.2.1 matrix format
# check the class
class(input.wide)

# It is data.frame. Convert to numeric matrix
input.wide <- as.matrix(input.wide)
class(input.wide)

# Visually no difference
head(input.wide)

### 3.2.3 Color bar and non-scaled data

group.color <- c(rep("red", 3), rep("orange", 3),
                 rep("yellow", 3), rep("blue", 3))

# remove row information (Protein IDs)
heatmap(input.wide, 
        ColSideColors=group.color,
        col=my.colors,
        labRow = FALSE,
        scale = 'none')

# no Protein IDs
heatmap(input.wide, 
        ColSideColors=group.color,
        col=my.colors,
        labRow = FALSE,
        scale = 'row')

# can change method for distance calculation
col_distance <- dist(t(input.wide), method = "euclidean")
# can change clustering method
col_cluster <- hclust(col_distance, method = "ward.D")

heatmap(input.wide, 
        col=my.colors, 
        ColSideColors=group.color,
        Colv = as.dendrogram(col_cluster),
        labRow = FALSE,
        scale = 'row')

# 4. Biological data 
## 4.1 Read subject quantification per protein and handle missing values.

crc <- read.csv("./data_day3/CRC_train.csv")
head(crc)
colnames(crc)
dim(crc)

crc.prot <- crc[,1:72]
crc.anno <- crc[,73:79]

## check missing value
sum(is.na(crc.prot))

# Deal with missing value
# First option: remove the samples with missing values
dim(na.omit(crc.prot))

# Second option: impute the missing values
random.imp <- function (a){
    missing <- is.na(a)
    n.missing <- sum(missing)
    a.obs <- a[!missing]
    imputed <- a
    # imputed[missing] <- 0
    # imputed[missing] <- sample(a.obs, n.missing, replace=TRUE)
    # imputed[missing] <- median(a.obs)
    imputed[missing] <- min(a.obs)
    return (imputed)
}

pMiss <- function(x){
    sum(is.na(x))/length(x)*100
}

# 1. Only keep the proteins with less than 5% missing values
proteinmissing <- apply(crc.prot, 2, pMiss)
selectedprotein <- which(proteinmissing <= 5) 
crc.prot <- crc.prot[, selectedprotein]

# 2. Only keep the samples with less than 5% missing values
samplemissing <- apply(crc.prot, 1, pMiss)
selectedSamples <- which(samplemissing <= 5) 
# Then, impute the missing value
imputed.crc.prot <- t(apply(crc.prot[selectedSamples,], 1, function(x) random.imp(x)))
imputed.crc.anno <- crc.anno[selectedSamples, ]
dim(imputed.crc.anno)

## 4.2 Heatmap of individual values of proteins

### 4.2.1 Default heatmap
library(genefilter)
## heatmap function in base stats package
heatmap(t(imputed.crc.prot))

### 4.2.2 Change color

# Non-scaled 
heatmap(t(imputed.crc.prot), cexRow = 0.3, cexCol = 0.2, margins = c(2, 2), Colv = NA, scale='none')

# scaled
heatmap(t(imputed.crc.prot), cexRow = 0.3, cexCol = 0.2, margins = c(2, 2), Colv = NA)

heatmap(t(imputed.crc.prot), cexRow = 0.3, cexCol = 0.2, margins = c(2, 2), Colv = NA,
        col=my.colors)

### 4.2.3 Color bar for group information
unique(imputed.crc.anno$Sub_group)
myColor <- rep("blue", nrow(imputed.crc.prot))
myColor[imputed.crc.anno$Sub_group == "CRC"] <- "red" 
myColor[imputed.crc.anno$Sub_group == "Benign"] <- "yellow" 

### 4.2.4 Standardized in a row (a protein)
# change color code: protein-standardized expressions
# quantile-spaced breaks in positive and negative directions
imputed_CRC_prot_centerScale <-  ( t(imputed.crc.prot) - rowMeans(t(imputed.crc.prot)) ) / rowSds(t(imputed.crc.prot))

exprs_brk <- c( quantile(imputed_CRC_prot_centerScale[imputed_CRC_prot_centerScale < 0], probs=seq(0, 1, length=10)), 0, quantile(imputed_CRC_prot_centerScale[imputed_CRC_prot_centerScale > 0], seq(0, 1, length=10)))
colors <- colorRampPalette(c("darkblue", "white", "darkred"))(length(exprs_brk)-1)

heatmap(t(imputed.crc.prot), 
        cexRow = 0.3, cexCol = 0.2, margins = c(2, 2), Rowv = NA,
        col = colors, breaks = exprs_brk, ColSideColors=myColor)

### 4.2.5 Different distance and clustering
# can change method for distance calculation
col_distance <- dist(imputed.crc.prot, method = "euclidean")
# can change clustering method
col_cluster <- hclust(col_distance, method = "ward.D")

heatmap(t(imputed.crc.prot),
        cexRow = 0.3, cexCol = 0.2, margins = c(2, 2), 
        col = colors, 
        breaks = exprs_brk, 
        ColSideColors = myColor,
        Rowv = NA, 
        Colv = as.dendrogram(col_cluster)) 

