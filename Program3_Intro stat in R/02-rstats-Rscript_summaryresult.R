### make summaryresult table from  Day 2

load('./data/iprg.rda')

# Let's start with one protein, named "sp|P44015|VAC2_YEAST"
oneproteindata <- iprg[iprg$Protein == "sp|P44015|VAC2_YEAST", ]

# there are 12 rows in oneproteindata
oneproteindata

### Calculate mean per groups
## splits 'oneproteindata' into subsets by 'Condition', 
## then, compute 'FUN=mean' of 'log2Int'
sub.mean <- aggregate(Log2Intensity ~ Condition,
                      data = oneproteindata,
                      FUN = mean)
sub.mean

## The same as mean calculation above. 'FUN' is changed to 'sd'.
sub.sd <- aggregate(Log2Intensity ~ Condition,
                    data = oneproteindata, FUN = sd)
sub.sd

## The same as mean calculation. 'FUN' is changed 'length'.
sub.len <- aggregate(Log2Intensity ~ Condition,
                     data = oneproteindata,
                     FUN = length)
sub.len

sub.se <- sqrt(sub.sd$Log2Intensity^2 / sub.len$Log2Intensity)
sub.se

## paste0 : concatenate vectors after convering to character.
grp <- paste0("Condition", 1:4)
## It is equivalent to paste("Condition", 1:4, sep="")
summaryresult <- data.frame(Group = grp,
                            mean = sub.mean$Log2Intensity,
                            sd = sub.sd$Log2Intensity, 
                            se = sub.se, 
                            length = sub.len$Log2Intensity)
summaryresult

summaryresult$ciw.lower.95 <- summaryresult$mean -
    qt(0.975, summaryresult$len - 1) * summaryresult$se
summaryresult$ciw.upper.95 <- summaryresult$mean +
    qt(0.975, summaryresult$len - 1) * summaryresult$se
summaryresult

save(summaryresult, file='summaryresult.rda')

