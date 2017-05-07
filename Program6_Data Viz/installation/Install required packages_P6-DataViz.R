########################################################################################
## Install required packages for Advanced R course section from CRAN and Bioconductor
required_packages <- c('reshape2','ggplot2','dplyr','gplots','RColorBrewer', 'limma', 'Biobase', 'bioDist', 'genefilter', 'VennDiagram', 'UpSetR' )

source("https://bioconductor.org/biocLite.R")
biocLite(required_packages)

########################################################################################
## Tests whether all required packages 
## Press the 'Source' button in the top right corner of this pane and check 
## whether the output in the Console pane confirms that all packages are installed

required_packages <- c(required_packages)

installed_packages <- required_packages %in% installed.packages()[,"Package"]
missing_packages <- required_packages[!installed_packages]
if ( length(missing_packages) > 0 ) {
  cat(sprintf('FOLLOWING PACKAGES NEED TO BE INSTALLED STILL:\n\t%s\n',
              paste(missing_packages, collapse=', ')))
} else{
  cat('ALL PACKAGES ARE INSTALLED, WE\'RE GOOD TO GO!\n')
}

