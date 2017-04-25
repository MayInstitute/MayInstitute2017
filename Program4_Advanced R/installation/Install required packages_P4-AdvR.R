########################################################################################
## Install required packages for Advanced R course section from CRAN and Bioconductor
required_packages_CRAN <- c("ggplot2", "reshape", "reshape2", "tidyverse", "pryr",
							"microbenchmark", "roxygen2", "RUnit", "testthat",
							"bigmemory", "ff", "foreach", "doParallel", "Rcpp")
required_packages_BioC <- c("MSstats", "limma","marray","preprocessCore",
							"MSnbase", "Cardinal", "matter", "HDF5Array")

install.packages(required_packages_CRAN)

source("http://bioconductor.org/biocLite.R")
biocLite(required_packages_BioC)

########################################################################################
## Tests whether all required packages 
## Press the 'Source' button in the top right corner of this pane and check 
## whether the output in the Console pane confirms that all packages are installed

required_packages <- c(required_packages_CRAN, required_packages_BioC)

installed_packages <- required_packages %in% installed.packages()[,"Package"]
missing_packages <- required_packages[!installed_packages]
if ( length(missing_packages) > 0 ) {
	cat(sprintf('FOLLOWING PACKAGES NEED TO BE INSTALLED STILL:\n\t%s\n',
		paste(missing_packages, collapse=', ')))
} else{
	cat('ALL PACKAGES ARE INSTALLED, WE\'RE GOOD TO GO!\n')
}

## NOTE: To run the examples that will be presented in the "Scalability" section:
## On Windows, please install Rtools (https://cran.r-project.org/bin/windows/Rtools/)
## On Mac, please install the Xcode Developer Tools (only Command Line Tools are necessary)
