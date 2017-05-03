
#' Proteomics Experiment
#'
#' Create a \code{ProtExp} object.
#'
#' \code{ProtExp} is a simple S3 class inheriting from \code{data.frame}
#' for storing proteomics experiments. It expects columns named
#' `protein', `feature', `run', `intensity', and `label'. Additional
#' columns are also allowed. A \code{ProtExp} object also has a
#' \code{is_log_exp} attribute for tracking whether the intensities
#' have been log2-transformed or not.
#'
#' @param protein A \code{character} vector of protein names.
#' @param feature A \code{character} vector of MS features.
#' @param run A \code{character} or \code{numeric} vector of MS runs.
#' @param intensity A \code{numeric} vector of MS intensities.
#' @param label A \code{character} or \code{factor} of MS labels.
#' @param ... Additional columns.
#' @param is_log_trans \code{TRUE} or \code{FALSE} indicating whether
#' the intensities has been log2-transformed or not.
#'
#' @return object A \code{ProtExp} object.
#'
#' @export
ProtExp <- function(protein, feature, run,
							intensity, label, ..., is_log_trans = FALSE)
{
	dots <- list(...)
	if ( length(dots) > 0 ) {
		object <- data.frame(protein=protein, feature=feature,
			run=run, intensity=intensity, label=label, list(...)) 
	} else {
		object <- data.frame(protein=protein, feature=feature,
			run=run, intensity=intensity, label=label) 
	}
	attr(object, "is_log_trans") <- is_log_trans
	class(object) <- c("ProtExp", "data.frame")
	object
}

#' Is a dataset log-transformed?
#'
#' Returns \code{TRUE} if the data has been log-transformed
#' and \code{FALSE} otherwise.
#'
#' This is an S3 generic function: methods can be defined for
#' it by writing functions with the naming convention
#' \code{is_log_trans.classname}. Functions should return
#' \code{TRUE} or \code{FALSE}. The default method gets and sets
#' the \code{is_log_trans} attribute.
#'
#' @param object An object with data that may have been log-transformed.
#' @param value \code{TRUE} is log-transformed, \code{FALSE} otherwise.
#'
#' @return \code{TRUE} if the data has been log-transformed
#' and \code{FALSE} otherwise.
#'
#' @examples
#' x <- structure("test", is_log_trans = TRUE)
#' is_log_trans(x)
#' is_log_trans(x) <- FALSE
#' is_log_trans(x)
#'
#' @rdname is_log_trans
#' @export
is_log_trans <- function(object) {
	UseMethod("is_log_trans")
}
#' @rdname is_log_trans
#' @export
`is_log_trans<-` <- function(object, value) {
	UseMethod("is_log_trans<-")
}
#' @rdname is_log_trans
#' @export
is_log_trans.default <- function(object) {
	attr(object, "is_log_trans")
}
#' @rdname is_log_trans
#' @export
`is_log_trans<-.default` <- function(object, value) {
	attr(object, "is_log_trans") <- value
	object
}


#' Normalize a dataset
#'
#' Performs normalization on an experimental dataset.
#'
#' This is an S3 generic function: methods can be defined for
#' it by writing functions with the naming convention
#' \code{normalize.classname}. Functions should be appropriate
#' for the experiment type.
#'
#' \code{normalize.ProtExp} performs normalization on a \code{ProtExp}
#' proteomics experiment dataset object using the median-of-medians method.
#' It also performs log2-transformation on the intensities if the data
#' has not already been log2-transformed.
#'
#' @param object An object with data to be normalized.
#' @param by The name of the label to use as the standard.
#' @param ... Ignored.
#'
#' @return Another object of the same class with the normalized data.
#'
#' @rdname normalize
#' @export
normalize <- function(object, ...) UseMethod("normalize")
#' @rdname normalize
#' @export
normalize.ProtExp <- function(object, ..., by) {
	if ( !is_log_trans(object) ) {
		object$intensity <- log2(object$intensity)
		is_log_trans(object) <- TRUE
	}
	wh <- which(object$label == by)
	medians <- tapply(object$intensity[wh], object$run[wh],
		median, na.rm=TRUE)
	gbl_median <- median(medians, na.rm=TRUE)
	object$intensity <- object$intensity - medians[object$run] + gbl_median
	object
}

#' @export
summary.ProtExp <- function(object, ...) {
	if ( is_log_trans(object) )
		object$intensity <- 2^object$intensity
	tapply(object$intensity,
		list(run=object$run,
			protein=object$protein,
			label=object$label),
		function(x) log2(sum(x, na.rm=TRUE)))
}


#' @export
plot.ProtExp <- function(x, y, ...) {
	if ( is_log_trans(x) ) {
		ylab <- "log2(intensity)"
	} else {
		ylab <- "intensity"
	}
	boxplot(x$intensity ~ x$run, xlab="run", ylab=ylab)
}

#' Twin DIA and SRM experiments
#'
#' DIA and SRM datasets from the same experiment. These datasets
#' are subsets of the original datasets and have 33 proteins by DIA,
#' and 39 proteins by SRM. The dataset includes 58 pairs of monozygotic (MZ)
#' and dizygotic (DZ) twins, each measured at 2 time points, yielding
#' 58 x 2 x 2 = 232 MS runs. The variables are as follows:
#'
#' \itemize{
#'   \item protein. (chr) protein name.
#'   \item feature. (chr) combination of peptide, precursor charge state,
#'			fragment ion, and product charge state, separated by _.
#'   \item run.  (chr) MS run identifier (R001-R232).
#'   \item pair. (int) pair identifier number (1-58).
#'   \item zygosity. (factor) zygosity (MZ, DZ).
#'   \item subject. (int) subject identifier number (1-116).
#'   \item visit. (int) time of visit (1, 2).
#'   \item intensity_h. (num) integrated feature intensity
#'						from light (L) channel.
#'   \item intensity_l. (num) integrated feature intensity
#'						from heavy (H, aka reference) channel.
#' }
#'
#' @name twin
#' @format Two data frames with 103936 rows and 9 variables for
#' 			\code{twin_dia} and 63104 rows and 9 variables for
#'			\code{twin_srm}.
#' @examples
#' data(twin)
#' head(twin_dia)
#' head(twin_srm)
NULL

