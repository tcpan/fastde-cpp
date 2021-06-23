#' fastde
#' 
#' Description of your package
#' 
#' @docType package
#' @author Tony Pan <tpan7@gatech.edu>
#' @import Rcpp tictoc Seurat BioQC
#' @importFrom Rcpp evalCpp
#' @import future
#' @useDynLib fastde
#' @name fastde
NULL  

# roxygen2 is needed to manage exports.
# follow https://gallery.rcpp.org/articles/documenting-rcpp-packages/
# and https://www.r-bloggers.com/2016/08/rcpp-and-roxygen2/
# each comment of funcs being documented has to have @name 
# @rdname xyz means xyz.Rd has to exist. so skip 1st use of rdname.