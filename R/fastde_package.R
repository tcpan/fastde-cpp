#' fastde
#' 
#' Description of your package
#' 
#' @docType package
#' @author Tony Pan <tcpan@emory.edu>
#' @import tictoc 
#' @import Seurat 
#' @import Matrix 
#' @import rhdf5 
#' @import BioQC
#' @import methods
#' @import future
#' @useDynLib fastde, .registration = TRUE
#' @name fastde
NULL  

# roxygen2 is needed to manage exports.
# follow https://gallery.rcpp.org/articles/documenting-rcpp-packages/
# and https://www.r-bloggers.com/2016/08/rcpp-and-roxygen2/
# each comment of funcs being documented has to have @name 
# @rdname xyz means xyz.Rd has to exist. so skip 1st use of rdname.