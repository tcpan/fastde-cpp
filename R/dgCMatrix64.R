library(Matrix)

# TODO: [] limit row and col dimensions to int.     need for long per axis?

#' modified dgCMatrix class (definitions.R) to support dimnames
#' 
#' added Dimnames, list of 2 char arrays.
#' @slot Dimnames list of 2 char arrays for the row and col names.
#' @importClassesFrom Matrix dsparseMatrix
#' @export
#' @export dgCMatrix64
#' @exportClass dgCMatrix64
dgCMatrix64 <- setClass("dgCMatrix64", contains = c("dsparseMatrix"),
	 slots = c(i = "numeric", p = "numeric"),
	 prototype = prototype(p = 0)
     )
# no validation.  use Dim - assume dimensions are at most 2B..

# Internal function to create new dgCMatrix64 or dgCMatrix64 object, depending on the size.
.newdgCMatrix64 <- function(x=NULL,
                     i=NULL,
                     p=0,
                     dim=c(0, 0),
                     dimnames=c(NULL, NULL)
                     )
{
    newx <- new("dgCMatrix64")
    slot(newx,"x",check=FALSE) <- x
    x <- NULL
    slot(newx,"i",check=FALSE) <- as.numeric(i)
    i <- NULL
    slot(newx,"p",check=FALSE) <- as.numeric(p)
    p <- NULL
    slot(newx,"Dim",check=FALSE) <- as.integer(dim)
    dim <- NULL
    slot(newx,"Dimnames", check=FALSE) <- dimnames
    dimnames <- NULL
    
    return(newx)
}

#' @export
"is.dgCMatrix64" <- function(x) is(x,"dgCMatrix64")

#' @exportMethod as.dgCMatrix64
"as.dgCMatrix64" <- function(x, eps = .Machine$double.eps)
    stop("coercion not defined form this class")


#' @export
as.dgCMatrix64.dgCMatrix64 <- function(x, eps = .Machine$double.eps)  {
    newx <- new("dgCMatrix64")
    slot(newx, "x", check=FALSE) <- x@x
    slot(newx, "i", check=FALSE) <- x@i
    slot(newx, "p", check=FALSE) <- x@p
    slot(newx, "Dim", check=FALSE) <- x@Dim
    slot(newx, "Dimnames", check=FALSE) <- x@Dimnames
    return(newx)
}

#' @export
as.dgCMatrix64.dgCMatrix <- function(x, eps = .Machine$double.eps)  {
    newx <- new("dgCMatrix64")
    slot(newx, "x", check=FALSE) <- x@x
    slot(newx, "i", check=FALSE) <- as.numeric(x@i)
    slot(newx, "p", check=FALSE) <- as.numeric(x@p)
    slot(newx, "Dim", check=FALSE) <- x@Dim
    slot(newx, "Dimnames", check=FALSE) <- x@Dimnames
    return(newx)
}

# do this in c++ (sparsify)
##' @export 
#as.dgCMatrix64.matrix <- function(x, eps = .Machine$double.eps) {

setGeneric("as.dgCMatrix64")
setMethod("as.dgCMatrix64","dgCMatrix64", as.dgCMatrix64.dgCMatrix64)
setMethod("as.dgCMatrix64","dgCMatrix", as.dgCMatrix64.dgCMatrix)
#setMethod("as.dgCMatrix64","matrix", as.dgCMatrix64.matrix)

