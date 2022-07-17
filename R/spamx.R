library(Matrix)
library(spam)
library(spam64)



#' modified spam class (definitions.R) to support dimnames
#' 
#' added Dimnames, list of 2 char arrays.
#' @slot Dimnames list of 2 char arrays for the row and col names.
#' @importClassesFrom spam spam
#' @importClassesFrom Matrix Matrix
#' @export
#' @export spamx
#' @exportClass spamx
spamx <- setClass("spamx", contains=c("spam", "Matrix"))


#------------ helper, from spam/rep_len64.R

rep_len64 <- function(x, length.out, NAOK = getOption("spam.NAOK")){
    if(getOption("spam.force64") || length.out > 2147483647){
        .format64()
        return(.C64("rep_len64_c",
                    SIGNATURE = c("double", "int64", "int64", "double"),
                    
                    x = x,
                    lx = length(x),
                    length = length.out,
                    
                    out = numeric_dc(length.out),
                    
                    INTENT = c("r","r","r","w"),
                    NAOK = NAOK,
                    PACKAGE = "spam64"
                    )$out)
    }
    return(rep_len(x = x, length.out = length.out))
}


#------------ helpers, from spam/profile.R
.format32 <- list(
    name      = "32-bit",
    type      = "integer", 
    signature = "integer",  
    package   = "spam")

.format64 <- list(
    name      = "64-bit",
    type      = "numeric",
    signature = "int64",
    package   = "spam64")

.format64 <- function(){
    if (!isNamespaceLoaded("spam64")) {
        stop("Large (64-bit) sparse matrices detected. Please load the required package 'spam64' and see the help page '?large_matrix'.")
    }
    list(
        name      = "64-bit",
        type      = "numeric",
        signature = "int64",
        package   = "spam64")
}

.format.spam <- function(x, ... , validate = getOption("spam.validate") ) {
    objects <- c(list(x), list(...))

    if (validate) for(o in objects)  stopifnot(spam::validate_spam(o))
    
    for(o in objects){
        ## If both pointer vectors are of the same type,
        ## use this type to determine the format
        if(identical(typeof(o@colindices), typeof(o@rowpointers))) {
            if(identical(typeof(o@colindices), "double")){
                return(.format64())
            }
            next
        }
        
        ## As fallback use the length of the entries vector and the dimension
        if(nrow(o) > 2147483647 || ncol(o) > 2147483647 || 
           length(o@entries) > 2147483647){
            return(.format64())
        }
    }   
    return(.format32)
}
#------------ END helpers, from spam/profile.R


# Internal function to create new spamx or spamx64 object, depending on the size.
.newSpamx <- function(entries=0,
                     colindices=1,
                     rowpointers=NULL,
                     dimension=c(0, 0),
                     Dim=c(0, 0),
                     dimnames=c(NULL, NULL),
                     force64=getOption("spam.force64")
                     )
{
    if(is.null(rowpointers)) {
        rowpointers <- c(1, rep_len64(2,dimension[1]))
    }

    if(force64 || length(colindices) > 2147483647 ||
       length(colindices) > 2147483647 || max(dimension) > 2147483647 ){
        newx <- new("spamx")
        slot(newx,"entries",check=FALSE) <- entries
        entries <- NULL
        slot(newx,"colindices",check=FALSE) <- as.numeric(colindices)
        colindices <- NULL
        slot(newx,"rowpointers",check=FALSE) <- as.numeric(rowpointers)
        rowpointers <- NULL
        slot(newx,"dimension",check=FALSE) <- as.numeric(dimension)
        slot(newx,"Dim",check=FALSE) <- as.integer(c(-1, -1))
    } else {
        newx <- new("spamx")
        slot(newx,"entries",check=FALSE) <- entries
        entries <- NULL
        slot(newx,"colindices",check=FALSE) <- as.integer(colindices)
        colindices <- NULL
        slot(newx,"rowpointers",check=FALSE) <- as.integer(rowpointers)
        rowpointers <- NULL
        slot(newx,"dimension",check=FALSE) <- as.integer(dimension)
        slot(newx,"Dim",check=FALSE) <- as.integer(dimension)
    }

    slot(newx, "Dimnames", check=FALSE) <- dimnames
    return(newx)
}

# skip 
# .print.spam, print_nnzpos, print.spam, show, print, length<-,  length, c, diag, t, as.matrix, as.vector
# accesor functions [], assign, %*%, upper.tri, lower.tri, %+%, all.equal, isSymmetric

#' @export
"is.spamx" <- function(x) is(x,"spamx")

#' @exportMethod as.spamx
"as.spamx" <- function(x,  eps = getOption("spam.eps"))
    stop("coercion not defined form this class")

#' @exportMethod spamx
"spamx" <- function(x, nrow = 1, ncol = 1, eps = getOption("spam.eps"))
    stop("argument 'x' should be of mode 'numeric' (or 'spamx')")

#' @export
as.spamx.spam <- function(x, eps = getOption("spam.eps")) {
    force64 <- getOption("spam.force64")

    if(force64)
        SS <- .format64()
    else
        SS <- .format.spam(x)

    if (eps < .Machine$double.eps)
        stop("'eps' should not be smaller than machine precision", call.=FALSE)
    dimx <- x@dimension

    z <- .C64("cleanspam",
              SIGNATURE = c(SS$signature, "double", SS$signature, SS$signature, "double"),

              nrow=dimx[1],
              entries=x@entries,
              colindices=x@colindices,
              rowpointers=x@rowpointers,
              eps=eps,

              INTENT=c("r", "rw", "rw", "rw", "r"),
              NAOK = getOption("spam.NAOK"),
              PACKAGE = SS$package
              )
    nz <- z$rowpointers[dimx[1]+1]-1
    if(nz==0) {
        return(.newSpamx(
            dimension=dimx,
            force64 = force64
        ))
    }

    return(.newSpamx(
        entries=z$entries[1:nz],
        colindices=z$colindices[1:nz],
        rowpointers=z$rowpointers[1:(dimx[1]+1)],
        dimension=dimx,
        force64=force64
    ))
}

#' @export
as.spamx.spamx <- function(x, eps = getOption("spam.eps"))  {
    newx <- as.spamx.spam(x, eps)
    slot(newx, "Dimnames", check=FALSE) <- x@Dimnames
    return(newx)

}

#' @export
"cleanup" <- function(x, eps = getOption("spam.eps")) {
  if (is.spamx(x)) as.spamx.spamx(x,eps) 
  else if (is.spam(x)) as.spamx.spam(x, eps)
  else x
}

#' @export
as.spamx.matrix <- function(x, eps = getOption("spam.eps")) {
    force64 <- getOption("spam.force64")
    if (eps<.Machine$double.eps) stop("'eps' should not be smaller than machine precision",call.=FALSE)

    dimx <- dim(x)
    nonz <- abs(x) > eps
    nz <- sum(nonz, na.rm = TRUE) + sum(!is.finite(nonz))

    #only zero entries
    if(nz==0) {
        return(.newSpamx(
            dimension=dimx,
            force64=force64
        ))
    }

    # EXPLICIT-STORAGE-FORMAT: Depending on the length of x, use 32 or 64-bit:
    SS <- if( force64 || nz > 2147483647 || max(dimx) > 2147483647 ){
              .format64()
          }else{
              .format32
          }

    z <- .C64("spamdnscsr",
              SIGNATURE=c(SS$signature, SS$signature, "double", SS$signature,
                          "double", SS$signature, SS$signature, "double"),

              nrow = dimx[1],
              ncol = dimx[2],
              x = x,
              dimx[1],

              entries = vector_dc("double",nz),
              colindices = vector_dc(SS$type,nz),
              rowpointers = vector_dc(SS$type,dimx[1]+1),
              eps = eps,

              INTENT=c("r", "r", "r", "r",
                       "w", "w", "w", "r"),
              NAOK = getOption("spam.NAOK"),
              PACKAGE = SS$package
              )

    return(.newSpamx( entries=z$entries,
                    colindices=z$colindices,
                    rowpointers=z$rowpointers,
                    dimension=dimx,
                    dimnames=list(rownames(x), colnames(x)),
                    force64 = force64))
}

#' @export
"as.spamx.numeric" <- function(x, eps = getOption("spam.eps")) {
    if (eps<.Machine$double.eps) stop("'eps' should not be smaller than machine precision",call.=FALSE)

    force64 <- getOption("spam.force64")

    poselements <- (abs(x)>eps)
    if (any(!is.finite(x))) {
        poselements[!is.finite(x)] <- TRUE
    }
    lx <- length(x)
    nz <- sum(as.numeric(poselements))

    # empty matrix
    if (identical(nz,0)){
        return(
            .newSpamx(
                # rowpointers = c(1,rep_len64(2,lx)),
                dimension = c(lx,1),
                force64 = force64
            )
        )
    }

    return(
        .newSpamx(
            entries = as.double(x[poselements]),
            colindices = rep_len64(1, nz),
            rowpointers = cumsum(c(1, poselements)),
            dimension = c(lx,1),
            force64 = force64
        )
    )
}

#' @export
as.spamx.dist <- function(x, eps = getOption("spam.eps")) {
    if (eps<.Machine$double.eps) stop("'eps' should not be smaller than machine precision",call.=FALSE)
    if (any(!is.finite(x))) {
        warning("'NA/NaN/Inf' coerced to zero")
        x[!is.finite(x)] <- 0
    }
    dimx <- attr(x,"Size")
    size <- dimx*(dimx-1)/2

    force64 <- getOption("spam.force64")

    if(force64 || dimx^2 > 2147483647)
        SS <- .format64()
    else
        SS <- .format32

    z <- .C64("disttospam",
              SIGNATURE = c(SS$signature, "double",
                            "double", SS$signature, SS$signature,
                            "double"),

              nrow=dimx, #r
              x=as.vector(x,"double"), #r

              entries=vector_dc("double",size), #w
              colindices=vector_dc("integer",size), #w
              rowpointers=vector_dc("integer",dimx+1), #w

              eps=eps, #r

              INTENT = c("r", "r",
                         "w", "w", "w",
                         "r"),
              NAOK=getOption("spam.NAOK"),
              PACKAGE = SS$package
              )
    nz <- z$rowpointers[dimx+1]-1
    if(nz==0) return( .newSpamx(
                              # rowpointers = c( 1, rep( 2, dimx)),
                              dimension = c( dimx, dimx),
                              force64 = force64))


    return(.newSpamx(
        entries=z$entries[1:nz],
        colindices=z$colindices[1:nz],
        rowpointers=z$rowpointers[1:(dimx+1)],
        dimension=c(dimx,dimx),
        force64 = force64
    ))
}


#' @export
spamx.numeric <- function(x, nrow = 1, ncol = 1, eps = getOption("spam.eps"))  {
    force64 <- getOption("spam.force64")
    if (eps<.Machine$double.eps) stop("'eps' should not be smaller than machine precision",call.=FALSE)
    if (any(!is.finite(x))) {
        warning("'NA/NaN/Inf' coerced to zero")
        x[!is.finite(x)] <- 0
    }
    lenx <- length( x)
    if (missing(nrow))
        nrow <- ceiling( lenx/ncol)
    else if (missing(ncol))
        ncol <- ceiling( lenx/nrow)
    dimx <- c(nrow, ncol)
    if (lenx != prod(nrow,  ncol)) {
        if(lenx==1 && abs(x)<eps) {
            return(.newSpamx(
                dimension=dimx,
                force64=force64
            ))
        }
        else if(prod(nrow,ncol)%%lenx!=0)
            warning("ncol*nrow indivisable by length(x)")

        x <- rep_len64(x, prod(nrow,ncol))
    }


    dimx <- c(nrow, ncol)

    nz <- sum(as.numeric(abs(x) > eps))

    if(is.na(nz) | is.nan(nz))
        stop("NA or NaN in matrix.")

    if(nz==0) {
        return(.newSpamx(
            dimension=dimx,
            force64=force64
        ))
    }

    # EXPLICIT-STORAGE-FORMAT: Depending on the length of x, use 32 or 64-bit:
    if(force64 || nz > 2147483647 || max(dimx) > 2147483647 )
        SS <- .format64()
    else
        SS <- .format32

    z <- .C64("spamdnscsr",
              SIGNATURE=c(SS$signature, SS$signature, "double", SS$signature,
                          "double", SS$signature, SS$signature, "double"),

              nrow = dimx[1],
              ncol = dimx[2],
              x = x,
              dimx[1],

              entries = vector_dc("double",nz),
              colindices = vector_dc(SS$type,nz),
              rowpointers = vector_dc(SS$type,dimx[1]+1),
              eps = eps,

              INTENT=c("r", "r", "r", "r",
                       "w", "w", "w", "r"),
              NAOK = getOption("spam.NAOK"),
              PACKAGE = SS$package
              )

    entries <- z$entries
    colindices <- z$colindices
    rowpointers <- z$rowpointers
    z <- NULL

    return( .newSpamx( entries = entries,
                     colindices = colindices,
                     rowpointers = rowpointers,
                     dimension = dimx,
                     force64 = force64))
}

# setOldClass(c("dist", "numeric"))






# CAN ALSO use the syntax, but by default this is 64bit 
#obj <- new("employee",name="Steven", id=1002, contact="West Avenue")


#------------ from spam/spamlist.R
#' @export
as.spamx.list <- function(x, eps = getOption("spam.eps")) {
    spamx.list( x,  eps=eps)
}

#' @export
spamx.list <-  function(x, nrow, ncol, eps = getOption("spam.eps")) {
    force64 <- getOption("spam.force64")
    
    if (eps<.Machine$double.eps) stop("'eps' should not be smaller than machine precision",call.=FALSE)
    if (!is.list(x)|(length(x)<2)|(length(x)>3))
        stop("Argument 'x' needs to be a list with two or three elements")
                                        # two cases: list of length
                                        # -  two (matrix with two columns called ind* and the elements)
                                        # -  three (each one column called i*, j*.   

    if (identical(length(x),2L)) {
        indnr <- pmatch("ind",names(x)) 
        if (is.na(indnr)) stop("Argument 'x' needs an element called 'indices'")
        elenr <- ifelse( identical( indnr,1L), 2L, 1L)
        
        nz <- length( x[[elenr]])

        dimx <- dim(x[[indnr]])
        if (is.null(dimx)||(dimx[2] != 2))  stop("Indices should have two columns")
        if (dimx[1] != nz) stop("Number of indices does not match with number of elements")
        
        ir <- as.integer(x[[indnr]][,1])
        jc <- as.integer(x[[indnr]][,2])

        if(force64 || length(x[[elenr]]) > 2147483646)
            SS <- .format64()
        else
            SS <- .format32
        
    } else {
        inr <- pmatch("i",names(x)) 
        jnr <- pmatch("j",names(x))
        
        if (is.na(inr)||is.na(jnr)) stop("Argument 'x' needs elements called 'i' and 'j'")
        elenr <- c(1:3)[-c(inr,jnr)]
        nz <- length( x[[elenr]])
        
        ir <- as.integer(x[[inr]])
        jc <- as.integer(x[[jnr]])

        if ((length(ir) != nz)||(length(jc) != nz))
            stop("Number of indices does not match with number of elements")

        if(force64 || length(x[[elenr]]) > 2147483646)
            SS <- .format64()
        else
            SS <- .format32
    }
    
    if (nz == 0)
        return(.newSpamx(
            # rowpointers = c(1,rep_len64(2, nrow)),
            dimension = c(nrow,ncol),
            force64 = force64))
    if (any( ir <= 0) || any( jc <= 0))
        stop("Indices need to be positive")
    if (any(!is.finite(x[[elenr]]))) {
        warning("'NA/NaN/Inf' coerced to zero")
        x[[elenr]][!is.finite(x[[elenr]])] <- 0
    }
    nrow <- as.integer(ifelse(missing(nrow),max(ir),nrow))
    ncol <- as.integer(ifelse(missing(ncol),max(jc),ncol))
    ## z <- .Fortran(ifelse(toupper(getOption("spam.listmethod")=="PE"),"triplet3csr","triplet2csr"),
    ##               nrow=as.integer(nrow), ncol=as.integer(ncol),
    ##               nz=as.integer(nz),
    ##               as.double(x[[elenr]]),as.integer(ir),as.integer(jc),
    ##               entries=vector("double",nz),
    ##               colindices=vector("integer",nz),
    ##               rowpointers=vector("integer",nrow+1), as.double(eps),
    ##               NAOK=TRUE, PACKAGE = "spam"
    ##               )
    
    z <- .C64(ifelse(toupper(getOption("spam.listmethod")=="PE"),"triplet3csr","triplet2csr"),
              ## subroutine triplet3csr(nrow,ncol,nnz,a,ir,jc,ao,jao,iao,eps)
              ## subroutine triplet2csr(nrow,ncol,nnz,a,ir,jc,ao,jao,iao,eps)
              SIGNATURE = c( SS$signature, SS$signature, SS$signature,
                  "double", SS$signature, SS$signature,
                  "double", SS$signature, SS$signature,
                  "double"),
              
              nrow = nrow,
              ncol = ncol,
              nz = nz,
              
              x[[elenr]],
              ir,
              jc,
              
              entries = vector_dc( "double", nz),
              colindices = vector_dc( SS$type, nz),
              rowpointers = vector_dc( SS$type, nrow+1),
              
              eps,

              INTENT = c("r", "r", "rw",
                  "rw", "rw", "rw",
                  "rw", "rw", "rw",
                  "r"),
              NAOK=TRUE,
              PACKAGE = SS$package )

    

                                        #  print(z)
    if (z$nz == 0){
    ## if (identical(z$nz, 0)){
        ## print("special case")
        return(.newSpamx(
            # rowpointers = c(1, rep_len64(2,nrow)),
            dimension = c(nrow, ncol),
            force64 = force64))
         ## return(new("spam",rowpointers=c(1L,rep.int(2L,nrow)), dimension=c(nrow,ncol)))
    }
   
    ## newx <- new("spam")
    ## slot(newx,"entries",check=FALSE) <- z$entries[1:z$nz]
    ## slot(newx,"colindices",check=FALSE) <- z$colindices[1:z$nz]
    ## slot(newx,"rowpointers",check=FALSE) <- z$rowpointers
    ## slot(newx,"dimension",check=FALSE) <- c(nrow,ncol)
    ## return(newx)
    return(.newSpamx(
        entries = z$entries[1:z$nz],
        colindices = z$colindices[1:z$nz],
        rowpointers = z$rowpointers,
        dimension = c(nrow,ncol),
        force64 = force64))
}


#------------ END from spam/spamlist.R



#' @export
"as.spamx.list" <-  function(x, eps = getOption("spam.eps")) {
    spamx.list(x,eps)
}


setGeneric("as.spamx")
setMethod("as.spamx","spamx",  as.spamx.spamx)
setMethod("as.spamx","spam",   as.spamx.spam)
setMethod("as.spamx","matrix", as.spamx.matrix)
setMethod("as.spamx","numeric",as.spamx.numeric)
setMethod("as.spamx","dist",   as.spamx.dist)
setMethod("as.spamx","list",   as.spamx.list) #  { function(x,eps) spam.list(x,eps=eps)})

setGeneric("spamx")
setMethod("spamx","numeric",spamx.numeric)
setMethod("spamx","spam", as.spamx)
setMethod("spamx","spamx",as.spamx)
setMethod("spamx","list", spamx.list)
