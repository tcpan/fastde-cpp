library(tictoc)

#' R Sparse transpose
#'
#' This implementation directly constructs the new transposed matrix.  
#'     There is random memory writes.
#' 
#' @rdname sp_transpose
#' @param spmat a sparse matrix, of the form dgCMatrix
#' @return matrix dense matrix.
#' @name sp_transpose
#' @export
sp_transpose <- function(spmat, threads = 1) {
    if (is(spmat, 'dgCMatrix')) {
        tic("cpp transpose")
        mlist <- cpp11_sp_transpose(spmat@x, spmat@i, spmat@p, spmat@Dim[1], spmat@Dim[2], threads)
        toc()
        tic("create obj")
        out <- new("dgCMatrix", x=mlist$x, i=mlist$i, p=mlist$p, Dim=c(spmat@Dim[2], spmat@Dim[1]), Dimnames=list(colnames(spmat), rownames(spmat)))
        toc()
        return (out)
    } else if (is(spmat, 'dgCMatrix64')) {
        tic("cpp transpose")
        mlist <- cpp11_sp64_transpose(spmat@x, spmat@i, spmat@p, spmat@Dim[1], spmat@Dim[2], threads)
        toc()
        tic("create obj")
        out <- new("dgCMatrix64", x=mlist$x, i=mlist$i, p=mlist$p, Dim=c(spmat@Dim[2], spmat@Dim[1]), Dimnames=list(colnames(spmat), rownames(spmat)))
        toc()
        return(out)
    } else {
        print("USING R DEFAULT")
        return(t(spmat))
    }
    # # rownames(out) <- rownames(x)
    # # colnames(out) <- colnames(x)
    # return(out)
}


#' R Sparse To Dense Matrix
#'
#' This implementation directly constructs the new dense matrix.  
#'     There is random memory writes.
#' 
#' @rdname sp_to_dense
#' @param spmat a sparse matrix, of the form dgCMatrix
#' @return matrix dense matrix.
#' @name sp_to_dense
#' @export
sp_to_dense <- function(spmat, threads = 1) {
    if (is(spmat, 'dgCMatrix')) {
        m <- cpp11_sp_to_dense(spmat@x, spmat@i, spmat@p, spmat@Dim[1], spmat@Dim[2], threads)
        rownames(m) <- rownames(spmat)
        colnames(m) <- colnames(spmat)
        return(m)
    } else if (is(spmat, 'dgCMatrix64')) {
        m <- cpp11_sp64_to_dense(spmat@x, spmat@i, spmat@p, spmat@Dim[1], spmat@Dim[2], threads)
        rownames(m) <- rownames(spmat)
        colnames(m) <- colnames(spmat)
        return(m)
    } else {
        print("USING R DEFAULT")
        return(as.matrix(spmat))
    }
    # # rownames(out) <- rownames(x)
    # # colnames(out) <- colnames(x)
    # return(out)
}

#' R Sparse To Dense Matrix
#'
#' This implementation directly constructs the new dense matrix.  
#'     There is random memory writes.
#' 
#' @rdname sp_to_dense_transposed
#' @param spmat a sparse matrix, of the form dgCMatrix
#' @return matrix dense matrix.
#' @name sp_to_dense_transposed
#' @export
sp_to_dense_transposed <- function(spmat, threads = 1) {
    if (is(spmat, 'dgCMatrix')) {
        m <- cpp11_sp_to_dense_transposed(spmat@x, spmat@i, spmat@p, spmat@Dim[1], spmat@Dim[2], threads)
        rownames(m) <- colnames(spmat)
        colnames(m) <- rownames(spmat)
        return(m)
    } else if (is(spmat, 'dgCMatrix64')) {
        m <- cpp11_sp64_to_dense_transposed(spmat@x, spmat@i, spmat@p, spmat@Dim[1], spmat@Dim[2], threads)
        rownames(m) <- colnames(spmat)
        colnames(m) <- rownames(spmat)
        return(m)
    } else {
        print("USING R DEFAULT")
        return(as.matrix(t(spmat)))
    }
    # # rownames(out) <- colnames(x)
    # # colnames(out) <- rownames(x)
    # return(out)
}
