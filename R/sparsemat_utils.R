#' R Sparse transpose
#'
#' This implementation directly constructs the new transposed matrix.  
#'     There is random memory writes.
#' 
#' @rdname sp_transpose
#' @param sp_matrix a sparse matrix, of the form dgCMatrix
#' @return matrix dense matrix.
#' @name sp_transpose
#' @export
sp_transpose <- function(x) {
    if (class(x) == "dgCMatrix")
        return(rc_sp_transpose(x))
    else if (class(x) == "dgCMatrix64") {
        return(rc_sp64_transpose(x))
    } else {
        return(t(x))
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
#' @param sp_matrix a sparse matrix, of the form dgCMatrix
#' @return matrix dense matrix.
#' @name sp_to_dense
#' @export
sp_to_dense <- function(x) {
    if (class(x) == "dgCMatrix")
        return(rc_sp_to_dense(x))
    else if (class(x) == "dgCMatrix64") {
        return(rc_sp64_to_dense(x))
    } else {
        return(as.matrix(x))
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
#' @param sp_matrix a sparse matrix, of the form dgCMatrix
#' @return matrix dense matrix.
#' @name sp_to_dense_transposed
#' @export
sp_to_dense_transposed <- function(x) {
    if (class(x) == "dgCMatrix")
        return(rc_sp_to_dense_transposed(x))
    else if (class(x) == "dgCMatrix64") {
        return(rc_sp64_to_dense_transposed(x))
    } else {
        return(as.matrix(t(x)))
    }
    # # rownames(out) <- colnames(x)
    # # colnames(out) <- rownames(x)
    # return(out)
}
