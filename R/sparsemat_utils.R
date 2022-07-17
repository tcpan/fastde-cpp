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
    out <- rc_sp_transpose(x)
    # rownames(out) <- rownames(x)
    # colnames(out) <- colnames(x)
    return(out)
}

#' R Spamx transpose
#'
#' This implementation directly constructs the new transposed matrix.  
#'     There is random memory writes.
#' 
#' @rdname spamx_transpose
#' @param sp_matrix a sparse matrix, of the form dgCMatrix
#' @return matrix dense matrix.
#' @name spamx_transpose
#' @export
spamx_transpose <- function(x) {
    if (class(x@dimension) == "integer") {
        out <- rc_spamx32_transpose(x)
    } else {
        out <- rc_spamx64_transpose(x)
    }
    # rownames(out) <- rownames(x)
    # colnames(out) <- colnames(x)
    return(out)
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
    out <- rc_sp_to_dense(x)
    # rownames(out) <- rownames(x)
    # colnames(out) <- colnames(x)
    return(out)
}

#' R Sparse To Dense Matrix
#'
#' This implementation directly constructs the new dense matrix.  
#'     There is random memory writes.
#' 
#' @rdname spamx_to_dense
#' @param sp_matrix a sparse matrix, of the form dgCMatrix
#' @return matrix dense matrix.
#' @name spamx_to_dense
#' @export
spamx_to_dense <- function(x) {
    if (class(x@dimension) == "integer") {
        out <- rc_spamx32_to_dense(x)
    } else {
        out <- rc_spamx64_to_dense(x)
    }
    # rownames(out) <- rownames(x)
    # colnames(out) <- colnames(x)
    return(out)
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
    out <- rc_sp_to_dense_transposed(x)
    # rownames(out) <- colnames(x)
    # colnames(out) <- rownames(x)
    return(out)
}

#' R Sparse To Dense Matrix
#'
#' This implementation directly constructs the new dense matrix.  
#'     There is random memory writes.
#' 
#' @rdname spamx_to_dense_transposed
#' @param sp_matrix a sparse matrix, of the form dgCMatrix
#' @return matrix dense matrix.
#' @name spamx_to_dense_transposed
#' @export
spamx_to_dense_transposed <- function(x) {
     if (class(x@dimension) == "integer") {
        out <- rc_spamx32_to_dense_transposed(x)
    } else {
        out <- rc_spamx64_to_dense_transposed(x)
    }
    # rownames(out) <- colnames(x)
    # colnames(out) <- rownames(x)
    return(out)
}
