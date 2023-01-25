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

#' R Sparse rbind
#'
#' This implementation allows production of very large sparse matrices.
#' 
#' @rdname sp_rbind
#' @param spmats a list of sparse matrice
#' @return a dgCMatrix64 object with the combined content.
#' @name sp_rbind
#' @export
sp_rbind <- function(spmats, threads = 1) {

    nmats = length(spmats)

    # create the data structures
    xs <- vector(mode = "list", length = nmats)
    iss <- vector(mode = "list", length = nmats)
    ps <- vector(mode = "list", length = nmats)
    nrs <- integer(nmats)
    ncs <- integer(nmats)
    nms <- character()

    tic("gathering data")
    for (i in 1:nmats) {
        xs[[i]] = spmats[[i]]@x
        iss[[i]] = spmats[[i]]@i
        ps[[i]] = spmats[[i]]@p
        nrs[[i]] = spmats[[i]]@Dim[1]
        ncs[[i]] = spmats[[i]]@Dim[2]

        nms <- c(nms, rownames(spmats[[i]]))
    }
    toc()

    tic("combine")
    # invoke the rbind
    if (is(spmats[[1]], 'dgCMatrix')) {
        m <- cpp11_sp_rbind(xs, iss, ps, nrs, ncs, threads)
    } else if (is(spmats[[1]], 'dgCMatrix64')) {
        m <- cpp11_sp64_rbind(xs, iss, ps, nrs, ncs, threads)
    } else {
        print("UNSUPPORTED.  Matrices must be sparse matrices")
        return(NULL)
    }
    toc()

    tic("build output")
    if (length(m[[3]]) <= .Machine$integer.max) {
        ip <- integer(length=length(m[[3]]))
        for (ii in 1:length(m[[3]])) {
            ip[[ii]] <- m[[3]][[ii]]
        }
        out <- new('dgCMatrix', x = m[[1]], i = m[[2]], p = ip, 
            Dim = c(m[[4]], m[[5]]))
    } else {
        out <- new('dgCMatrix64', x = m[[1]], i = m[[2]], p = m[[3]], 
            Dim = c(m[[4]], m[[5]]))
    }
    rownames(out) <- nms
    colnames(out) <- colnames(spmats[[1]])
    toc()

    return(out)

}


#' R Sparse cbind
#'
#' This implementation allows production of very large sparse matrices.
#' 
#' @rdname sp_cbind
#' @param spmats a list of sparse matrice
#' @return a dgCMatrix64 object with the combined content.
#' @name sp_cbind
#' @export
sp_cbind <- function(spmats, threads = 1) {

    nmats = length(spmats)

    # create the data structures
    xs <- vector(mode = "list", length = nmats)
    iss <- vector(mode = "list", length = nmats)
    ps <- vector(mode = "list", length = nmats)
    nrs <- integer(nmats)
    ncs <- integer(nmats)
    nms <- character()

    tic("gather data")
    for (i in 1:nmats) {
        xs[[i]] = spmats[[i]]@x
        iss[[i]] = spmats[[i]]@i
        ps[[i]] = spmats[[i]]@p
        nrs[[i]] = spmats[[i]]@Dim[1]
        ncs[[i]] = spmats[[i]]@Dim[2]

        nms <- c(nms, colnames(spmats[[i]]))
    }
    toc()

    tic("combine")
    # invoke the rbind
    if (is(spmats[[1]], 'dgCMatrix')) {
        m <- cpp11_sp_cbind(xs, iss, ps, nrs, ncs, threads)
    } else if (is(spmats[[1]], 'dgCMatrix64')) {
        m <- cpp11_sp64_cbind(xs, iss, ps, nrs, ncs, threads)
    } else {
        print("UNSUPPORTED.  Matrices must be sparse matrices")
        return(NULL)
    }
    toc()

    tic("create output")
    if (length(m[[3]]) <= .Machine$integer.max) {
        ip <- integer(length=length(m[[3]]))
        for (ii in 1:length(m[[3]])) {
            ip[[ii]] <- m[[3]][[ii]]
        }
        out <- new('dgCMatrix', x = m[[1]], i = m[[2]], p = ip, 
            Dim = c(m[[4]], m[[5]]))
    } else {
        out <- new('dgCMatrix64', x = m[[1]], i = m[[2]], p = m[[3]], 
            Dim = c(m[[4]], m[[5]]))
    }
    rownames(out) <- rownames(spmats[[1]])
    colnames(out) <- nms
    toc()

    return(out)

}


#' R Sparse rowSums
#'
#' This implementation allows production of very large sparse matrices.
#' 
#' @rdname sp_rowSums
#' @param spmats a list of sparse matrice
#' @return a dense array of row sums.
#' @name sp_rowSums
#' @export
sp_rowSums <- function(spmat, threads = 1) {
    tic("[TIME] Row sums")
    if (is(spmat, 'dgCMatrix')) {
        if (threads == 1) {
            m <- rowSums(spmat)
        } else {
            m <- cpp11_sp_rowSums(spmat@x, spmat@i, spmat@Dim[[2]], threads)
        }
    } else if (is(spmat, 'dgCMatrix64')) {
        m <- cpp11_sp_rowSums(spmat@x, spmat@i, spmat@Dim[[2]], threads)
    } else {
        print("UNSUPPORTED.  Matrices must be sparse matrices")
        return(NULL)
    }
    names(m) <- rownames(spmat)
    toc()
    return(m)
}

#' R Sparse rowSums
#'
#' This implementation allows production of very large sparse matrices.
#' 
#' @rdname sp_rowSums
#' @param spmats a list of sparse matrice
#' @return a dense array of row sums.
#' @name sp_rowSums
#' @export
sp_colSums <- function(spmat, threads = 1) {
    tic("[TIME] Row sums")
    if (is(spmat, 'dgCMatrix')) {
        if (threads == 1) {
            m <- colSums(spmat)
        } else {
            m <- cpp11_sp_colSums(spmat@x, spmat@p, threads)
        }
    } else if (is(spmat, 'dgCMatrix64')) {
        m <- cpp11_sp64_colSums(spmat@x, spmat@p, threads)
    } else {
        print("UNSUPPORTED.  Matrices must be sparse matrices")
        return(NULL)
    }
    names(m) <- colnames(spmat)

    toc()
    return(m)
}