
#' Fast t-Test for dense matrix
#'
#' This implementation uses normal approximation, which works reasonably well if sample size is large (say N>=20)
#' 
#' @rdname ttest_fast
#' @param mat an expression matrix, COLUMN-MAJOR, each col is a feature, each row a sample
#' @param labels an integer vector, each element indicating the group to which a sample belongs.
#' @param alternative 
#' \itemize{
#' \item{0} : p(two.sided)
#' \item{1} : p(less)
#' \item{2} : p(greater)
#' }
#' @param var_equal TRUE/FALSE to indicate the variance is expected to be equal
#' @param as_dataframe TRUE/FALSE - TRUE returns a dataframe, FALSE returns a matrix
#' @param threads  number of concurrent threads.
#' @return array or dataframe.  for each gene/feature, the rows for the clusters are ordered by id.
#' @name ttest_fast
#' @export
ttest_fast <- function(mat, labels,
    alternative, var_equal, as_dataframe, threads) {
        out <- cpp11_dense_ttest(mat, colnames(mat), labels, 
        alternative, as.logical(var_equal), as.logical(as_dataframe), threads)
        
    if (!as_dataframe) {
        L <- unique(sort(labels))
        colnames(out) <- colnames(mat)
        rownames(out) <- L
    }

    return(out)
}

#' Fast T-Test for sparse matrix.  2 sample t-test.
#'
#' This implementation uses normal approximation, which works reasonably well if sample size is large (say N>=20)
#' 
#' @rdname sparse_ttest_fast
#' @param mat an expression matrix, COLUMN-MAJOR, each col is a feature, each row a sample
#' @param labels an integer vector, each element indicating the group to which a sample belongs.
#' @param alternative 
#' \itemize{
#' \item{0} : p(two.sided)
#' \item{1} : p(less)
#' \item{2} : p(greater)
#' }
#' @param var_equal TRUE/FALSE to indicate the variance is expected to be equal
#' @param as_dataframe TRUE/FALSE - TRUE returns a dataframe, FALSE returns a matrix
#' @param threads  number of concurrent threads.
#' @return array or dataframe.  for each gene/feature, the rows for the clusters are ordered by id.
#' @name sparse_ttest_fast
#' @export
sparse_ttest_fast <- function(mat, labels,
    features_as_rows, alternative, var_equal, as_dataframe, threads) {
    if (features_as_rows) 
        fnames <- rownames(mat)
    else 
        fnames <- colnames(mat)


    if (is(mat, 'dgCMatrix64')) {
        out <- cpp11_sparse64_ttest(mat@x, mat@i, mat@p, 
            fnames, nrow(mat), ncol(mat),
            labels, as.logical(features_as_rows), alternative, 
            as.logical(var_equal), as.logical(as_dataframe), threads)

    } else {
        out <- cpp11_sparse_ttest(mat@x, mat@i, mat@p, 
            fnames, nrow(mat), ncol(mat),
            labels, as.logical(features_as_rows), alternative, 
            as.logical(var_equal), as.logical(as_dataframe), threads)
    }
    if (!as_dataframe) {
        L <- unique(sort(labels))

        colnames(out) <- fnames
        rownames(out) <- L
    }

    return(out)

}
