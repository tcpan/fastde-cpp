
#' Fast Wilcoxon-Mann-Whitney Test for dense matrix
#'
#' This implementation uses normal approximation, which works reasonably well if sample size is large (say N>=20)
#' 
#' @rdname wmw_fast
#' @param matrix an expression matrix, COLUMN-MAJOR, each col is a feature, each row a sample
#' @param labels an integer vector, each element indicating the group to which a sample belongs.
#' @param rtype 
#' \itemize{
#' \item{0} : p(less)
#' \item{1} : p(greater)
#' \item{2} : p(twoSided)
#' \item{3} : U
#' }
#' @param continuity_correction TRUE/FALSE for continuity correction
#' @param as_dataframe TRUE/FALSE - TRUE returns a dataframe, FALSE returns a matrix
#' @param threads  number of concurrent threads.
#' @return array or dataframe.  for each gene/feature, the rows for the clusters are ordered by id.
#' @name wmw_fast
#' @export
wmw_fast <- function(input, labels,
    rtype, continuity_correction, as_dataframe, threads) {

        out <- cpp11_dense_wmw(input, colnames(input), labels, 
        rtype, continuity_correction, as_dataframe, threads)
        
    if (!as_dataframe) {
        L <- unique(sort(labels))
        colnames(out) <- colnames(input)
        rownames(out) <- L
    }

    return(out)
}

#' Fast Wilcoxon-Mann-Whitney Test for dense matrix
#'
#' This implementation uses normal approximation, which works reasonably well if sample size is large (say N>=20)
#' 
#' @rdname sparse_wmw_fast
#' @param matrix an expression matrix, COLUMN-MAJOR, each col is a feature, each row a sample
#' @param labels an integer vector, each element indicating the group to which a sample belongs.
#' @param rtype 
#' \itemize{
#' \item{0} : p(less)
#' \item{1} : p(greater)
#' \item{2} : p(twoSided)
#' \item{3} : U
#' }
#' @param continuity_correction TRUE/FALSE for continuity_correction correction
#' @param as_dataframe TRUE/FALSE - TRUE returns a dataframe, FALSE returns a matrix
#' @param threads  number of concurrent threads.
#' @return array or dataframe.  for each gene/feature, the rows for the clusters are ordered by id.
#' @name sparse_wmw_fast
#' @export
sparse_wmw_fast <- function(matrix, labels,
    features_as_rows, rtype, continuity_correction, as_dataframe, threads) {
    if (features_as_rows) 
        fnames <- rownames(matrix)
    else 
        fnames <- colnames(matrix)


    if (is(matrix, 'dgCMatrix64')) {
        out <- cpp11_sparse64_wmw(matrix@x, matrix@i, matrix@p, 
            fnames, nrow(matrix), ncol(matrix),
            labels, features_as_rows, rtype, continuity_correction, as_dataframe, threads)

    } else {
        out <- cpp11_sparse_wmw(matrix@x, matrix@i, matrix@p, 
            fnames, nrow(matrix), ncol(matrix),
            labels, features_as_rows, rtype, continuity_correction, as_dataframe, threads)
    }
    if (!as_dataframe) {
        L <- unique(sort(labels))

        colnames(out) <- fnames
        rownames(out) <- L
    }

    return(out)

}
