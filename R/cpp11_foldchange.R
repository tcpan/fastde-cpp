#' Fold Change
#' 
#' https://stackoverflow.com/questions/38338270/how-to-return-a-named-vecsxp-when-writing-r-extensions
#' 
#' @rdname ComputeFoldChange
#' @param mat an expression matrix, COLUMN-MAJOR, each row is a sample, each column a gene
#' @param labels an integer vector, each element indicating the group to which a sample belongs.
#' @param calc_percents  a boolean to indicate whether to compute percents or not.
#' @param fc_name column name to use for the fold change results 
#' @param use_expm1 for "data", use expm1
#' @param min_threshold minimum threshold to count towards pct.1 and pct.2 percentages.
#' @param use_log for "data" and default log type, indicate log of the sum is to be used.
#' @param log_base base for the log
#' @param use_pseudocount for "data" and default log type, add pseudocount after log.
#' @param as_dataframe TRUE/FALSE.  TRUE = return a linearized dataframe.  FALSE = return matrices.
#' @param threads number of threads to use
#' @return array or dataframe
#' @name ComputeFoldChange
#' @export
ComputeFoldChange <- function(mat, labels, 
    calc_percents, fc_name, use_expm1, min_threshold, 
    use_log, log_base, use_pseudocount, as_dataframe, threads) {

    out <- cpp11_ComputeFoldChange(mat, colnames(mat), 
        labels, as.logical(calc_percents), fc_name, as.logical(use_expm1),
         min_threshold, 
        as.logical(use_log), log_base, as.logical(use_pseudocount), 
        as.logical(as_dataframe), threads)


    if (!as_dataframe) {
        L <- unique(sort(labels))
        if (calc_percents) {
            colnames(out[[fc_name]]) <- colnames(mat)
            rownames(out[[fc_name]]) <- L
        } else {
            colnames(out) <- colnames(mat)
            rownames(out) <- L
        }
    }

    return(out)
}

#' Fold Change
#' 
#' https://stackoverflow.com/questions/38338270/how-to-return-a-named-vecsxp-when-writing-r-extensions
#' 
#' @rdname ComputeFoldChangeSparse
#' @param mat an expression matrix, COLUMN-MAJOR, each row is a sample, each column a gene
#' @param labels an integer vector, each element indicating the group to which a sample belongs.
#' @param calc_percents  a boolean to indicate whether to compute percents or not.
#' @param fc_name column name to use for the fold change results 
#' @param use_expm1 for "data", use expm1
#' @param min_threshold minimum threshold to count towards pct.1 and pct.2 percentages.
#' @param use_log for "data" and default log type, indicate log of the sum is to be used.
#' @param log_base base for the log
#' @param use_pseudocount for "data" and default log type, add pseudocount after log.
#' @param as_dataframe TRUE/FALSE.  TRUE = return a linearized dataframe.  FALSE = return matrices.
#' @param threads number of threads to use
#' @return array or dataframe
#' @name ComputeFoldChangeSparse
#' @export
ComputeFoldChangeSparse <- function(mat, labels, 
    features_as_rows,
    calc_percents, fc_name, use_expm1, min_threshold, 
    use_log, log_base, use_pseudocount, as_dataframe, threads) {

    if (features_as_rows) 
        fnames <- rownames(mat)
    else 
        fnames <- colnames(mat)

    compute <- if (is(mat, 'dgCMatrix64')) {
        cpp11_ComputeFoldChangeSparse64
    } else {
        cpp11_ComputeFoldChangeSparse
    }
    out <- compute(x=mat@x, i=mat@i, p=mat@p, 
            features = fnames, rows = nrow(mat), cols = ncol(mat),
            labels = labels, features_as_rows = as.logical(features_as_rows),
            calc_percents = as.logical(calc_percents), 
            fc_name= fc_name, use_expm1=as.logical(use_expm1), 
            min_threshold=min_threshold, 
            use_log=as.logical(use_log), log_base=log_base, 
            use_pseudocount=as.logical(use_pseudocount), 
            as_dataframe=as.logical(as_dataframe), threads= threads)

    if (!as_dataframe) {
        L <- unique(sort(labels))
        if (calc_percents) {
            colnames(out[[fc_name]]) <- fnames
            rownames(out[[fc_name]]) <- L
        } else {
            colnames(out) <- fnames
            rownames(out) <- L
        }
    }

    return(out)

}



#' Filter based on FoldChange
#' 
#'  https://stackoverflow.com/questions/38338270/how-to-return-a-named-vecsxp-when-writing-r-extensions
#' 
#' @rdname FilterFoldChange
#' @param fc foldchange values, either as a vector or a matrix.  if matrix: 1 gene per col, 1 clust per row
#' @param pct1 percent greater than threshold (0) in class 1.
#' @param pct2 percent greater than threshold (0) in class 2.
#' @param init_mask  initial mask, based on an external "features" vector.
#' @param min_pct minimum threshold for max pct1, pct2
#' @param min_diff_pct minimum threshold for difference between max and min {pct1, pct2}
#' @param logfc_threshold if not scaled.data, then compare to the logfc.
#' @param only_pos keep only positive fc value, and not use abs when thresholding.
#' @param not_count  not scaled.data
#' @param threads number of threads to use
#' @return array of same shape as fc (dense)
#' @name FilterFoldChange
#' @export
FilterFoldChange <- function(fc, pct1, pct2, init_mask, 
    min_pct, min_diff_pct, logfc_threshold, only_pos, not_count, threads) {
    compute <- if (is.matrix(fc)) {
        cpp11_FilterFoldChangeMat
    } else {
        cpp11_FilterFoldChange
    }
    if (is.null(init_mask)) {
        imask <- rep(TRUE, length(pct1))
    } else {
        imask <- init_mask
    }
    return(compute(fc = fc, pct1 = pct1, pct2 = pct2, init_mask = imask, 
            min_pct = min_pct, min_diff_pct = min_diff_pct, 
            logfc_threshold = logfc_threshold, only_pos = as.logical(only_pos),
            not_count = as.logical(not_count), threads = threads))
}
