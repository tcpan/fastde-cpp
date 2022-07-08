// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "fastde_types.hpp"
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rttest_dgCMatrix
Rcpp::dgCMatrix rttest_dgCMatrix(Rcpp::dgCMatrix& mat);
RcppExport SEXP _fastde_rttest_dgCMatrix(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::dgCMatrix& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(rttest_dgCMatrix(mat));
    return rcpp_result_gen;
END_RCPP
}
// rttest_spamx32
Rcpp::spamx32 rttest_spamx32(Rcpp::spamx32& mat);
RcppExport SEXP _fastde_rttest_spamx32(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::spamx32& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(rttest_spamx32(mat));
    return rcpp_result_gen;
END_RCPP
}
// rttest_spamx64
Rcpp::spamx64 rttest_spamx64(Rcpp::spamx64& mat);
RcppExport SEXP _fastde_rttest_spamx64(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::spamx64& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(rttest_spamx64(mat));
    return rcpp_result_gen;
END_RCPP
}
// rttest_spam32
Rcpp::spam32 rttest_spam32(Rcpp::spam32& mat);
RcppExport SEXP _fastde_rttest_spam32(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::spam32& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(rttest_spam32(mat));
    return rcpp_result_gen;
END_RCPP
}
// rttest_spam64
Rcpp::spam64 rttest_spam64(Rcpp::spam64& mat);
RcppExport SEXP _fastde_rttest_spam64(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::spam64& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(rttest_spam64(mat));
    return rcpp_result_gen;
END_RCPP
}
// ComputeFoldChange
extern SEXP ComputeFoldChange(Rcpp::NumericMatrix matrix, Rcpp::IntegerVector labels, SEXP calc_percents, SEXP fc_name, SEXP use_expm1, SEXP min_threshold, SEXP use_log, SEXP log_base, SEXP use_pseudocount, SEXP as_dataframe, SEXP threads);
RcppExport SEXP _fastde_ComputeFoldChange(SEXP matrixSEXP, SEXP labelsSEXP, SEXP calc_percentsSEXP, SEXP fc_nameSEXP, SEXP use_expm1SEXP, SEXP min_thresholdSEXP, SEXP use_logSEXP, SEXP log_baseSEXP, SEXP use_pseudocountSEXP, SEXP as_dataframeSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type matrix(matrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type calc_percents(calc_percentsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type fc_name(fc_nameSEXP);
    Rcpp::traits::input_parameter< SEXP >::type use_expm1(use_expm1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type min_threshold(min_thresholdSEXP);
    Rcpp::traits::input_parameter< SEXP >::type use_log(use_logSEXP);
    Rcpp::traits::input_parameter< SEXP >::type log_base(log_baseSEXP);
    Rcpp::traits::input_parameter< SEXP >::type use_pseudocount(use_pseudocountSEXP);
    Rcpp::traits::input_parameter< SEXP >::type as_dataframe(as_dataframeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(ComputeFoldChange(matrix, labels, calc_percents, fc_name, use_expm1, min_threshold, use_log, log_base, use_pseudocount, as_dataframe, threads));
    return rcpp_result_gen;
END_RCPP
}
// ComputeSparseFoldChange
extern SEXP ComputeSparseFoldChange(Rcpp::dgCMatrix matrix, Rcpp::IntegerVector labels, SEXP calc_percents, SEXP fc_name, SEXP use_expm1, SEXP min_threshold, SEXP use_log, SEXP log_base, SEXP use_pseudocount, SEXP as_dataframe, SEXP threads);
RcppExport SEXP _fastde_ComputeSparseFoldChange(SEXP matrixSEXP, SEXP labelsSEXP, SEXP calc_percentsSEXP, SEXP fc_nameSEXP, SEXP use_expm1SEXP, SEXP min_thresholdSEXP, SEXP use_logSEXP, SEXP log_baseSEXP, SEXP use_pseudocountSEXP, SEXP as_dataframeSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::dgCMatrix >::type matrix(matrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type calc_percents(calc_percentsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type fc_name(fc_nameSEXP);
    Rcpp::traits::input_parameter< SEXP >::type use_expm1(use_expm1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type min_threshold(min_thresholdSEXP);
    Rcpp::traits::input_parameter< SEXP >::type use_log(use_logSEXP);
    Rcpp::traits::input_parameter< SEXP >::type log_base(log_baseSEXP);
    Rcpp::traits::input_parameter< SEXP >::type use_pseudocount(use_pseudocountSEXP);
    Rcpp::traits::input_parameter< SEXP >::type as_dataframe(as_dataframeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(ComputeSparseFoldChange(matrix, labels, calc_percents, fc_name, use_expm1, min_threshold, use_log, log_base, use_pseudocount, as_dataframe, threads));
    return rcpp_result_gen;
END_RCPP
}
// FilterFoldChange
extern SEXP FilterFoldChange(SEXP fc, SEXP pct1, SEXP pct2, SEXP init_mask, SEXP min_pct, SEXP min_diff_pct, SEXP logfc_threshold, SEXP only_pos, SEXP not_count, SEXP threads);
RcppExport SEXP _fastde_FilterFoldChange(SEXP fcSEXP, SEXP pct1SEXP, SEXP pct2SEXP, SEXP init_maskSEXP, SEXP min_pctSEXP, SEXP min_diff_pctSEXP, SEXP logfc_thresholdSEXP, SEXP only_posSEXP, SEXP not_countSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type fc(fcSEXP);
    Rcpp::traits::input_parameter< SEXP >::type pct1(pct1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type pct2(pct2SEXP);
    Rcpp::traits::input_parameter< SEXP >::type init_mask(init_maskSEXP);
    Rcpp::traits::input_parameter< SEXP >::type min_pct(min_pctSEXP);
    Rcpp::traits::input_parameter< SEXP >::type min_diff_pct(min_diff_pctSEXP);
    Rcpp::traits::input_parameter< SEXP >::type logfc_threshold(logfc_thresholdSEXP);
    Rcpp::traits::input_parameter< SEXP >::type only_pos(only_posSEXP);
    Rcpp::traits::input_parameter< SEXP >::type not_count(not_countSEXP);
    Rcpp::traits::input_parameter< SEXP >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(FilterFoldChange(fc, pct1, pct2, init_mask, min_pct, min_diff_pct, logfc_threshold, only_pos, not_count, threads));
    return rcpp_result_gen;
END_RCPP
}
// sp_transpose
extern SEXP sp_transpose(SEXP sp_matrix);
RcppExport SEXP _fastde_sp_transpose(SEXP sp_matrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type sp_matrix(sp_matrixSEXP);
    rcpp_result_gen = Rcpp::wrap(sp_transpose(sp_matrix));
    return rcpp_result_gen;
END_RCPP
}
// sp_to_dense
extern SEXP sp_to_dense(SEXP sp_matrix);
RcppExport SEXP _fastde_sp_to_dense(SEXP sp_matrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type sp_matrix(sp_matrixSEXP);
    rcpp_result_gen = Rcpp::wrap(sp_to_dense(sp_matrix));
    return rcpp_result_gen;
END_RCPP
}
// ttest_fast
extern SEXP ttest_fast(SEXP matrix, SEXP labels, SEXP alternative, SEXP var_equal, SEXP as_dataframe, SEXP threads);
RcppExport SEXP _fastde_ttest_fast(SEXP matrixSEXP, SEXP labelsSEXP, SEXP alternativeSEXP, SEXP var_equalSEXP, SEXP as_dataframeSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type matrix(matrixSEXP);
    Rcpp::traits::input_parameter< SEXP >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type alternative(alternativeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type var_equal(var_equalSEXP);
    Rcpp::traits::input_parameter< SEXP >::type as_dataframe(as_dataframeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(ttest_fast(matrix, labels, alternative, var_equal, as_dataframe, threads));
    return rcpp_result_gen;
END_RCPP
}
// sparse_ttest_fast
extern SEXP sparse_ttest_fast(SEXP matrix, SEXP labels, SEXP alternative, SEXP var_equal, SEXP as_dataframe, SEXP threads);
RcppExport SEXP _fastde_sparse_ttest_fast(SEXP matrixSEXP, SEXP labelsSEXP, SEXP alternativeSEXP, SEXP var_equalSEXP, SEXP as_dataframeSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type matrix(matrixSEXP);
    Rcpp::traits::input_parameter< SEXP >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type alternative(alternativeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type var_equal(var_equalSEXP);
    Rcpp::traits::input_parameter< SEXP >::type as_dataframe(as_dataframeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(sparse_ttest_fast(matrix, labels, alternative, var_equal, as_dataframe, threads));
    return rcpp_result_gen;
END_RCPP
}
// wmw_fast
extern SEXP wmw_fast(SEXP matrix, SEXP labels, SEXP rtype, SEXP continuity_correction, SEXP as_dataframe, SEXP threads);
RcppExport SEXP _fastde_wmw_fast(SEXP matrixSEXP, SEXP labelsSEXP, SEXP rtypeSEXP, SEXP continuity_correctionSEXP, SEXP as_dataframeSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type matrix(matrixSEXP);
    Rcpp::traits::input_parameter< SEXP >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type rtype(rtypeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type continuity_correction(continuity_correctionSEXP);
    Rcpp::traits::input_parameter< SEXP >::type as_dataframe(as_dataframeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(wmw_fast(matrix, labels, rtype, continuity_correction, as_dataframe, threads));
    return rcpp_result_gen;
END_RCPP
}
// sparse_wmw_fast
extern SEXP sparse_wmw_fast(SEXP matrix, SEXP labels, SEXP rtype, SEXP continuity_correction, SEXP as_dataframe, SEXP threads);
RcppExport SEXP _fastde_sparse_wmw_fast(SEXP matrixSEXP, SEXP labelsSEXP, SEXP rtypeSEXP, SEXP continuity_correctionSEXP, SEXP as_dataframeSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type matrix(matrixSEXP);
    Rcpp::traits::input_parameter< SEXP >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type rtype(rtypeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type continuity_correction(continuity_correctionSEXP);
    Rcpp::traits::input_parameter< SEXP >::type as_dataframe(as_dataframeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(sparse_wmw_fast(matrix, labels, rtype, continuity_correction, as_dataframe, threads));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fastde_rttest_dgCMatrix", (DL_FUNC) &_fastde_rttest_dgCMatrix, 1},
    {"_fastde_rttest_spamx32", (DL_FUNC) &_fastde_rttest_spamx32, 1},
    {"_fastde_rttest_spamx64", (DL_FUNC) &_fastde_rttest_spamx64, 1},
    {"_fastde_rttest_spam32", (DL_FUNC) &_fastde_rttest_spam32, 1},
    {"_fastde_rttest_spam64", (DL_FUNC) &_fastde_rttest_spam64, 1},
    {"_fastde_ComputeFoldChange", (DL_FUNC) &_fastde_ComputeFoldChange, 11},
    {"_fastde_ComputeSparseFoldChange", (DL_FUNC) &_fastde_ComputeSparseFoldChange, 11},
    {"_fastde_FilterFoldChange", (DL_FUNC) &_fastde_FilterFoldChange, 10},
    {"_fastde_sp_transpose", (DL_FUNC) &_fastde_sp_transpose, 1},
    {"_fastde_sp_to_dense", (DL_FUNC) &_fastde_sp_to_dense, 1},
    {"_fastde_ttest_fast", (DL_FUNC) &_fastde_ttest_fast, 6},
    {"_fastde_sparse_ttest_fast", (DL_FUNC) &_fastde_sparse_ttest_fast, 6},
    {"_fastde_wmw_fast", (DL_FUNC) &_fastde_wmw_fast, 6},
    {"_fastde_sparse_wmw_fast", (DL_FUNC) &_fastde_sparse_wmw_fast, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_fastde(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
