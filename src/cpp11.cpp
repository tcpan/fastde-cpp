// Generated by cpp11: do not edit by hand
// clang-format off


#include "cpp11/declarations.hpp"
#include <R_ext/Visibility.h>

// cpp11_foldchange.cpp
extern cpp11::sexp cpp11_ComputeFoldChange(cpp11::doubles_matrix<cpp11::by_column> const & matrix, cpp11::strings const & features, cpp11::integers const & labels, bool calc_percents, std::string fc_name, bool use_expm1, double min_threshold, bool use_log, double log_base, bool use_pseudocount, bool as_dataframe, int threads);
extern "C" SEXP _fastde_cpp11_ComputeFoldChange(SEXP matrix, SEXP features, SEXP labels, SEXP calc_percents, SEXP fc_name, SEXP use_expm1, SEXP min_threshold, SEXP use_log, SEXP log_base, SEXP use_pseudocount, SEXP as_dataframe, SEXP threads) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp11_ComputeFoldChange(cpp11::as_cpp<cpp11::decay_t<cpp11::doubles_matrix<cpp11::by_column> const &>>(matrix), cpp11::as_cpp<cpp11::decay_t<cpp11::strings const &>>(features), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(labels), cpp11::as_cpp<cpp11::decay_t<bool>>(calc_percents), cpp11::as_cpp<cpp11::decay_t<std::string>>(fc_name), cpp11::as_cpp<cpp11::decay_t<bool>>(use_expm1), cpp11::as_cpp<cpp11::decay_t<double>>(min_threshold), cpp11::as_cpp<cpp11::decay_t<bool>>(use_log), cpp11::as_cpp<cpp11::decay_t<double>>(log_base), cpp11::as_cpp<cpp11::decay_t<bool>>(use_pseudocount), cpp11::as_cpp<cpp11::decay_t<bool>>(as_dataframe), cpp11::as_cpp<cpp11::decay_t<int>>(threads)));
  END_CPP11
}
// cpp11_foldchange.cpp
extern cpp11::sexp cpp11_ComputeFoldChangeSparse(cpp11::doubles const & x, cpp11::integers const & i, cpp11::integers const & p, cpp11::strings const & features, int const & rows, int const & cols, cpp11::integers const & labels, bool features_as_rows, bool calc_percents, std::string fc_name, bool use_expm1, double min_threshold, bool use_log, double log_base, bool use_pseudocount, bool as_dataframe, int threads);
extern "C" SEXP _fastde_cpp11_ComputeFoldChangeSparse(SEXP x, SEXP i, SEXP p, SEXP features, SEXP rows, SEXP cols, SEXP labels, SEXP features_as_rows, SEXP calc_percents, SEXP fc_name, SEXP use_expm1, SEXP min_threshold, SEXP use_log, SEXP log_base, SEXP use_pseudocount, SEXP as_dataframe, SEXP threads) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp11_ComputeFoldChangeSparse(cpp11::as_cpp<cpp11::decay_t<cpp11::doubles const &>>(x), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(i), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(p), cpp11::as_cpp<cpp11::decay_t<cpp11::strings const &>>(features), cpp11::as_cpp<cpp11::decay_t<int const &>>(rows), cpp11::as_cpp<cpp11::decay_t<int const &>>(cols), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(labels), cpp11::as_cpp<cpp11::decay_t<bool>>(features_as_rows), cpp11::as_cpp<cpp11::decay_t<bool>>(calc_percents), cpp11::as_cpp<cpp11::decay_t<std::string>>(fc_name), cpp11::as_cpp<cpp11::decay_t<bool>>(use_expm1), cpp11::as_cpp<cpp11::decay_t<double>>(min_threshold), cpp11::as_cpp<cpp11::decay_t<bool>>(use_log), cpp11::as_cpp<cpp11::decay_t<double>>(log_base), cpp11::as_cpp<cpp11::decay_t<bool>>(use_pseudocount), cpp11::as_cpp<cpp11::decay_t<bool>>(as_dataframe), cpp11::as_cpp<cpp11::decay_t<int>>(threads)));
  END_CPP11
}
// cpp11_foldchange.cpp
extern cpp11::sexp cpp11_ComputeFoldChangeSparse64(cpp11::doubles const & x, cpp11::integers const & i, cpp11::doubles const & p, cpp11::strings const & features, int const & rows, int const & cols, cpp11::integers const & labels, bool features_as_rows, bool calc_percents, std::string fc_name, bool use_expm1, double min_threshold, bool use_log, double log_base, bool use_pseudocount, bool as_dataframe, int threads);
extern "C" SEXP _fastde_cpp11_ComputeFoldChangeSparse64(SEXP x, SEXP i, SEXP p, SEXP features, SEXP rows, SEXP cols, SEXP labels, SEXP features_as_rows, SEXP calc_percents, SEXP fc_name, SEXP use_expm1, SEXP min_threshold, SEXP use_log, SEXP log_base, SEXP use_pseudocount, SEXP as_dataframe, SEXP threads) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp11_ComputeFoldChangeSparse64(cpp11::as_cpp<cpp11::decay_t<cpp11::doubles const &>>(x), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(i), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles const &>>(p), cpp11::as_cpp<cpp11::decay_t<cpp11::strings const &>>(features), cpp11::as_cpp<cpp11::decay_t<int const &>>(rows), cpp11::as_cpp<cpp11::decay_t<int const &>>(cols), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(labels), cpp11::as_cpp<cpp11::decay_t<bool>>(features_as_rows), cpp11::as_cpp<cpp11::decay_t<bool>>(calc_percents), cpp11::as_cpp<cpp11::decay_t<std::string>>(fc_name), cpp11::as_cpp<cpp11::decay_t<bool>>(use_expm1), cpp11::as_cpp<cpp11::decay_t<double>>(min_threshold), cpp11::as_cpp<cpp11::decay_t<bool>>(use_log), cpp11::as_cpp<cpp11::decay_t<double>>(log_base), cpp11::as_cpp<cpp11::decay_t<bool>>(use_pseudocount), cpp11::as_cpp<cpp11::decay_t<bool>>(as_dataframe), cpp11::as_cpp<cpp11::decay_t<int>>(threads)));
  END_CPP11
}
// cpp11_foldchange.cpp
extern cpp11::writable::logicals cpp11_FilterFoldChange(cpp11::doubles const & fc, cpp11::doubles const & pct1, cpp11::doubles const & pct2, cpp11::logicals const & init_mask, double min_pct, double min_diff_pct, double logfc_threshold, bool only_pos, bool not_count, int threads);
extern "C" SEXP _fastde_cpp11_FilterFoldChange(SEXP fc, SEXP pct1, SEXP pct2, SEXP init_mask, SEXP min_pct, SEXP min_diff_pct, SEXP logfc_threshold, SEXP only_pos, SEXP not_count, SEXP threads) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp11_FilterFoldChange(cpp11::as_cpp<cpp11::decay_t<cpp11::doubles const &>>(fc), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles const &>>(pct1), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles const &>>(pct2), cpp11::as_cpp<cpp11::decay_t<cpp11::logicals const &>>(init_mask), cpp11::as_cpp<cpp11::decay_t<double>>(min_pct), cpp11::as_cpp<cpp11::decay_t<double>>(min_diff_pct), cpp11::as_cpp<cpp11::decay_t<double>>(logfc_threshold), cpp11::as_cpp<cpp11::decay_t<bool>>(only_pos), cpp11::as_cpp<cpp11::decay_t<bool>>(not_count), cpp11::as_cpp<cpp11::decay_t<int>>(threads)));
  END_CPP11
}
// cpp11_foldchange.cpp
extern cpp11::writable::logicals_matrix<cpp11::by_column> cpp11_FilterFoldChangeMat(cpp11::doubles_matrix<cpp11::by_column> const & fc, cpp11::doubles_matrix<cpp11::by_column> const & pct1, cpp11::doubles_matrix<cpp11::by_column> const & pct2, cpp11::logicals_matrix<cpp11::by_column> const & init_mask, double min_pct, double min_diff_pct, double logfc_threshold, bool only_pos, bool not_count, int threads);
extern "C" SEXP _fastde_cpp11_FilterFoldChangeMat(SEXP fc, SEXP pct1, SEXP pct2, SEXP init_mask, SEXP min_pct, SEXP min_diff_pct, SEXP logfc_threshold, SEXP only_pos, SEXP not_count, SEXP threads) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp11_FilterFoldChangeMat(cpp11::as_cpp<cpp11::decay_t<cpp11::doubles_matrix<cpp11::by_column> const &>>(fc), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles_matrix<cpp11::by_column> const &>>(pct1), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles_matrix<cpp11::by_column> const &>>(pct2), cpp11::as_cpp<cpp11::decay_t<cpp11::logicals_matrix<cpp11::by_column> const &>>(init_mask), cpp11::as_cpp<cpp11::decay_t<double>>(min_pct), cpp11::as_cpp<cpp11::decay_t<double>>(min_diff_pct), cpp11::as_cpp<cpp11::decay_t<double>>(logfc_threshold), cpp11::as_cpp<cpp11::decay_t<bool>>(only_pos), cpp11::as_cpp<cpp11::decay_t<bool>>(not_count), cpp11::as_cpp<cpp11::decay_t<int>>(threads)));
  END_CPP11
}
// cpp11_sparsemat.cpp
extern cpp11::writable::list cpp11_sp_transpose(cpp11::doubles const & x, cpp11::integers const & i, cpp11::integers const & p, int const & nrow, int const & ncol, int const & threads);
extern "C" SEXP _fastde_cpp11_sp_transpose(SEXP x, SEXP i, SEXP p, SEXP nrow, SEXP ncol, SEXP threads) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp11_sp_transpose(cpp11::as_cpp<cpp11::decay_t<cpp11::doubles const &>>(x), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(i), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(p), cpp11::as_cpp<cpp11::decay_t<int const &>>(nrow), cpp11::as_cpp<cpp11::decay_t<int const &>>(ncol), cpp11::as_cpp<cpp11::decay_t<int const &>>(threads)));
  END_CPP11
}
// cpp11_sparsemat.cpp
extern cpp11::writable::list cpp11_sp64_transpose(cpp11::doubles const & x, cpp11::integers const & i, cpp11::doubles const & p, int const & nrow, int const & ncol, int const & threads);
extern "C" SEXP _fastde_cpp11_sp64_transpose(SEXP x, SEXP i, SEXP p, SEXP nrow, SEXP ncol, SEXP threads) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp11_sp64_transpose(cpp11::as_cpp<cpp11::decay_t<cpp11::doubles const &>>(x), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(i), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles const &>>(p), cpp11::as_cpp<cpp11::decay_t<int const &>>(nrow), cpp11::as_cpp<cpp11::decay_t<int const &>>(ncol), cpp11::as_cpp<cpp11::decay_t<int const &>>(threads)));
  END_CPP11
}
// cpp11_sparsemat.cpp
extern cpp11::writable::doubles_matrix<cpp11::by_column> cpp11_sp_to_dense(cpp11::doubles const & x, cpp11::integers const & i, cpp11::integers const & p, int const & nrow, int const & ncol, int const & threads);
extern "C" SEXP _fastde_cpp11_sp_to_dense(SEXP x, SEXP i, SEXP p, SEXP nrow, SEXP ncol, SEXP threads) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp11_sp_to_dense(cpp11::as_cpp<cpp11::decay_t<cpp11::doubles const &>>(x), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(i), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(p), cpp11::as_cpp<cpp11::decay_t<int const &>>(nrow), cpp11::as_cpp<cpp11::decay_t<int const &>>(ncol), cpp11::as_cpp<cpp11::decay_t<int const &>>(threads)));
  END_CPP11
}
// cpp11_sparsemat.cpp
extern cpp11::writable::doubles_matrix<cpp11::by_column> cpp11_sp64_to_dense(cpp11::doubles const & x, cpp11::integers const & i, cpp11::doubles const & p, int const & nrow, int const & ncol, int const & threads);
extern "C" SEXP _fastde_cpp11_sp64_to_dense(SEXP x, SEXP i, SEXP p, SEXP nrow, SEXP ncol, SEXP threads) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp11_sp64_to_dense(cpp11::as_cpp<cpp11::decay_t<cpp11::doubles const &>>(x), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(i), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles const &>>(p), cpp11::as_cpp<cpp11::decay_t<int const &>>(nrow), cpp11::as_cpp<cpp11::decay_t<int const &>>(ncol), cpp11::as_cpp<cpp11::decay_t<int const &>>(threads)));
  END_CPP11
}
// cpp11_sparsemat.cpp
extern cpp11::writable::doubles_matrix<cpp11::by_column> cpp11_sp_to_dense_transposed(cpp11::doubles const & x, cpp11::integers const & i, cpp11::integers const & p, int const & nrow, int const & ncol, int const & threads);
extern "C" SEXP _fastde_cpp11_sp_to_dense_transposed(SEXP x, SEXP i, SEXP p, SEXP nrow, SEXP ncol, SEXP threads) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp11_sp_to_dense_transposed(cpp11::as_cpp<cpp11::decay_t<cpp11::doubles const &>>(x), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(i), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(p), cpp11::as_cpp<cpp11::decay_t<int const &>>(nrow), cpp11::as_cpp<cpp11::decay_t<int const &>>(ncol), cpp11::as_cpp<cpp11::decay_t<int const &>>(threads)));
  END_CPP11
}
// cpp11_sparsemat.cpp
extern cpp11::writable::doubles_matrix<cpp11::by_column> cpp11_sp64_to_dense_transposed(cpp11::doubles const & x, cpp11::integers const & i, cpp11::doubles const & p, int const & nrow, int const & ncol, int const & threads);
extern "C" SEXP _fastde_cpp11_sp64_to_dense_transposed(SEXP x, SEXP i, SEXP p, SEXP nrow, SEXP ncol, SEXP threads) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp11_sp64_to_dense_transposed(cpp11::as_cpp<cpp11::decay_t<cpp11::doubles const &>>(x), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(i), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles const &>>(p), cpp11::as_cpp<cpp11::decay_t<int const &>>(nrow), cpp11::as_cpp<cpp11::decay_t<int const &>>(ncol), cpp11::as_cpp<cpp11::decay_t<int const &>>(threads)));
  END_CPP11
}
// cpp11_ttest.cpp
extern cpp11::sexp cpp11_dense_ttest(cpp11::doubles_matrix<cpp11::by_column> const & input, cpp11::strings const & features, cpp11::integers const & labels, int alternative, bool var_equal, bool as_dataframe, int threads);
extern "C" SEXP _fastde_cpp11_dense_ttest(SEXP input, SEXP features, SEXP labels, SEXP alternative, SEXP var_equal, SEXP as_dataframe, SEXP threads) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp11_dense_ttest(cpp11::as_cpp<cpp11::decay_t<cpp11::doubles_matrix<cpp11::by_column> const &>>(input), cpp11::as_cpp<cpp11::decay_t<cpp11::strings const &>>(features), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(labels), cpp11::as_cpp<cpp11::decay_t<int>>(alternative), cpp11::as_cpp<cpp11::decay_t<bool>>(var_equal), cpp11::as_cpp<cpp11::decay_t<bool>>(as_dataframe), cpp11::as_cpp<cpp11::decay_t<int>>(threads)));
  END_CPP11
}
// cpp11_ttest.cpp
extern cpp11::sexp cpp11_sparse_ttest(cpp11::doubles const & x, cpp11::integers const & i, cpp11::integers const & p, cpp11::strings const & features, int const & rows, int const & cols, cpp11::integers const & labels, bool features_as_rows, int alternative, bool var_equal, bool as_dataframe, int threads);
extern "C" SEXP _fastde_cpp11_sparse_ttest(SEXP x, SEXP i, SEXP p, SEXP features, SEXP rows, SEXP cols, SEXP labels, SEXP features_as_rows, SEXP alternative, SEXP var_equal, SEXP as_dataframe, SEXP threads) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp11_sparse_ttest(cpp11::as_cpp<cpp11::decay_t<cpp11::doubles const &>>(x), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(i), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(p), cpp11::as_cpp<cpp11::decay_t<cpp11::strings const &>>(features), cpp11::as_cpp<cpp11::decay_t<int const &>>(rows), cpp11::as_cpp<cpp11::decay_t<int const &>>(cols), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(labels), cpp11::as_cpp<cpp11::decay_t<bool>>(features_as_rows), cpp11::as_cpp<cpp11::decay_t<int>>(alternative), cpp11::as_cpp<cpp11::decay_t<bool>>(var_equal), cpp11::as_cpp<cpp11::decay_t<bool>>(as_dataframe), cpp11::as_cpp<cpp11::decay_t<int>>(threads)));
  END_CPP11
}
// cpp11_ttest.cpp
extern cpp11::sexp cpp11_sparse64_ttest(cpp11::doubles const & x, cpp11::integers const & i, cpp11::doubles const & p, cpp11::strings const & features, int const & rows, int const & cols, cpp11::integers const & labels, bool features_as_rows, int alternative, bool var_equal, bool as_dataframe, int threads);
extern "C" SEXP _fastde_cpp11_sparse64_ttest(SEXP x, SEXP i, SEXP p, SEXP features, SEXP rows, SEXP cols, SEXP labels, SEXP features_as_rows, SEXP alternative, SEXP var_equal, SEXP as_dataframe, SEXP threads) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp11_sparse64_ttest(cpp11::as_cpp<cpp11::decay_t<cpp11::doubles const &>>(x), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(i), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles const &>>(p), cpp11::as_cpp<cpp11::decay_t<cpp11::strings const &>>(features), cpp11::as_cpp<cpp11::decay_t<int const &>>(rows), cpp11::as_cpp<cpp11::decay_t<int const &>>(cols), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(labels), cpp11::as_cpp<cpp11::decay_t<bool>>(features_as_rows), cpp11::as_cpp<cpp11::decay_t<int>>(alternative), cpp11::as_cpp<cpp11::decay_t<bool>>(var_equal), cpp11::as_cpp<cpp11::decay_t<bool>>(as_dataframe), cpp11::as_cpp<cpp11::decay_t<int>>(threads)));
  END_CPP11
}
// cpp11_wmwtest.cpp
extern cpp11::sexp cpp11_dense_wmw(cpp11::doubles_matrix<cpp11::by_column> const & input, cpp11::strings const & features, cpp11::integers const & labels, int rtype, bool continuity_correction, bool as_dataframe, int threads);
extern "C" SEXP _fastde_cpp11_dense_wmw(SEXP input, SEXP features, SEXP labels, SEXP rtype, SEXP continuity_correction, SEXP as_dataframe, SEXP threads) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp11_dense_wmw(cpp11::as_cpp<cpp11::decay_t<cpp11::doubles_matrix<cpp11::by_column> const &>>(input), cpp11::as_cpp<cpp11::decay_t<cpp11::strings const &>>(features), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(labels), cpp11::as_cpp<cpp11::decay_t<int>>(rtype), cpp11::as_cpp<cpp11::decay_t<bool>>(continuity_correction), cpp11::as_cpp<cpp11::decay_t<bool>>(as_dataframe), cpp11::as_cpp<cpp11::decay_t<int>>(threads)));
  END_CPP11
}
// cpp11_wmwtest.cpp
extern cpp11::sexp cpp11_sparse_wmw(cpp11::doubles const & x, cpp11::integers const & i, cpp11::integers const & p, cpp11::strings const & features, int const & rows, int const & cols, cpp11::integers const & labels, bool features_as_rows, int rtype, bool continuity_correction, bool as_dataframe, int threads);
extern "C" SEXP _fastde_cpp11_sparse_wmw(SEXP x, SEXP i, SEXP p, SEXP features, SEXP rows, SEXP cols, SEXP labels, SEXP features_as_rows, SEXP rtype, SEXP continuity_correction, SEXP as_dataframe, SEXP threads) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp11_sparse_wmw(cpp11::as_cpp<cpp11::decay_t<cpp11::doubles const &>>(x), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(i), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(p), cpp11::as_cpp<cpp11::decay_t<cpp11::strings const &>>(features), cpp11::as_cpp<cpp11::decay_t<int const &>>(rows), cpp11::as_cpp<cpp11::decay_t<int const &>>(cols), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(labels), cpp11::as_cpp<cpp11::decay_t<bool>>(features_as_rows), cpp11::as_cpp<cpp11::decay_t<int>>(rtype), cpp11::as_cpp<cpp11::decay_t<bool>>(continuity_correction), cpp11::as_cpp<cpp11::decay_t<bool>>(as_dataframe), cpp11::as_cpp<cpp11::decay_t<int>>(threads)));
  END_CPP11
}
// cpp11_wmwtest.cpp
extern cpp11::sexp cpp11_sparse64_wmw(cpp11::doubles const & x, cpp11::integers const & i, cpp11::doubles const & p, cpp11::strings const & features, int const & rows, int const & cols, cpp11::integers const & labels, bool features_as_rows, int rtype, bool continuity_correction, bool as_dataframe, int threads);
extern "C" SEXP _fastde_cpp11_sparse64_wmw(SEXP x, SEXP i, SEXP p, SEXP features, SEXP rows, SEXP cols, SEXP labels, SEXP features_as_rows, SEXP rtype, SEXP continuity_correction, SEXP as_dataframe, SEXP threads) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp11_sparse64_wmw(cpp11::as_cpp<cpp11::decay_t<cpp11::doubles const &>>(x), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(i), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles const &>>(p), cpp11::as_cpp<cpp11::decay_t<cpp11::strings const &>>(features), cpp11::as_cpp<cpp11::decay_t<int const &>>(rows), cpp11::as_cpp<cpp11::decay_t<int const &>>(cols), cpp11::as_cpp<cpp11::decay_t<cpp11::integers const &>>(labels), cpp11::as_cpp<cpp11::decay_t<bool>>(features_as_rows), cpp11::as_cpp<cpp11::decay_t<int>>(rtype), cpp11::as_cpp<cpp11::decay_t<bool>>(continuity_correction), cpp11::as_cpp<cpp11::decay_t<bool>>(as_dataframe), cpp11::as_cpp<cpp11::decay_t<int>>(threads)));
  END_CPP11
}

extern "C" {
static const R_CallMethodDef CallEntries[] = {
    {"_fastde_cpp11_ComputeFoldChange",         (DL_FUNC) &_fastde_cpp11_ComputeFoldChange,         12},
    {"_fastde_cpp11_ComputeFoldChangeSparse",   (DL_FUNC) &_fastde_cpp11_ComputeFoldChangeSparse,   17},
    {"_fastde_cpp11_ComputeFoldChangeSparse64", (DL_FUNC) &_fastde_cpp11_ComputeFoldChangeSparse64, 17},
    {"_fastde_cpp11_FilterFoldChange",          (DL_FUNC) &_fastde_cpp11_FilterFoldChange,          10},
    {"_fastde_cpp11_FilterFoldChangeMat",       (DL_FUNC) &_fastde_cpp11_FilterFoldChangeMat,       10},
    {"_fastde_cpp11_dense_ttest",               (DL_FUNC) &_fastde_cpp11_dense_ttest,                7},
    {"_fastde_cpp11_dense_wmw",                 (DL_FUNC) &_fastde_cpp11_dense_wmw,                  7},
    {"_fastde_cpp11_sp64_to_dense",             (DL_FUNC) &_fastde_cpp11_sp64_to_dense,              6},
    {"_fastde_cpp11_sp64_to_dense_transposed",  (DL_FUNC) &_fastde_cpp11_sp64_to_dense_transposed,   6},
    {"_fastde_cpp11_sp64_transpose",            (DL_FUNC) &_fastde_cpp11_sp64_transpose,             6},
    {"_fastde_cpp11_sp_to_dense",               (DL_FUNC) &_fastde_cpp11_sp_to_dense,                6},
    {"_fastde_cpp11_sp_to_dense_transposed",    (DL_FUNC) &_fastde_cpp11_sp_to_dense_transposed,     6},
    {"_fastde_cpp11_sp_transpose",              (DL_FUNC) &_fastde_cpp11_sp_transpose,               6},
    {"_fastde_cpp11_sparse64_ttest",            (DL_FUNC) &_fastde_cpp11_sparse64_ttest,            12},
    {"_fastde_cpp11_sparse64_wmw",              (DL_FUNC) &_fastde_cpp11_sparse64_wmw,              12},
    {"_fastde_cpp11_sparse_ttest",              (DL_FUNC) &_fastde_cpp11_sparse_ttest,              12},
    {"_fastde_cpp11_sparse_wmw",                (DL_FUNC) &_fastde_cpp11_sparse_wmw,                12},
    {NULL, NULL, 0}
};
}

extern "C" attribute_visible void R_init_fastde(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
