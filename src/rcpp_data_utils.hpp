#pragma once

// ------- function declaration

#include <Rcpp.h>
#include <vector>

#include "fastde_types.hpp"



// from https://gallery.rcpp.org/articles/sparse-matrix-class/
// namespace Rcpp {

//     // specialization of Rcpp::as
//     template <> dgCMatrix as(SEXP mat);

//     // specialization of Rcpp::wrap
//     template <> SEXP wrap(const dgCMatrix& sm);


//     // NOTE : spam is row major. here we are column major, so it's a transpose of original.

//     // no specialization of Rcpp::as - no way to 
//     template <> spam32 as(SEXP mat);
//     template <> spam64 as(SEXP mat);

//     // // specialization of Rcpp::wrap
//     // template <> SEXP wrap(const spam64& sm);
//     // specialization of Rcpp::wrap
//     template <> SEXP wrap(const spam64& sm);
//     template <> SEXP wrap(const spam32& sm);

// }

//' round trip testing for dgCMatrix
//'
//' @rdname rttest_dgCMatrix
//' @param mat an R sparse compressed column matrix
//' @return R sparse compressed column matrix
//' @name rttest_dgCMatrix
//' @export
//[[Rcpp::export]]
Rcpp::dgCMatrix rttest_dgCMatrix(Rcpp::dgCMatrix& mat);

//' round trip testing for spam32
//'
//' @rdname rttest_spam32
//' @param mat an spam sparse compressed column matrix
//' @return spam sparse compressed column matrix
//' @name rttest_spam32
//' @export
//[[Rcpp::export]]
Rcpp::spam32 rttest_spam32(Rcpp::spam32& mat);

//' round trip testing for spam64
//'
//' @rdname rttest_spam64
//' @param mat an spam sparse compressed column matrix
//' @return spam sparse compressed column matrix
//' @name rttest_spam64
//' @export
//[[Rcpp::export]]
Rcpp::spam64 rttest_spam64(Rcpp::spam64& mat);


Rcpp::StringVector copy_rmatrix_to_cppvector(Rcpp::LogicalMatrix _matrix, std::vector<bool> & mat,
    size_t & nrow, size_t & ncol, size_t & nelem);

Rcpp::StringVector copy_rmatrix_to_cppvector(Rcpp::NumericMatrix _matrix, std::vector<double> & mat,
    size_t & nrow, size_t & ncol, size_t & nelem);

size_t copy_rvector_to_cppvector(Rcpp::LogicalVector _vector, std::vector<bool> & vec, 
    size_t const & length = std::numeric_limits<size_t>::max(), size_t const & offset = 0);

size_t copy_rvector_to_cppvector(Rcpp::NumericVector _vector, std::vector<double> & vec, 
    size_t const & length = std::numeric_limits<size_t>::max(), size_t const & offset = 0);

size_t copy_rvector_to_cppvector(Rcpp::IntegerVector _vector, std::vector<int> & vec, 
    size_t const & length = std::numeric_limits<size_t>::max(), size_t const & offset = 0);

size_t copy_rvector_to_cppvector(Rcpp::NumericVector _vector, std::vector<long> & vec, 
    size_t const & length = std::numeric_limits<size_t>::max(), size_t const & offset = 0);


Rcpp::StringVector copy_rsparsematrix_to_cppvectors(Rcpp::dgCMatrix matrix, 
    std::vector<double> & x,
    std::vector<int> & i,
    std::vector<int> & p,
    size_t & nrow, size_t & ncol, size_t & nelem);


Rcpp::StringVector copy_rsparsematrix_to_cppvectors(Rcpp::spam64 matrix, 
    std::vector<double> & x,
    std::vector<long> & i,
    std::vector<long> & p,
    size_t & nrow, size_t & ncol, size_t & nelem);

Rcpp::StringVector copy_rsparsematrix_to_cppvectors(Rcpp::spam32 matrix, 
    std::vector<double> & x,
    std::vector<int> & i,
    std::vector<int> & p,
    size_t & nrow, size_t & ncol, size_t & nelem);

Rcpp::StringVector copy_rsparsematrix_to_cppvectors(
    Rcpp::NumericVector _x, Rcpp::IntegerVector _i, Rcpp::IntegerVector _p, 
    std::vector<double> & x,
    std::vector<int> & i,
    std::vector<int> & p,
    size_t & nrow, size_t & ncol, size_t & nelem);


Rcpp::StringVector copy_rsparsematrix_to_cppvectors(
    Rcpp::NumericVector _x, Rcpp::NumericVector _i, Rcpp::NumericVector _p, 
    std::vector<double> & x,
    std::vector<long> & i,
    std::vector<long> & p,
    size_t & nrow, size_t & ncol, size_t & nelem);


//------ these may not be needed - Rcpp can wrap?
void import_r_common_params(SEXP as_dataframe, SEXP threads,
    bool & _as_dataframe, int & nthreads);
void import_de_common_params(SEXP rtype,
    SEXP bool_param, int & type, bool & bool_val);
void import_fc_common_params(SEXP calc_percents, 
    SEXP min_threshold, SEXP use_expm1,
    SEXP use_log, SEXP log_base, SEXP use_pseudocount,
    bool & perc, double & min_thresh, bool & _use_expm1, 
    bool & _use_log, double & _log_base, bool & _use_pseudocount);

void import_filterfc_common_params(SEXP min_pct, 
    SEXP min_diff_pct, SEXP logfc_threshold, SEXP only_pos, 
    double & _min_pct, double & _min_diff_pct, double & _logfc_thresh, bool & _only_pos);
//-----------------

Rcpp::StringVector populate_feature_names(Rcpp::StringVector const & features,
    size_t const & nfeatures);

Rcpp::NumericMatrix export_de_to_r_matrix(
    std::vector<double> const & pv,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    Rcpp::StringVector const & features
);

Rcpp::DataFrame export_de_to_r_dataframe(
    std::vector<double> const & pv, std::string const & name,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    Rcpp::StringVector const & features
);

Rcpp::DataFrame export_fc_to_r_dataframe(
    std::vector<double> const & fc, std::string const & fcname,
    std::vector<double> const & p1, std::string const & p1name,
    std::vector<double> const & p2, std::string const & p2name,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    Rcpp::StringVector const & features
);
Rcpp::List export_fc_to_r_matrix(
    std::vector<double> const & fc, std::string const & fcname,
    std::vector<double> const & p1, std::string const & p1name,
    std::vector<double> const & p2, std::string const & p2name,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    Rcpp::StringVector const & features
);
Rcpp::List export_fc_to_r_matrix(
    std::vector<double> const & fc, std::string const & fcname,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    Rcpp::StringVector const & features
);
