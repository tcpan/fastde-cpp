#pragma once

// ------- function declaration

#include <Rcpp.h>
#include <vector>

#include "fastde_types.hpp"



//' round trip testing for dgCMatrix
//'
//' @rdname rttest_dgCMatrix
//' @param mat an R sparse compressed column matrix
//' @return R sparse compressed column matrix
//' @name rttest_dgCMatrix
//' @export
//[[Rcpp::export]]
Rcpp::dgCMatrix rttest_dgCMatrix(Rcpp::dgCMatrix const & mat);

//' round trip testing for dgCMatrix
//'
//' @rdname rttest_dgCMatrix64
//' @param mat an R sparse compressed column matrix
//' @return R sparse compressed column matrix
//' @name rttest_dgCMatrix64
//' @export
//[[Rcpp::export]]
Rcpp::dgCMatrix64 rttest_dgCMatrix64(Rcpp::dgCMatrix64 const & mat);




Rcpp::StringVector copy_rmatrix_to_cppvector(Rcpp::LogicalMatrix const & _matrix, std::vector<bool> & mat,
    size_t & nrow, size_t & ncol, size_t & nelem);

Rcpp::StringVector copy_rmatrix_to_cppvector(Rcpp::NumericMatrix const & _matrix, std::vector<double> & mat,
    size_t & nrow, size_t & ncol, size_t & nelem);

size_t copy_rvector_to_cppvector(Rcpp::LogicalVector const & _vector, std::vector<bool> & vec, 
    size_t const & length = std::numeric_limits<size_t>::max(), size_t const & offset = 0);

size_t copy_rvector_to_cppvector(Rcpp::NumericVector const & _vector, std::vector<double> & vec, 
    size_t const & length = std::numeric_limits<size_t>::max(), size_t const & offset = 0);

size_t copy_rvector_to_cppvector(Rcpp::IntegerVector const & _vector, std::vector<int> & vec, 
    size_t const & length = std::numeric_limits<size_t>::max(), size_t const & offset = 0);

size_t copy_rvector_to_cppvector(Rcpp::NumericVector const & _vector, std::vector<long> & vec, 
    size_t const & length = std::numeric_limits<size_t>::max(), size_t const & offset = 0);


Rcpp::StringVector copy_rsparsematrix_to_cppvectors(Rcpp::dgCMatrix const & matrix, 
    std::vector<double> & x,
    std::vector<int> & i,
    std::vector<int> & p,
    size_t & nrow, size_t & ncol, size_t & nelem);

Rcpp::StringVector copy_rsparsematrix_to_cppvectors(Rcpp::dgCMatrix64 const & matrix, 
    std::vector<double> & x,
    std::vector<long> & i,
    std::vector<long> & p,
    size_t & nrow, size_t & ncol, size_t & nelem);


Rcpp::StringVector copy_rsparsematrix_to_cppvectors(
    Rcpp::NumericVector const & _x, Rcpp::IntegerVector const & _i, Rcpp::IntegerVector const & _p, 
    std::vector<double> & x,
    std::vector<int> & i,
    std::vector<int> & p,
    size_t & nrow, size_t & ncol, size_t & nelem);


Rcpp::StringVector copy_rsparsematrix_to_cppvectors(
    Rcpp::NumericVector const & _x, Rcpp::NumericVector const & _i, Rcpp::NumericVector const & _p, 
    std::vector<double> & x,
    std::vector<long> & i,
    std::vector<long> & p,
    size_t & nrow, size_t & ncol, size_t & nelem);


// //------ these may not be needed - Rcpp can wrap?
// void import_r_common_params(SEXP as_dataframe, SEXP threads,
//     bool & _as_dataframe, int & nthreads);
// void import_de_common_params(SEXP rtype,
//     SEXP bool_param, int & type, bool & bool_val);
// void import_fc_common_params(SEXP calc_percents, 
//     SEXP min_threshold, SEXP use_expm1,
//     SEXP use_log, SEXP log_base, SEXP use_pseudocount,
//     bool & perc, double & min_thresh, bool & _use_expm1, 
//     bool & _use_log, double & _log_base, bool & _use_pseudocount);

// void import_filterfc_common_params(SEXP min_pct, 
//     SEXP min_diff_pct, SEXP logfc_threshold, SEXP only_pos, 
//     double & _min_pct, double & _min_diff_pct, double & _logfc_thresh, bool & _only_pos);
// //-----------------

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

