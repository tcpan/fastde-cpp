#pragma once

// ------- function declaration

#include <Rcpp.h>
#include <vector>

Rcpp::StringVector rmatrix_to_vector(SEXP matrix, std::vector<bool> & mat,
    size_t & nrow, size_t & ncol, size_t & nelem);

Rcpp::StringVector rmatrix_to_vector(SEXP matrix, std::vector<double> & mat,
    size_t & nrow, size_t & ncol, size_t & nelem);

void rvector_to_vector(SEXP vect, std::vector<bool> & vec);

void rvector_to_vector(SEXP vect, std::vector<double> & vec);

void rvector_to_vector(SEXP vect, std::vector<int> & vec);

size_t rvector_to_vector(SEXP vect, std::vector<bool> & vec, size_t const & length, size_t const & offset = 0);

size_t rvector_to_vector(SEXP vect, std::vector<double> & vec, size_t const & length, size_t const & offset = 0);

size_t rvector_to_vector(SEXP vect, std::vector<int> & vec, size_t const & length, size_t const & offset = 0);


Rcpp::StringVector rsparsematrix_to_vectors(SEXP matrix, 
    std::vector<double> & x,
    std::vector<int> & i,
    std::vector<int> & p,
    size_t & nrow, size_t & ncol, size_t & nelem);

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
    SEXP min_diff_pct, SEXP logfc_threshold,
    SEXP only_pos, double _min_pct, double _min_diff_pct, double _logfc_thresh, bool _only_pos);


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
