#pragma once

// ------- function declaration

#include <Rcpp.h>
#include <vector>

Rcpp::StringVector rmatrix_to_vector(SEXP matrix, std::vector<double> & mat,
    size_t & nrow, size_t & ncol, size_t & nelem);

void rvector_to_vector(SEXP vect, std::vector<double> & vec);

void rvector_to_vector(SEXP vect, std::vector<int> & vec);

size_t rvector_to_vector(SEXP vect, std::vector<double> & vec, size_t const & length, size_t const & offset = 0);

size_t rvector_to_vector(SEXP vect, std::vector<int> & vec, size_t const & length, size_t const & offset = 0);


Rcpp::StringVector rsparsematrix_to_vectors(SEXP matrix, 
    std::vector<double> & x,
    std::vector<int> & i,
    std::vector<int> & p,
    size_t & nrow, size_t & ncol, size_t & nelem);

void import_r_common_params(SEXP alternative,
    SEXP var_equal, SEXP as_dataframe, SEXP threads,
    int & type, bool & var_eq, bool & _as_dataframe, int & nthreads);

Rcpp::StringVector populate_feature_names(Rcpp::StringVector const & features,
    size_t const & nfeatures);

Rcpp::NumericMatrix export_de_to_r_matrix(
    std::vector<double> const & pv,
    std::vector<int> const & sorted_labels,
    Rcpp::StringVector const & features
);

Rcpp::DataFrame export_de_to_r_dataframe(
    std::vector<double> const & pv,
    std::vector<int> const & sorted_labels,
    Rcpp::StringVector const & features
);