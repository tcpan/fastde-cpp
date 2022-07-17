#pragma once

// from https://gallery.rcpp.org/articles/sparse-matrix-class/
namespace Rcpp {
    class dgCMatrix {  // column major
        public:
            using INDEX = int;

            Rcpp::IntegerVector i;   // zero-based, row id
            Rcpp::IntegerVector p;   // zero-based, column pointers.
            Rcpp::IntegerVector Dim;
            Rcpp::NumericVector x;
            Rcpp::List Dimnames;

            // constructor
            dgCMatrix(Rcpp::S4 mat);
            dgCMatrix(int const & nrows, int const & ncols, int const & nelem);

            // column iterator -- see https://gallery.rcpp.org/articles/sparse-matrix-class/
        
    };

    // specialization of Rcpp::as
    template <> dgCMatrix as(SEXP mat);

    // specialization of Rcpp::wrap
    template <> SEXP wrap(const dgCMatrix& sm);

}
