#pragma once

// from https://gallery.rcpp.org/articles/sparse-matrix-class/
namespace Rcpp {
    class dgCMatrix64 {  // column major
        public:
            using INDEX = long;

            Rcpp::NumericVector i;   // zero-based, row id
            Rcpp::NumericVector p;   // zero-based, column pointers.
            Rcpp::NumericVector dimension;
            Rcpp::IntegerVector Dim;
            Rcpp::NumericVector x;
            Rcpp::List Dimnames;

            // constructor
            dgCMatrix64(Rcpp::S4 mat);
            dgCMatrix64(long const & nrows, long const & ncols, long const & nelem);

            // column iterator -- see https://gallery.rcpp.org/articles/sparse-matrix-class/        
    };

    // specialization of Rcpp::as
    template <> dgCMatrix64 as(SEXP mat);

    // specialization of Rcpp::wrap
    template <> SEXP wrap(const dgCMatrix64& sm);

}
