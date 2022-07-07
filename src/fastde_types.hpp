// This file contains defined Rcpp types that are automatically included in RcppExports...
#include <Rcpp.h>


// from https://gallery.rcpp.org/articles/sparse-matrix-class/
namespace Rcpp {
    class dgCMatrix {
    public:
        Rcpp::IntegerVector i, p, Dim;
        Rcpp::NumericVector x;
        Rcpp::List Dimnames;

        // constructor
        dgCMatrix(Rcpp::S4 mat);

        // column iterator -- see https://gallery.rcpp.org/articles/sparse-matrix-class/
        
    };



    // NOTE : spam is row major. here we are column major, so it's a transpose of original.

    class spam32 {
        public:
            Rcpp::IntegerVector i, p, Dim;
            Rcpp::NumericVector x;

            // constructor
            spam32(Rcpp::S4 mat);
    };

    class spam64 {
        public:
            Rcpp::NumericVector i, p, Dim;
            Rcpp::NumericVector x;

            // constructor
            spam64(Rcpp::S4 mat);
    };



    // specialization of Rcpp::as
    template <> dgCMatrix as(SEXP mat);

    // specialization of Rcpp::wrap
    template <> SEXP wrap(const dgCMatrix& sm);


    // NOTE : spam is row major. here we are column major, so it's a transpose of original.

    // no specialization of Rcpp::as - no way to 
    template <> spam32 as(SEXP mat);
    template <> spam64 as(SEXP mat);

    // // specialization of Rcpp::wrap
    // template <> SEXP wrap(const spam64& sm);
    // specialization of Rcpp::wrap
    template <> SEXP wrap(const spam64& sm);
    template <> SEXP wrap(const spam32& sm);


}
