#pragma once

// This file contains defined Rcpp types that are automatically included in RcppExports...
#include <Rcpp.h>

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
            dgCMatrix(int nrows, int ncols, int nelem);

            // column iterator -- see https://gallery.rcpp.org/articles/sparse-matrix-class/
        
    };



    // NOTE : spam is row major. here we are column major, so it's a transpose of original.
    // NOTE : spam has 1-based index, and dgCMatrix has 0-based index.

    class spamx32 {  // row major
        public:
            using INDEX = int;
        
            Rcpp::IntegerVector i;   // one-based? col id
            Rcpp::IntegerVector p;   // one-based, row pointers
            Rcpp::IntegerVector Dim;
            Rcpp::NumericVector x;
            Rcpp::List Dimnames;

            // constructor
            spamx32(Rcpp::S4 mat);
            spamx32(int nrows, int ncols, int nelem);
    };

    class spamx64 {   // row major
        public:
            using INDEX = long;
        
            Rcpp::NumericVector i;   // one-based? col id
            Rcpp::NumericVector p;   // one-based, row pointer
            Rcpp::NumericVector Dim;
            Rcpp::NumericVector x;
            Rcpp::List Dimnames;

            // constructor
            spamx64(Rcpp::S4 mat);
            spamx64(long nrows, long ncols, long nelem);
    };


    class spam32 {  // row major
        public:
            using INDEX = int;

            Rcpp::IntegerVector i;   // one-based? col id
            Rcpp::IntegerVector p;   // one-based, row pointer
            Rcpp::IntegerVector Dim;
            Rcpp::NumericVector x;

            // constructor
            spam32(Rcpp::S4 mat);
            spam32(int nrows, int ncols, int nelem);
    };

    class spam64 {  // row major
        public:
            using INDEX = long;
            
            Rcpp::NumericVector i;   // one-based? col id
            Rcpp::NumericVector p;   // one-based, row pointer.
            Rcpp::NumericVector Dim;
            Rcpp::NumericVector x;

            // constructor
            spam64(Rcpp::S4 mat);
            spam64(long nrows, long ncols, long nelem);
    };


    // specialization of Rcpp::as
    template <> dgCMatrix as(SEXP mat);

    // specialization of Rcpp::wrap
    template <> SEXP wrap(const dgCMatrix& sm);


    // NOTE : spam is row major. here we are column major, so it's a transpose of original.

    // no specialization of Rcpp::as - no way to 
    template <> spamx32 as(SEXP mat);
    template <> spamx64 as(SEXP mat);

    // // specialization of Rcpp::wrap
    // template <> SEXP wrap(const spam64& sm);
    // specialization of Rcpp::wrap
    template <> SEXP wrap(const spamx64& sm);
    template <> SEXP wrap(const spamx32& sm);


    // no specialization of Rcpp::as - no way to 
    template <> spam32 as(SEXP mat);
    template <> spam64 as(SEXP mat);

    // // specialization of Rcpp::wrap
    // template <> SEXP wrap(const spam64& sm);
    // specialization of Rcpp::wrap
    template <> SEXP wrap(const spam64& sm);
    template <> SEXP wrap(const spam32& sm);

}
