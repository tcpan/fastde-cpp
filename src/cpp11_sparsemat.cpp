#include "utils_sparsemat.hpp"
#include "cpp11/doubles.hpp"
#include "cpp11/integers.hpp"
#include "cpp11/list.hpp"
#include "cpp11/matrix.hpp"

#include "cpp11/r_vector.hpp"
#include "cpp11/list_of.hpp"
#include <R.h>
#include <vector>

// adapters for cpp11 r_vector to c++ implementation

[[cpp11::register]]
extern cpp11::writable::list cpp11_sp_transpose(cpp11::doubles const & x,
    cpp11::integers const & i, cpp11::integers const & p, int const & nrow, int const & ncol, int const & threads) {
    if (threads == 1) return _sp_transpose(x, i, p, nrow, ncol, threads);
    else return _sp_transpose_par(x, i, p, nrow, ncol, threads);
}


[[cpp11::register]]
extern cpp11::writable::list cpp11_sp64_transpose(cpp11::doubles const & x,
    cpp11::integers const & i, cpp11::doubles const & p, int const & nrow, int const & ncol, int const & threads) {
    if (threads == 1) return _sp_transpose(x, i, p, nrow, ncol, threads);
    else return _sp_transpose_par(x, i, p, nrow, ncol, threads);
}



[[cpp11::register]]
extern cpp11::writable::doubles_matrix<cpp11::by_column> cpp11_sp_to_dense(cpp11::doubles const & x,
    cpp11::integers const & i, cpp11::integers const & p, int const & nrow, int const & ncol, int const & threads) {
    return _sp_to_dense<cpp11::writable::doubles_matrix<cpp11::by_column>>(x, i, p, nrow, ncol, threads);
}


[[cpp11::register]]
extern cpp11::writable::doubles_matrix<cpp11::by_column> cpp11_sp64_to_dense(cpp11::doubles const & x,
    cpp11::integers const & i, cpp11::doubles const & p, int const & nrow, int const & ncol, int const & threads) {
    return _sp_to_dense<cpp11::writable::doubles_matrix<cpp11::by_column>>(x, i, p, nrow, ncol, threads);
}



[[cpp11::register]]
extern cpp11::writable::doubles_matrix<cpp11::by_column> cpp11_sp_to_dense_transposed(
    cpp11::doubles const & x,
    cpp11::integers const & i, 
    cpp11::integers const & p, 
    int const & nrow, int const & ncol,
    int const & threads
) {
    return _sp_to_dense_transposed<cpp11::writable::doubles_matrix<cpp11::by_column>>(x, i, p, nrow, ncol, threads);
}

[[cpp11::register]]
extern cpp11::writable::doubles_matrix<cpp11::by_column> cpp11_sp64_to_dense_transposed(
    cpp11::doubles const & x,
    cpp11::integers const & i, 
    cpp11::doubles const & p, 
    int const & nrow, int const & ncol,
    int const & threads
) {
    return _sp_to_dense_transposed<cpp11::writable::doubles_matrix<cpp11::by_column>>(x, i, p, nrow, ncol, threads);
}


// this is a test to see how to use SEXP directly.
// namespace cpp11 {

// template <typename T>
// typename std::enable_if<std::is_same<T, std::vector<cpp11::r_vector<double>>>::value, T>::type as_cpp(SEXP from) {
//     T out;
//     if (TYPEOF(from) == VECSXP) {
//         R_xlen_t nlists = XLENGTH(from);
//         out.reserve(nlists);

//         for (R_xlen_t i = 0; i < nlists; ++i) {
//             // get list items in "from", convert to r_vector, and insert into output
//             out.emplace_back(cpp11::r_vector<double>(VECTOR_ELT(from, i)));
//         }
//     }
//     return out;
// }

// template <typename T>
// typename std::enable_if<std::is_same<T, std::vector<cpp11::r_vector<int>>>::value, T>::type as_cpp(SEXP from) {
//     T out;
//     if (TYPEOF(from) == VECSXP) {
//         R_xlen_t nlists = XLENGTH(from);
//         out.reserve(nlists);

//         for (R_xlen_t i = 0; i < nlists; ++i) {
//             // get list items in "from", convert to r_vector, and insert into output
//             out.emplace_back(cpp11::r_vector<int>(VECTOR_ELT(from, i)));
//         }
//     }
//     return out;
// }

// }


// this is much easier way to map.
// question - do we want to merge THEN compute using existing APIs?   would be easier....

// THIS IS A GOOD COMPILE TEST.
// [[cpp11::register]]
// extern int merge_sparse_mats(
//     cpp11::list_of<cpp11::doubles> const & xvecs,
//     cpp11::list_of<cpp11::integers> const & ivecs,
//     cpp11::list_of<cpp11::doubles> const & pvecs,
//     cpp11::integers const & nrows,
//     cpp11::integers const & ncols
// ) {
//     // have to do this because list_of.hpp's operator[] are not set up to work on const, and would throw an error.
//     // return cpp11::as_integers(ivecs.at(0))[0];

//     // below requries const operator[]
//     return ivecs[0][0];
// }

// "list_of" does note react well to "const &" - reports SEXPREC declared but not defined.
// this is because list is an r_vector<SEXP> and SEXP is an opaque pointer to SEXPREC, so it's not know that it is doubles array or integers array.
// data is in column major
[[cpp11::register]]
extern cpp11::writable::list cpp11_sp_rbind(
    cpp11::list_of<cpp11::doubles> const & xvecs,
    cpp11::list_of<cpp11::integers> const & ivecs,
    cpp11::list_of<cpp11::integers> const & pvecs,
    cpp11::integers const & nrows,
    cpp11::integers const & ncols,
    int const & threads
) {
    return _sp_rbind(xvecs, ivecs, pvecs, nrows, ncols, threads);
}

[[cpp11::register]]
extern cpp11::writable::list cpp11_sp64_rbind(
    cpp11::list_of<cpp11::doubles> const & xvecs,
    cpp11::list_of<cpp11::integers> const & ivecs,
    cpp11::list_of<cpp11::doubles> const & pvecs,
    cpp11::integers const & nrows,
    cpp11::integers const & ncols,
    int const & threads
) {
    return _sp_rbind(xvecs, ivecs, pvecs, nrows, ncols, threads);
}

[[cpp11::register]]
extern cpp11::writable::list cpp11_sp_cbind(
    cpp11::list_of<cpp11::doubles> const & xvecs,
    cpp11::list_of<cpp11::integers> const & ivecs,
    cpp11::list_of<cpp11::integers> const & pvecs,
    cpp11::integers const & nrows,
    cpp11::integers const & ncols,
    int const & threads
) {
    return _sp_cbind(xvecs, ivecs, pvecs, nrows, ncols, threads);
}

[[cpp11::register]]
extern cpp11::writable::list cpp11_sp64_cbind(
    cpp11::list_of<cpp11::doubles> const & xvecs,
    cpp11::list_of<cpp11::integers> const & ivecs,
    cpp11::list_of<cpp11::doubles> const & pvecs,
    cpp11::integers const & nrows,
    cpp11::integers const & ncols,
    int const & threads
) {
    return _sp_cbind(xvecs, ivecs, pvecs, nrows, ncols, threads);
}




[[cpp11::register]]
extern cpp11::writable::doubles cpp11_sp_colSums(
    cpp11::doubles const & x,
    cpp11::integers const & p,
    int const & threads
) {
    return _sp_colsums(x, p, static_cast<int>(p.size() - 1), threads);
}

[[cpp11::register]]
extern cpp11::writable::doubles cpp11_sp64_colSums(
    cpp11::doubles const & x,
    cpp11::doubles const & p,
    int const & threads
) {
    return _sp_colsums(x, p, static_cast<int>(p.size() - 1), threads);
}


[[cpp11::register]]
extern cpp11::writable::doubles cpp11_sp_rowSums(
    cpp11::doubles const & x,
    cpp11::integers const & i,
    int const & nrows,
    int const & threads
) {
    return _sp_rowsums(x, i, nrows, static_cast<size_t>(x.size()), threads);
}
