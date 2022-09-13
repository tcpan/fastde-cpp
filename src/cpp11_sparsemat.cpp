#include "utils_sparsemat.hpp"
#include "cpp11/doubles.hpp"
#include "cpp11/integers.hpp"
#include "cpp11/list.hpp"
#include "cpp11/matrix.hpp"


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