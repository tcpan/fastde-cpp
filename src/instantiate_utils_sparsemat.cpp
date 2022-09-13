#include "utils_sparsemat.tpp"
#include "cpp11/doubles.hpp"
#include "cpp11/integers.hpp"



// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// template cpp11::writable::list _sp_transpose_sort(
//     cpp11::doubles const & x, 
//     cpp11::integers const & i, 
//     cpp11::integers const & p, int const & nrow, int const & ncol, int const & threads);


// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// template cpp11::writable::list _sp_transpose_sort(
//     cpp11::doubles const & x, 
//     cpp11::integers const & i, 
//     cpp11::doubles const & p, int const & nrow, int const & ncol, int const & threads);


// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
template cpp11::writable::list _sp_transpose_par(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::integers const & p, 
    int const & nrow, int const & ncol, int const & threads);

// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
template cpp11::writable::list _sp_transpose_par(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::doubles const & p, 
    int const & nrow, int const & ncol, int const & threads);



// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
template cpp11::writable::list _sp_transpose(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::integers const & p, int const & nrow, int const & ncol, int const & threads) ;

// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
template cpp11::writable::list _sp_transpose(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::doubles const & p, int const & nrow, int const & ncol, int const & threads) ;

// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template void _sp_transpose_par(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::integers const & p, 
    int const & nrow, int const & ncol,
    std::vector<double> & tx, 
    std::vector<int> & ti, 
    std::vector<int> & tp, 
    int const & threads);

// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template void _sp_transpose_par(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::doubles const & p, 
    int const & nrow, int const & ncol,
    std::vector<double> & tx, 
    std::vector<int> & ti, 
    std::vector<long> & tp, 
    int const & threads);

// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template void _sp_transpose_par(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::integers const & p, 
    int const & nrow, int const & ncol,
    double * tx, 
    int * ti, 
    int * tp, 
    int const & threads);

// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template void _sp_transpose_par(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::doubles const & p, 
    int const & nrow, int const & ncol,
    double * tx, 
    int * ti, 
    long * tp, 
    int const & threads);


// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template void _sp_transpose(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::integers const & p, 
    int const & nrow, int const & ncol, 
    std::vector<double> & tx, 
    std::vector<int> & ti, 
    std::vector<int> & tp, 
    int const & threads);

// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template void _sp_transpose(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::doubles const & p, 
    int const & nrow, int const & ncol, 
    std::vector<double> & tx, 
    std::vector<int> & ti, 
    std::vector<long> & tp, 
    int const & threads);

// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template void _sp_transpose(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::integers const & p, 
    int const & nrow, int const & ncol, 
    double * tx, 
    int * ti, 
    int * tp, 
    int const & threads);

// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template void _sp_transpose(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::doubles const & p, 
    int const & nrow, int const & ncol, 
    double * tx, 
    int * ti, 
    long * tp, 
    int const & threads);

// no names.
template cpp11::writable::doubles_matrix<cpp11::by_column> _sp_to_dense(cpp11::doubles const & x, 
    cpp11::integers const & i, cpp11::integers const & p, 
    int const & nrow, int const & ncol, int const & threads);

template cpp11::writable::doubles_matrix<cpp11::by_column> _sp_to_dense(cpp11::doubles const & x, 
    cpp11::integers const & i, cpp11::doubles const & p, 
    int const & nrow, int const & ncol, int const & threads);

template cpp11::writable::doubles_matrix<cpp11::by_column> _sp_to_dense_transposed(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::integers const & p, int const & nrow, int const & ncol, int const & threads);

template cpp11::writable::doubles_matrix<cpp11::by_column> _sp_to_dense_transposed(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::doubles const & p, int const & nrow, int const & ncol, int const & threads);