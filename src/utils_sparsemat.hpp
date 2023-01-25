#pragma once 

#include "cpp11/r_vector.hpp"
#include "cpp11/list.hpp"
#include "cpp11/list_of.hpp"

// Todo: remove all cpp11 references.  actually, using ITER makes it slower?  move this code to fastde core for now.

#include <vector>

/*
 * wrapper for R dgCMatrix
 *
 */

// adapters for cpp11 r_vector to c++ implementation


// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
template <typename XT, typename IT, typename PT, typename IT2>
extern cpp11::writable::list _sp_transpose(
    cpp11::r_vector<XT> const & x, 
    cpp11::r_vector<IT> const & i, 
    cpp11::r_vector<PT> const & p, 
    IT2 const & nrow, IT2 const & ncol, int const & threads) ;


// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template <typename XT, typename IT, typename PT, typename IT2, typename PT2>
extern void _sp_transpose(
    cpp11::r_vector<XT> const & x, 
    cpp11::r_vector<IT> const & i, 
    cpp11::r_vector<PT> const & p, 
    IT2 const & nrow, IT2 const & ncol, 
    std::vector<XT> & tx, 
    std::vector<IT2> & ti, 
    std::vector<PT2> & tp, 
    int const & threads);

// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template <typename XT, typename IT, typename PT, typename IT2, typename PT2>
extern void _sp_transpose(
    cpp11::r_vector<XT> const & x, 
    cpp11::r_vector<IT> const & i, 
    cpp11::r_vector<PT> const & p, 
    IT2 const & nrow, IT2 const & ncol, 
    XT * tx, 
    IT2 * ti, 
    PT2 * tp, 
    int const & threads);

// no names.
template <typename OUT, typename XT, typename IT, typename PT, typename IT2>
extern OUT _sp_to_dense(
    cpp11::r_vector<XT> const & x, 
    cpp11::r_vector<IT> const & i, 
    cpp11::r_vector<PT> const & p, 
    IT2 const & nrow, IT2 const & ncol, int const & threads);

template <typename OUT, typename XT, typename IT, typename PT, typename IT2>
extern OUT _sp_to_dense_transposed(
    cpp11::r_vector<XT> const & x, 
    cpp11::r_vector<IT> const & i, 
    cpp11::r_vector<PT> const & p, 
    IT2 const & nrow, IT2 const & ncol, int const & threads);


template <typename XT, typename IT, typename PT>
extern cpp11::writable::list _sp_rbind(
    cpp11::list_of<cpp11::r_vector<XT>> const & xvecs, 
    cpp11::list_of<cpp11::r_vector<IT>> const & ivecs, 
    cpp11::list_of<cpp11::r_vector<PT>> const & pvecs, 
    cpp11::r_vector<IT> const & nrows, 
    cpp11::r_vector<IT> const & ncols, int const & threads);

template <typename XT, typename IT, typename PT>
extern cpp11::writable::list _sp_cbind(
    cpp11::list_of<cpp11::r_vector<XT>> const & xvecs, 
    cpp11::list_of<cpp11::r_vector<IT>> const & ivecs, 
    cpp11::list_of<cpp11::r_vector<PT>> const & pvecs, 
    cpp11::r_vector<IT> const & nrows, 
    cpp11::r_vector<IT> const & ncols, int const & threads);


// csc
template <typename XT, typename PT, typename IT>
extern cpp11::r_vector<XT> _sp_colsums(
    cpp11::r_vector<XT> const & x, 
    cpp11::r_vector<PT> const & p, 
    IT const & ncol, 
    int const & threads);

// csc
template <typename XT, typename IT, typename IT2>
extern cpp11::r_vector<XT> _sp_rowsums(
    cpp11::r_vector<XT> const & x, 
    cpp11::r_vector<IT> const & i, 
    IT const & nrow, IT2 const & nzcount, 
    int const & threads);
