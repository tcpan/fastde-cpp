#pragma once 

// pure c/c++ versions

// Todo: remove all cpp11 references.  actually, using ITER makes it slower?  move this code to fastde core for now.

#include <vector>

/*
 * wrapper for R dgCMatrix
 *
 */


// transpose by sorting VERY SLOW
// // NOTe:  there is no formal definition of sparse matrix.
// // input is column major, so i has the row ids, and p is per column.
// template <typename XITER, typename IITER, typename PITER, typename XITER2, typename IITER2, typename PITER2>
// extern void _sp_transpose_sort(
//     XITER x, 
//     IITER i, 
//     PITER p, 
//     IT2 const & nrow, IT2 const & ncol, 
    // XITER2 tx, 
    // IITER2 ti, 
    // PITER2 tp, 
    // int const & threads);

// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template <typename XITER, typename IITER, typename PITER, typename XITER2, typename IITER2, typename PITER2>
extern void _sp_transpose_par(
    XITER x, 
    IITER i, 
    PITER p, 
    IT2 const & nrow, IT2 const & ncol,
    XITER2 tx, 
    IITER2 ti, 
    PITER2 tp, 
    int const & threads);


// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template <typename XITER, typename IITER, typename PITER, 
    typename XITER2, typename IITER2, typename PITER2>
extern void _sp_transpose(
    XITER x, 
    IITER i, 
    PITER p, 
    IT2 const & nrow, IT2 const & ncol,
    XITER2 tx, 
    IITER2 ti, 
    PITER2 tp, 
    int const & threads);

// no names.
template <typename OITER, typename XITER, typename IITER, typename PITER, typename IT2>
extern void _sp_to_dense(
    XITER x, 
    IITER i,
    PITER p, 
    IT2 const & nrow, IT2 const & ncol, 
    OITER out,
    int const & threads);

template <typename OITER, typename XITER, typename IITER, typename PITER, typename IT2>
extern OUT _sp_to_dense_transposed(
    XITER x, 
    IITER i, 
    PITER p, 
    IT2 const & nrow, IT2 const & ncol, 
    OITER out,
    int const & threads);