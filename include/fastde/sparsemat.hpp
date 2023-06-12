#pragma once 

// pure c/c++ versions.
// NOTE: using iterators is potentially because of multiple constructor calls for iterators.
//       may be better to leverage "container" templates with operator[] to speed up.
// when comparing using cpp11::writable::doubles::iterator, cpp11::doubles::const_iterator, 
//    std::vector::iterator, and containers directly.  2800 colsums.  R built in ; 0.011s.
//  1. cpp11 const_iterator for read and cpp11 iterator for write (with proxy) - 1T 0.017s, 4T 0.007s
//  2. cpp11 const_iterator for read, std vec iterator for write (and copy to cpp11 doubles):  1T 0.016s, 4T 0.008s
//  3. cpp11 doubles containers in and out:  1T 0.005s, 4T 0.004s  (0.002s for wrapper including col names copy)
//  trying to template the containers directly, and the same colsum:  1T: 0.005s 4T: 0.003s.  R 0.01s.

// USE containers.

// Todo: remove all cpp11 references.  
// actually, using ITER makes it slower?  move this code to fastde core for now.

#include <stddef.h>

#include <vector>

/*
 * wrapper for R dgCMatrix
 *
 */


// transpose by sorting VERY SLOW
// // NOTe:  there is no formal definition of sparse matrix.
// // input is column major, so i has the row ids, and p is per column.
// template <typename XITER, typename IITER, typename PITER, typename XITER2, typename IITER2, typename PITER2>
// extern void csc_transpose_sort(
//     XITER x, 
//     IITER i, 
//     PITER p, 
//     IT2 const & nrow, IT2 const & ncol, 
    // XITER2 tx, 
    // IITER2 ti, 
    // PITER2 tp, 
    // int const & threads);

// NOTE:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template <typename XITER, typename IITER, typename PITER, typename IT2,
    typename XITER2, typename IITER2, typename PITER2>
extern void csc_transpose_csc(
    XITER x, 
    IITER i, 
    PITER p, 
    IT2 const & nrow, IT2 const & ncol,
    XITER2 tx, 
    IITER2 ti, 
    PITER2 tp, 
    int const & threads);

template <typename XVEC, typename IVEC, typename PVEC, typename IT2,
    typename XVEC2, typename IVEC2, typename PVEC2>
extern void csc_transpose_csc_vec(
    XVEC const & x, 
    IVEC const & i, 
    PVEC const & p, 
    IT2 const & nrow, IT2 const & ncol,
    XVEC2 & tx, 
    IVEC2 & ti, 
    PVEC2 & tp, 
    int const & threads);



// no names.
template <typename OITER, typename XITER, typename IITER,
    typename PITER, typename IT2>
extern void csc_to_dense_c(
    XITER x, 
    IITER i,
    PITER p, 
    IT2 const & nrow, IT2 const & ncol, 
    OITER out,
    int const & threads);


// no names.
template <typename OVEC, typename XVEC, typename IVEC,
    typename PVEC, typename IT2>
extern void csc_to_dense_c_vec(
    XVEC const & x, 
    IVEC const & i,
    PVEC const & p, 
    IT2 const & nrow, IT2 const & ncol, 
    OVEC & out,
    int const & threads);

template <typename OITER, typename XITER, typename IITER,
    typename PITER, typename IT2>
extern void csc_to_dense_transposed_c(
    XITER x, 
    IITER i, 
    PITER p, 
    IT2 const & nrow, IT2 const & ncol, 
    OITER out,
    int const & threads);


template <typename OVEC, typename XVEC, typename IVEC,
    typename PVEC, typename IT2>
extern void csc_to_dense_transposed_c_vec(
    XVEC const & x, 
    IVEC const & i, 
    PVEC const & p, 
    IT2 const & nrow, IT2 const & ncol, 
    OVEC & out,
    int const & threads);


// return total number of rows.
template <typename XITER, typename IITER, typename PITER, 
    typename IT2,
    typename XITER2, typename IITER2, typename PITER2>
extern IT2 csc_rbind(
    std::vector<XITER> const & xvecs, 
    std::vector<IITER> const & ivecs, 
    std::vector<PITER> const & pvecs, 
    std::vector<IT2> const & nrows, 
    std::vector<IT2> const & ncols, 
    XITER2 ox,
    IITER2 oi,
    PITER2 op,
    int const & threads);


// return total number of rows.
template <template <typename...> class VEC, 
    typename XVEC, typename IVEC, typename PVEC, 
    typename COUNTVEC,
    typename XVEC2, typename IVEC2, typename PVEC2>
extern int csc_rbind_vec(
    VEC<XVEC> const & xvecs, 
    VEC<IVEC> const & ivecs, 
    VEC<PVEC> const & pvecs, 
    COUNTVEC const & nrows, 
    COUNTVEC const & ncols, 
    XVEC2 & ox,
    IVEC2 & oi,
    PVEC2 & op,
    int const & threads);

// return total number of columns.
template <typename XITER, typename IITER, typename PITER, typename IT2,
    typename XITER2, typename IITER2, typename PITER2>
extern IT2 csc_cbind(
    std::vector<XITER> const & xvecs, 
    std::vector<IITER> const & ivecs, 
    std::vector<PITER> const & pvecs, 
    std::vector<IT2> const & nrows, 
    std::vector<IT2> const & ncols, 
    XITER2 ox,
    IITER2 oi,
    PITER2 op,
    int const & threads);

// return total number of columns.
// XVEC and others only need to have [] operator, so will need the nrow and ncol opertors.
template <template < typename... > class VEC, typename XVEC, typename IVEC, typename PVEC, typename COUNTVEC,
    typename XVEC2, typename IVEC2, typename PVEC2>
extern int csc_cbind_vec(
    VEC<XVEC> const & xvecs, 
    VEC<IVEC> const & ivecs, 
    VEC<PVEC> const & pvecs, 
    COUNTVEC const & nrows, 
    COUNTVEC const & ncols, 
    XVEC2 & ox,
    IVEC2 & oi,
    PVEC2 & op,
    int const & threads);

// csc
template <typename OITER, typename XITER, typename PITER, typename IT2>
extern void csc_colsums_iter(
    XITER x, 
    PITER p, 
    IT2 const & ncol, 
    OITER out,
    int const & threads);

template <typename OVEC, typename XVEC, typename PVEC, typename IT2>
extern void csc_colsums_vec(
    XVEC const & x, 
    PVEC const & p, 
    IT2 const & ncol, 
    OVEC & out,
    int const & threads);


// csc
template <typename OITER, typename XITER, typename IITER, typename IT2, typename IT>
extern void csc_rowsums_iter(
    XITER x, 
    IITER i, 
    IT2 const & nrow, IT const & nzcount, 
    OITER out,
    int const & threads);

    
// csc
template <typename OVEC, typename XVEC, typename IVEC, typename IT2>
extern void csc_rowsums_vec(
    XVEC const & x, 
    IVEC const & i, 
    IT2 const & nrow, 
    OVEC & out,
    int const & threads);