#pragma once

/*
 *  copied and modified from wave/include/transform/wilcox.hpp
 *
 *  Created on: Jun 8, 2021
 *  Author: Tony Pan
 *  Affiliation: Institute for Data Engineering and Science
 *  			Georgia Institute of Technology, Atlanta, GA 30332
 */

// direct write to matrix may not be fast for cpp11:  proxy object creation and iterator creation....

#include <utility>
#include <unordered_map>
#include <vector>


template <typename IT, typename LABEL>
void sum_rank(std::vector<std::pair<IT, LABEL>> const & temp, 
  std::unordered_map<LABEL, size_t> const & z_cl_counts,
  size_t const & count, IT const & zero_val,
  std::unordered_map<LABEL, size_t> & rank_sums,
  double & tie_sum
);

double ranksum_to_pval(double const & val, size_t const & n1, 
  double const & tie_sum, size_t const & count, 
  int const & test_type, bool const & continuity);

// ids points to start of the column's positive element row ids
// x  points to the start fo the columns positive element values
// count is hte number of positive elements in the column.
template <typename IT_ITER, typename IDX_ITER, typename LABEL_ITER,
  typename IT,
  typename LABEL>
void spmat_sort(IT_ITER in, IDX_ITER ids,
  size_t const & nz_count, 
  LABEL_ITER labels, size_t const & count, IT const & zero_val,
  std::vector<std::pair<LABEL, size_t> > const & cl_counts,
  std::vector<std::pair<IT, LABEL> > & temp,
  std::unordered_map<LABEL, size_t> & z_cl_counts) ;

// rank_sums output:  map cluster to rank_sum.
// counts zeros and skip sorting of those.  this is faster than sorting all.
template <typename IT_ITER, typename LABEL_ITER,
  typename IT,
  typename LABEL>
void pseudosparse_sort(
  IT_ITER in, LABEL_ITER labels, size_t const & count, IT const & zero_val,
  std::vector<std::pair<LABEL, size_t> > const & cl_counts,
  std::vector<std::pair<IT, LABEL> > & temp,
  std::unordered_map<LABEL, size_t> & z_cl_counts);

// // rank_sums output:  map cluster to rank_sum.
// // old, and slower than pseudosparse.  do not use.
// template <typename IT_ITER, typename LABEL_ITER,
//   typename IT = typename std::iterator_traits<IT_ITER>::value_type,
//   typename LABEL = typename std::iterator_traits<LABEL_ITER>::value_type>
// void dense_wmw_summary(
//   IT_ITER in, LABEL_ITER labels, size_t const & count,
//   std::unordered_map<LABEL, size_t> & rank_sums, double & tie_sum);

// types:  0 = less, 1 = greater, 2 = twosided (default), 3 = U2
template <typename LABEL, typename OT_ITER>
void wmw(
  std::vector<std::pair<LABEL, size_t> > const & cl_counts, 
  std::unordered_map<LABEL, size_t> const & rank_sums, 
  double const & tie_sum,
  size_t const & count,
  OT_ITER out,
  int const & test_type, bool const & continuity) ;


template < typename XIT, typename LIT, typename LABEL>
void omp_dense_wmw(
    XIT mat, size_t const & nsamples, size_t const & nfeatures,
    LIT lab, 
    int rtype, 
    bool continuity_correction, 
    std::vector<double> &pv,
    std::vector<std::pair<LABEL, size_t> > &sorted_cluster_counts,
    int threads);


// =================================
template <typename XIT, typename IIT, typename PIT, 
    typename LIT, typename LABEL>
void omp_sparse_wmw(
    XIT x, IIT i, PIT p, size_t nsamples, size_t nfeatures,
    LIT lab,
    int rtype, 
    bool continuity_correction, 
    std::vector<double> &pv,
    std::vector<std::pair<LABEL, size_t> > &sorted_cluster_counts,
    int threads);


//  VEC version -------------------------



// ids points to start of the column's positive element row ids
// x  points to the start fo the columns positive element values
// count is hte number of positive elements in the column.
template <typename XVEC, typename IVEC, typename LVEC,
  typename XT,
  typename LABEL>
void spmat_sort_vec(XVEC const & in, IVEC const & ids,
  size_t const & nz_offset, 
  size_t const & nz_count, 
  LVEC const & labels, size_t const & count, XT const & zero_val,
  std::vector<std::pair<LABEL, size_t> > const & cl_counts,
  std::vector<std::pair<XT, LABEL> > & temp,
  std::unordered_map<LABEL, size_t> & z_cl_counts) ;


// rank_sums output:  map cluster to rank_sum.
// counts zeros and skip sorting of those.  this is faster than sorting all.
template <typename XVEC, typename LVEC,
  typename XT,
  typename LABEL>
void pseudosparse_sort_vec(
  XVEC const & in, size_t const & offset, 
  LVEC const & labels, size_t const & count, XT const & zero_val,
  std::vector<std::pair<LABEL, size_t> > const & cl_counts,
  std::vector<std::pair<XT, LABEL> > & temp,
  std::unordered_map<LABEL, size_t> & z_cl_counts) ;


template <typename LABEL, typename PVVEC>
void wmw_vec(
  std::vector<std::pair<LABEL, size_t> > const & cl_counts, 
  std::unordered_map<LABEL, size_t> const & rank_sums, 
  double const & tie_sum,
  size_t const & count,
  PVVEC & out, size_t const & offset,
  int const & test_type, bool const & continuity) ;


template < typename XVEC, typename LVEC, typename PVVEC, typename LABEL>
void vecsc_wmw_vecsc(
    XVEC const & mat, size_t const & nsamples, size_t const & nfeatures,
    LVEC const & lab, 
    int const & rtype, 
    bool const & continuity_correction, 
    PVVEC & pv,
    std::vector<std::pair<LABEL, size_t> > const & sorted_cluster_counts,
    int const & threads);

template < typename XMAT, typename LVEC, typename PVVEC, typename LABEL>
void matc_wmw_vecsc(
    XMAT const & mat, size_t const & nsamples, size_t const & nfeatures,
    LVEC const & lab, 
    int const & rtype, 
    bool const & continuity_correction, 
    PVVEC & pv,
    std::vector<std::pair<LABEL, size_t> > const & sorted_cluster_counts,
    int const & threads);


// // assume mat has [] operator to get a slice, then [] on slice to get element.
// template < typename XMAT, typename LVEC, typename PVMAT, typename LABEL = decltype(std::declval<LVEC>()[0])>
// void matc_wmw_matc(
//     XMAT const & mat, size_t const & nsamples, size_t const & nfeatures,
//     LVEC const & lab, 
//     int const & rtype, 
//     bool const & continuity_correction, 
//     PVMAT & pv,
//     std::vector<std::pair<LABEL, size_t> > const & sorted_cluster_counts,
//     int const & threads);


// =================================
template <typename XVEC, typename IVEC, typename PVEC, 
    typename LVEC, typename PVVEC, typename LABEL>
void csc_wmw_vecsc(
    XVEC const & x, IVEC const & i, PVEC const & p, size_t nsamples, size_t nfeatures,
    LVEC const & lab,
    int const & rtype, 
    bool const & continuity_correction, 
    PVVEC & pv,
    std::vector<std::pair<LABEL, size_t> > const & sorted_cluster_counts,
    int const & threads);

// template <typename XVEC, typename IVEC, typename PVEC, 
//     typename LVEC, typename PVMAT, typename LABEL = decltype(std::declval<LVEC>()[0])>
// void csc_wmw_matc(
//     XVEC const & x, IVEC const & i, PVEC const & p, size_t nsamples, size_t nfeatures,
//     LVEC const & lab,
//     int const & rtype, 
//     bool const & continuity_correction, 
//     PVMAT & pv,
//     std::vector<std::pair<LABEL, size_t> > const & sorted_cluster_counts,
//     int const & threads);
