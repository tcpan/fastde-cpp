#pragma once

/*
 *  copied and modified from wave/include/transform/wilcox.hpp
 *
 *  Created on: Jun 8, 2021
 *  Author: Tony Pan
 *  Affiliation: Institute for Data Engineering and Science
 *  			Georgia Institute of Technology, Atlanta, GA 30332
 */

#include <utility>
#include <unordered_map>
#include <vector>


template <typename IT, typename LABEL>
void sparse_sum_rank(std::pair<IT, LABEL> * temp, size_t const & nzc, 
  size_t const & count, IT const & zero_val,
  std::unordered_map<LABEL, size_t> const & z_cl_counts,
  std::unordered_map<LABEL, size_t> & rank_sums,
  double & tie_sum
);

// ids points to start of the column's positive element row ids
// x  points to the start fo the columns positive element values
// count is hte number of positive elements in the column.
template <typename IT_ITER, typename IDX_ITER, typename LABEL_ITER,
  typename IT = typename std::iterator_traits<IT_ITER>::value_type,
  typename IDX = typename std::iterator_traits<IDX_ITER>::value_type,
  typename LABEL = typename std::iterator_traits<LABEL_ITER>::value_type>
void sparse_wmw_summary(IT_ITER in, IDX_ITER ids,
  size_t const & nz_count, 
  LABEL_ITER labels, size_t const & count, IT const & zero_val,
  std::vector<std::pair<LABEL, size_t> > const & cl_counts,
  std::unordered_map<LABEL, size_t> & rank_sums, double & tie_sum) ;

// rank_sums output:  map cluster to rank_sum.
// counts zeros and skip sorting of those.  this is faster than sorting all.
template <typename IT_ITER, typename LABEL_ITER,
  typename IT = typename std::iterator_traits<IT_ITER>::value_type,
  typename LABEL = typename std::iterator_traits<LABEL_ITER>::value_type>
void pseudosparse_wmw_summary(
  IT_ITER in, LABEL_ITER labels, size_t const & count, IT const & zero_val,
  std::vector<std::pair<LABEL, size_t> > const & cl_counts,
  std::unordered_map<LABEL, size_t> & rank_sums, double & tie_sum);

// // rank_sums output:  map cluster to rank_sum.
// // old, and slower than pseudosparse.  do not use.
// template <typename IT_ITER, typename LABEL_ITER,
//   typename IT = typename std::iterator_traits<IT_ITER>::value_type,
//   typename LABEL = typename std::iterator_traits<LABEL_ITER>::value_type>
// void dense_wmw_summary(
//   IT_ITER in, LABEL_ITER labels, size_t const & count,
//   std::unordered_map<LABEL, size_t> & rank_sums, double & tie_sum);

// types:  0 = less, 1 = greater, 2 = twosided (default), 3 = U2
template <typename LABEL, typename OT_ITER,
  typename OT = typename std::iterator_traits<OT_ITER>::value_type>
void wmw(
  std::vector<std::pair<LABEL, size_t> > const & cl_counts, 
  std::unordered_map<LABEL, size_t> const & rank_sums, 
  double const & tie_sum,
  size_t const & count,
  OT_ITER out,
  int const & test_type, bool const & continuity) ;


template < typename XIT, typename LIT>
void omp_dense_wmw(
    XIT mat, size_t const & nsamples, size_t const & nfeatures,
    LIT lab, 
    int rtype, 
    bool continuity_correction, 
    std::vector<double> &pv,
    std::vector<std::pair<int, size_t> > &sorted_cluster_counts,
    int threads);


// =================================
template <typename XIT, typename IIT, typename PIT, 
    typename LIT>
void omp_sparse_wmw(
    XIT x, IIT i, PIT p, size_t nsamples, size_t nfeatures,
    LIT lab,
    int rtype, 
    bool continuity_correction, 
    std::vector<double> &pv,
    std::vector<std::pair<int, size_t> > &sorted_cluster_counts,
    int threads);


//  VEC version -------------------------



// ids points to start of the column's positive element row ids
// x  points to the start fo the columns positive element values
// count is hte number of positive elements in the column.
template <typename XVEC, typename IVEC, typename LVEC,
  typename XT = decltype(std::declval<XVEC>()[0]),
  typename LT = decltype(std::declval<LVEC>()[0])>
void sparse_wmw_summary_vec(XVEC const & in, IVEC const & ids,
  size_t const & nz_offset, 
  size_t const & nz_count, 
  LVEC const & labels, size_t const & count, XT const & zero_val,
  std::vector<std::pair<LT, size_t> > const & cl_counts,
  std::unordered_map<LT, size_t> & rank_sums, double & tie_sum) ;

// rank_sums output:  map cluster to rank_sum.
// counts zeros and skip sorting of those.  this is faster than sorting all.
template <typename XVEC, typename LVEC,
  typename XT = decltype(std::declval<XVEC>()[0]),
  typename LT = decltype(std::declval<LVEC>()[0])>
void pseudosparse_wmw_summary_vec(
  XVEC const & in, size_t const & offset, 
  LVEC const & labels, size_t const & count, XT const & zero_val,
  std::vector<std::pair<LT, size_t> > const & cl_counts,
  std::unordered_map<LT, size_t> & rank_sums, double & tie_sum);


template <typename LT, typename PVVEC,
  typename OT = decltype(std::declval<PVVEC>()[0])>
void wmw_vec(
  std::vector<std::pair<LT, size_t> > const & cl_counts, 
  std::unordered_map<LT, size_t> const & rank_sums, 
  double const & tie_sum,
  size_t const & count,
  PVVEC & out, size_t const & offset, size_t const & label_count,
  int const & test_type, bool const & continuity) ;


template < typename XVEC, typename LVEC, typename PVVEC, typename LT = decltype(std::declval<LVEC>()[0])>
void csc_dense_wmw_vec(
    XVEC const & mat, size_t const & nsamples, size_t const & nfeatures,
    LVEC const & lab, 
    int const & rtype, 
    bool const & continuity_correction, 
    PVVEC & pv,
    std::vector<std::pair<LT, size_t> > &sorted_cluster_counts,
    int const & threads);


// =================================
template <typename XVEC, typename IVEC, typename PVEC, 
    typename LVEC, typename PVVEC, typename LT = decltype(std::declval<LVEC>()[0])>
void csc_sparse_wmw_vec(
    XVEC const & x, IVEC const & i, PVEC const & p, size_t nsamples, size_t nfeatures,
    LVEC const & lab,
    int const & rtype, 
    bool const & continuity_correction, 
    PVVEC & pv,
    std::vector<std::pair<LT, size_t> > &sorted_cluster_counts,
    int const & threads);
