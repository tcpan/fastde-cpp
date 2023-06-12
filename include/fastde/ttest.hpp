#pragma once

// NOTE: column major R data structure.  getting a whole column is faster?
// INPUT: samples are in rows, all samples (in rows) for a feature (a column) are contiguous.
// OUTPUT: labels in rows, all label (in rows) for a feature (a column) are contiguous

/*
 *  Created on: Jan 31, 2022
 *  Author: Tony Pan
 *  Affiliation: Institute for Data Engineering and Science
 *  			Georgia Institute of Technology, Atlanta, GA 30332
 */

#include <stddef.h>

#include <utility>
#include <unordered_map>  // unordered_map
#include <vector>

#include "fastde/cluster_utils.hpp"

// #ifndef INV_SQRT2
// #define INV_SQRT2 0.70710678118
// #endif

// sum, sum of square, and 
template <typename IT>
struct gaussian_stats {
  size_t count;
  size_t zeros;
  IT sum;
  IT sum_of_square;
};

// ids points to start of the column's positive element row ids
// x  points to the start fo the columns positive element values
// count is hte number of positive elements in the column.
template <typename ITER, typename IDX_ITER, typename LABEL_ITER,
  typename IT = typename std::iterator_traits<ITER>::value_type,
  typename IDX = typename std::iterator_traits<IDX_ITER>::value_type,
  typename LABEL = typename std::iterator_traits<LABEL_ITER>::value_type>
void sparse_ttest_summary(ITER in, IDX_ITER ids, size_t const & nz_count, 
  LABEL_ITER labels, size_t const & count, IT const & zero_val,
  std::vector<std::pair<LABEL, size_t> > const & clust_counts,
  std::unordered_map<LABEL, gaussian_stats<double> > & gaussian_sums);

// compute count, mean and variance for each subclass.
template <typename ITER, typename LABEL_ITER,
  typename IT = typename std::iterator_traits<ITER>::value_type,
  typename LABEL = typename std::iterator_traits<LABEL_ITER>::value_type>
void dense_ttest_summary(
  ITER in, LABEL_ITER labels, size_t const & count, IT const & zero_val,
  std::vector<std::pair<LABEL, size_t> > const & clust_counts,
  std::unordered_map<LABEL, gaussian_stats<double> > & gaussian_sums);

template <typename T, typename DEGREE, bool precompute>
class incomplete_beta;

template <typename T, typename DEGREE>
class incomplete_beta<T, DEGREE, false> {

  protected:
    static constexpr size_t ITERS = 200;   /// seems to be okay with the starting a, b.

  public:

    // regularized incomplete beta function.  from https://codeplea.com/incomplete-beta-function-c
    // call with a, b, and x as originally supplied by the user.  use "flip" to reverse them.
    double incbeta(DEGREE const & a, DEGREE const & b, T const & x, bool const & flip) ;

    // note that dof/2.0 for a AND b is not what wikipedia has.
    incomplete_beta() {};

    int get_convergence_type(DEGREE const & _a, DEGREE const & _b, T const & x) ;

    T unchecked_incbeta(DEGREE const & _a, DEGREE const & _b, T const & x) ;

    T operator()(DEGREE const & _a, DEGREE const & _b, T const & x) ;
};


template <typename T, typename DEGREE>
class incomplete_beta<T, DEGREE, true> {

  protected:
    const double a;
    const double b;
    const double x_thresh;
    const double lbeta;  // symmetric
    std::vector<double> numerator_ab;
    std::vector<double> numerator_ba;
    // static constexpr double STOP = 1.0e-8;
    // static constexpr double TINY = 1.0e-30;
    static constexpr size_t ITERS = 200;

    void precompute();

  public:
    double incbeta(T const & x, bool const & flip ) ;

    // note that dof/2.0 for a AND b is not what wikipedia has.
    incomplete_beta(DEGREE const & _a, DEGREE const & _b) ;

    int get_convergence_type(T const & x) ;

    T unchecked_incbeta(T const & x, bool const & flip) ;

    T operator()(T const & x);

};

template <typename T, typename DEGREE, bool precompute>
class t_distribution;

// based on nist's description of continued fractions.  http://dlmf.nist.gov/8.17#SS5.p1
// following Lentz' Algorihtm
// specialized for double for now.
template <typename T, typename DEGREE>
class t_distribution<T, DEGREE, false> {
  protected:
    incomplete_beta<T, DEGREE, false> incbeta;

  public:
    // note that dof/2.0 for a AND b is not what wikipedia has.
    t_distribution() {};

    // directly compute pval (extreme value, or error probability)
    T pval(T const & t, DEGREE const & _v, int const & test_type = PVAL_GREATER) ;


    T cdf(T const & t, DEGREE const & _v) ;
};

template <typename T, typename DEGREE>
class t_distribution<T, DEGREE, true> {
  protected:
    DEGREE _v;
    incomplete_beta<T, DEGREE, true> incbeta;

  public:
    // note that dof/2.0 for a AND b is not what wikipedia has.
    t_distribution(DEGREE const & dof);

    // directly compute pval to avoid repeated 1-x that causes epsilon loss of precision.
    T pval(T const & t, int const & test_type = PVAL_GREATER) ;

    // compute CDF directly
    T cdf(T const & t) ;

};

// types:  
template <typename LABEL, typename OT_ITER,
  typename OT = typename std::iterator_traits<OT_ITER>::value_type>
void two_sample_ttest(
  std::unordered_map<LABEL, gaussian_stats<double> > const & gaussian_sums,
  std::vector<std::pair<LABEL, size_t> > const & clust_counts,
  OT_ITER out, 
  int const & test_type = PVAL_TWO_SIDED, 
  bool const equal_variance = false) ;

// omp t-Test for dense matrix
template < typename XIT, typename LIT>
void omp_dense_ttest(
    XIT mat, size_t const & nsamples, size_t const & nfeatures,
    LIT lab, 
    int alternative, 
    bool var_equal, 
    std::vector<double> &pv,
    std::vector<std::pair<int, size_t> > &sorted_cluster_counts,
    int threads) ;

// =================================
template <typename XIT, typename IIT, typename PIT, 
    typename LIT>
void omp_sparse_ttest(
    XIT x, IIT i, PIT p, size_t nsamples, size_t nfeatures,
    LIT lab,
    int alternative, 
    bool var_equal, 
    std::vector<double> &pv,
    std::vector<std::pair<int, size_t> > &sorted_cluster_counts,
    int threads);
