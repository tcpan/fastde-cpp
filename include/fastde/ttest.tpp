#pragma once

#include "fastde/ttest.hpp"

// NOTE: column major R data structure.  getting a whole column is faster?
// INPUT: samples are in rows, all samples (in rows) for a feature (a column) are contiguous.
// OUTPUT: labels in rows, all label (in rows) for a feature (a column) are contiguous

/*
 *  Created on: Jan 31, 2022
 *  Author: Tony Pan
 *  Affiliation: Institute for Data Engineering and Science
 *  			Georgia Institute of Technology, Atlanta, GA 30332
 */

#include <cmath>  // sqrt
#include <algorithm>
#include <type_traits>
#include <cassert>

#include <omp.h>

// #ifndef INV_SQRT2
// #define INV_SQRT2 0.70710678118
// #endif

#include "fastde/cluster_utils.hpp"
#include "fastde/benchmark_utils.hpp"



// ids points to start of the column's positive element row ids
// x  points to the start fo the columns positive element values
// count is hte number of positive elements in the column.
template <typename ITER, typename IDX_ITER, typename LABEL_ITER,
  typename IT,
  typename IDX,
  typename LABEL>
void sparse_ttest_summary(ITER in, IDX_ITER ids, size_t const & nz_count, 
  LABEL_ITER labels, size_t const & count, IT const & zero_val,
  std::vector<std::pair<LABEL, size_t> > const & clust_counts,
  std::unordered_map<LABEL, gaussian_stats<double> > & gaussian_sums) {


  if (count == 0) return;

  // ======== initialize summary
  // first walk through labels, treat as if all are zeros.
  for (auto item : clust_counts) {
    gaussian_sums[item.first] = {0, item.second, 0.0, 0.0}; 
  }

  // =====  scan and create summary.
  LABEL l;
  // LABEL_ITER lab_iter = labels;
  // then walk through the input, and adjust the non-zeros.
  IT val;
  IDX id;
  ITER in_iter = in;
  IDX_ITER ids_iter = ids;
  for (size_t i = 0; i < nz_count; ++i, ++in_iter, ++ids_iter) {
    val = *in_iter;
    id = *ids_iter;
    l = *(labels + id);
    // check if label present in gaussian_sums
    if (val != zero_val) {
      ++gaussian_sums[l].count;
      gaussian_sums[l].sum += val;
      gaussian_sums[l].sum_of_square += val * val;
    }
  }

  // then walk through label. and adjust for the zeros.
  for (auto it = gaussian_sums.begin(); it != gaussian_sums.end(); ++it) {
    it->second.zeros -= it->second.count;
    if (zero_val != 0.0) {
      it->second.sum += it->second.zeros * zero_val;
      it->second.sum_of_square += it->second.zeros * zero_val * zero_val;
    }

    // Rprintf("label %d, count %d zeros %d, sums %f, sums squared %f\n", it->first, it->second.count, it->second.zeros, it->second.sum, it->second.sum_of_square);
  }
}


// compute count, mean and variance for each subclass.
template <typename ITER, typename LABEL_ITER,
  typename IT,
  typename LABEL>
void dense_ttest_summary(
  ITER in, LABEL_ITER labels, size_t const & count, IT const & zero_val,
  std::vector<std::pair<LABEL, size_t> > const & clust_counts,
  std::unordered_map<LABEL, gaussian_stats<double> > & gaussian_sums) {

  if (count == 0) return;

  // ======== initialize summary
  for (auto item : clust_counts) {
    gaussian_sums[item.first] = {0, item.second, 0.0, 0.0}; 
  }

  // =====  scan and create summary.

   // forward walk
  IT val;
  LABEL l;
  ITER in_iter = in;
  LABEL_ITER lab_iter = labels;
  for (size_t i = 0; i < count; ++i, ++in_iter, ++lab_iter) {
    val = *in_iter;
    l = *lab_iter;
    // check if label present in gaussian_sums
    if (val != zero_val) {
      ++gaussian_sums[l].count;
      gaussian_sums[l].sum += val;
      gaussian_sums[l].sum_of_square += val * val;
    }
  }

  // zeros correction.
  for (auto it = gaussian_sums.begin(); it != gaussian_sums.end(); ++it) {
    it->second.zeros -= it->second.count;
    if (zero_val != 0.0) {
      it->second.sum += it->second.zeros * zero_val;
      it->second.sum_of_square += it->second.zeros * zero_val * zero_val;
    }

    // Rprintf("label %d, count %d zeros %d, sums %f, sums squared %f\n", it->first, it->second.count, it->second.zeros, it->second.sum, it->second.sum_of_square);
  }
}


// TODO:
// [x]  test that t val is same as R's.  yes.
// [x]  output not identical from tstat to pval.  could it be because of the (v/2, 1/2) vs (v/2, v/2) used in the codeplea method?
//      fixed: 1-x reaching EPS limit, or 1-x with x at EPS limit.  fixed by reorg code to avoid 1-(1-x) = 0.
// [ ]  erfc-like function for t-distribution.


// alternative implementation from https://stats.stackexchange.com/questions/394978/how-to-approximate-the-student-t-cdf-at-a-point-without-the-hypergeometric-funct
// goal is to have double precision, but that imple seems to have some fortran constants and may be iterating more??


//F(t) = 1 - I_x(t) (dof/2, 1/2), x(t) = dof/ (t^2+dof), t > 0, dof is degree of freedom 
// I_x(t) (a, b) is the regularized incomplete beta function, 
// I_x(t) (a, b) = B(x; a, b)/ B(a, b), B(x; a,b) is the incomplete beta function, and B(a, b) is the beta function
// B(a,b) = \Gamma(a) * \Gamma(b) / \Gamma(a+b)
// <cmath> provides tgamma() and lgamma()
//  https://codeplea.com/incomplete-beta-function-c  for computing the regularized incomplete beta function and the cdf.
// the class below leverages the ideas and code from https://codeplea.com/incomplete-beta-function-c
// to produce the regularized incomplete beta function output given a degree of freedom (v), and an x value derived from t and v.


// regularized incomplete beta function.  from https://codeplea.com/incomplete-beta-function-c
// call with a, b, and x as originally supplied by the user.  use "flip" to reverse them.
template <typename T, typename DEGREE>
double incomplete_beta<T, DEGREE, false>::incbeta(DEGREE const & a, DEGREE const & b, 
      T const & x, bool const & flip) {

      /*Find the first part before the continued fraction.*/
      const double _lbeta_ab = lgamma(a)+lgamma(b)-lgamma(a + b);  // symmetric
      double front = exp(log(static_cast<double>(x))*a+log(1.0 - static_cast<double>(x))*b-_lbeta_ab);  // symmetric

      double _a = flip ? static_cast<double>(b) : static_cast<double>(a);
      double _b = flip ? static_cast<double>(a) : static_cast<double>(b);
      // if flip, then we need 1-x, but this can be delayed to computing numerator.

      front /= _a;   // divide by a, or b, depending on "flip"

      /*Use Lentz's algorithm to evaluate the continued fraction.*/
      //  a couple of changes vs wikipedia:  1.  check 1-cd instead of f_n - f_{n-1}.  2.  b_i = 1.0 always, so f_0 and c_0 are both 1.0
      double f = 1.0, c = 1.0, d = 0.0;

      size_t m;
      double dm, dm2, coeff, numerator, cd;
      for (size_t i = 0; i <= ITERS; ++i) {
          m = i >> 1;
          dm = static_cast<double>(m);
          dm2 = static_cast<double>(m << 1);

          if (i == 0) {
              numerator = 1.0; /*First numerator is 1.0.*/
          } else if ((i & 0x1) == 0) {
              coeff = (dm*(_b-dm))/((_a+dm2-1)*(_a+dm2)); // Even term.  need to multily by _x.
              numerator = flip ? (coeff - coeff * x) : (coeff * x);  // if flip, _x = (1-x).
          } else {
              coeff = -((_a+dm)*(_a+_b+dm))/((_a+dm2)*(_a+dm2+1)); /*Odd term.*/
              numerator = flip ? (coeff - coeff * x) : (coeff * x);
          }

          /*Do an iteration of Lentz's algorithm.*/
          d = 1.0 + numerator * d;
          if (fabs(d) == 0) d = std::numeric_limits<double>::min();   // to avoid division by 0.
          d = 1.0 / d;

          c = 1.0 + numerator / c;
          if (fabs(c) == 0) c = std::numeric_limits<double>::min();

          cd = c*d;
          f *= cd;

          /*Check for stop.*/
          if (fabs(1.0-cd) < std::numeric_limits<double>::epsilon()) {  // CD is within 1 epsilon of 1.0
              // Rprintf("stopping:  cd %.17g, front %.17g, f %.17g, numerator %.17g, iter %d\n", cd, front, f, numerator, i);
              // return front * (f-1.0);    worried about epsilon not precise enough
              return (front*f - front);
          }
      }

      return std::numeric_limits<double>::infinity(); /*Needed more loops, did not converge.*/
    }


template <typename T, typename DEGREE>
int incomplete_beta<T, DEGREE, false>::get_convergence_type(DEGREE const & _a, DEGREE const & _b, T const & x) {
      if (x < 0.0 || x > 1.0) return 0;

      /*The continued fraction converges nicely for x < (a+1)/(a+b+2)*/
      double _x_thresh = (static_cast<double>(_a) + 1.0) / (static_cast<double>(_a) + static_cast<double>(_b) + 2.0);
      if (static_cast<double>(x) >= _x_thresh) return -1;
      else return 1;
    }

template <typename T, typename DEGREE>
T incomplete_beta<T, DEGREE, false>::unchecked_incbeta(DEGREE const & _a, DEGREE const & _b, T const & x) {
      return incbeta(_a, _b, x, false);   // do as is.  no flip.
    }


template <typename T, typename DEGREE>
T incomplete_beta<T, DEGREE, false>::operator()(DEGREE const & _a, DEGREE const & _b, T const & x) {
      int type = get_convergence_type(_a, _b, x);

      if (type == 0) return std::numeric_limits<T>::infinity();
      else if (type > 0) return incbeta(_a, _b, x, false);   // forward.  no flip
      else return (1.0 - incbeta(_a, _b, x, true));  // inverse.  swap a and b, and complement x inside the function.
    }



template <typename T, typename DEGREE>
void incomplete_beta<T, DEGREE, true>::precompute() {
      numerator_ab.resize(ITERS + 1);
      numerator_ab[0] = 1.0;
      numerator_ba.resize(ITERS + 1);
      numerator_ba[0] = 1.0;
      size_t m;
      double dm, dm2;
      for (size_t i = 1; i <= ITERS; ++i) {
        m = i >> 1;
        dm = static_cast<double>(m);
        dm2 = static_cast<double>(m << 1);

        numerator_ab[i] = 
          (i & 0x1) ? 
            -((a+dm)*(a+b+dm))/((a+dm2)*(a+dm2+1))  // odd term
            : (dm*(b-dm))/((a+dm2-1)*(a+dm2));  // even term
        numerator_ba[i] = 
          (i & 0x1) ? 
            -((b+dm)*(b+a+dm))/((b+dm2)*(b+dm2+1))  // odd term
            : (dm*(a-dm))/((b+dm2-1)*(b+dm2));  // even term
      }
    }


template <typename T, typename DEGREE>
double incomplete_beta<T, DEGREE, true>::incbeta(T const & x, bool const & flip ) {

      double front = exp(log(static_cast<double>(x))*a+log(1.0-static_cast<double>(x))*b - lbeta);

      double _a = flip ? static_cast<double>(b) : static_cast<double>(a);
      // double _b = flip ? static_cast<double>(a) : static_cast<double>(b);
      std::vector<double> & numerators = flip ? numerator_ba : numerator_ab;

      /*Find the first part before the continued fraction.*/
      front /= _a;

      /*Use Lentz's algorithm to evaluate the continued fraction.*/
      double f = 1.0, c = 1.0, d = 0.0;

      size_t i;
      double numerator, cd;
      for ( i = 0; i <= ITERS; ++i) {

          if (i == 0) {
              numerator = 1.0; /*First numerator is 1.0.*/
          } else {
              numerator = flip ? (numerators[i] - numerators[i] * x) : (numerators[i] * x);
          }

          /*Do an iteration of Lentz's algorithm.*/
          d = 1.0 + numerator * d;
          if (fabs(d) == 0) d = std::numeric_limits<double>::min();   // to avoid division by 0.
          d = 1.0 / d;

          c = 1.0 + numerator / c;
          if (fabs(c) == 0) c = std::numeric_limits<double>::min();

          cd = c*d;
          f *= cd;

          /*Check for stop.*/
          if (fabs(1.0-cd) < std::numeric_limits<double>::epsilon()) {
              return (front * f - front);
          }
      }

      return std::numeric_limits<double>::infinity(); /*Needed more loops, did not converge.*/
    }

    // note that dof/2.0 for a AND b is not what wikipedia has.
template <typename T, typename DEGREE>
incomplete_beta<T, DEGREE, true>::incomplete_beta(DEGREE const & _a, DEGREE const & _b) :
      a(_a),
      b(_b), 
      x_thresh((a + 1.0) / (a + b + 2.0)),
      lbeta(lgamma(a)+lgamma(b)-lgamma(a+b)) {

      precompute();
    }



template <typename T, typename DEGREE>
int incomplete_beta<T, DEGREE, true>::get_convergence_type(T const & x) {
      if (x < 0.0 || x > 1.0) return 0;
      /*The continued fraction converges nicely for x < (a+1)/(a+b+2)*/
      if (static_cast<double>(x) >= x_thresh) return -1;
      else return 1;
    }


template <typename T, typename DEGREE>
T incomplete_beta<T, DEGREE, true>::unchecked_incbeta(T const & x, bool const & flip) {
      return incbeta(x, flip);   // flip is convergence_type == -1
    }


template <typename T, typename DEGREE>
T incomplete_beta<T, DEGREE, true>::operator()(T const & x) {
      int type = get_convergence_type(x);

      if (type == 0) return std::numeric_limits<T>::infinity();
      else if (type > 0) return incbeta(x, false);
      else return (1.0 - incbeta(x, true));  // note (1.0 - X) is subject epsilon limits for X close to 0 and close to 1.
    }



// based on nist's description of continued fractions.  http://dlmf.nist.gov/8.17#SS5.p1
// following Lentz' Algorihtm
// specialized for double for now.


    // directly compute pval (extreme value, or error probability)
template <typename T, typename DEGREE>
T t_distribution<T, DEGREE, false>::pval(T const & t, DEGREE const & _v, int const & test_type) {
      
      // first calcualte the cdf.
      // double _cdf;
      double o;
      if (_v == std::numeric_limits<T>::infinity()) {
        // _cdf = 0.5 * erfc(-t * M_SQRT1_2);

        // Rprintf("equal variance dof: %f \n", df);
        if (test_type == PVAL_TWO_SIDED)  // 2 sided
          // o = 2.0L * (1.0L - cdf( fabs(tstat)));   // this is correct.
          o = erfc(fabs(t) * M_SQRT1_2);   // the erfc part has range 1 to 2.
        else if (test_type == PVAL_LESS)  // less.
          o = 0.5 * erfc(-t * M_SQRT1_2); 
        else if (test_type == PVAL_GREATER)  // greater.  same as LESS for -t 
          // o = 1.0L - eqvar_t_cdf( tstat );
          o = 0.5 * erfc(t * M_SQRT1_2);
        else
          o = 1.0L;
      } else {
        double x1 = _v / (t * t + _v);
        int converge_flag = incbeta.get_convergence_type(_v * 0.5, 0.5, x1);

        if (converge_flag == 0) return std::numeric_limits<T>::infinity();

        if (converge_flag > 0) {
          // out1 =  1.0 - 0.5 * incbeta.incbeta(_v * 0.5, 0.5, x1, false);
          // Rprintf("equal variance dof: %f \n", df);
          double y = incbeta.incbeta(_v * 0.5, 0.5, x1, false);
          if (test_type == PVAL_TWO_SIDED)  // 2 sided
            // o = 2.0L * (1.0L - eqvar_t_cdf( fabs(tstat)));   // this is correct.
            o = y;
          else if (test_type == PVAL_LESS)  // less.
            // o = eqvar_t_cdf( tstat ); 
            o = 1.0 - 0.5 * y;
          else if (test_type == PVAL_GREATER)  // greater.  same as LESS for -t 
            // o = 1.0L - eqvar_t_cdf( tstat);
            o = 0.5 * y;
          else
            o = 1.0L;
        } else {
          // out1 =  0.5 + 0.5 * incbeta.incbeta(_v * 0.5, 0.5, x1, true);
          // Rprintf("equal variance dof: %f \n", df);
          double y = incbeta.incbeta(_v * 0.5, 0.5, x1, true);
          if (test_type == PVAL_TWO_SIDED)  // 2 sided
            // o = 2.0L * (1.0L - eqvar_t_cdf( fabs(tstat)));   // this is correct.
            o = 1.0 - y;
          else if (test_type == PVAL_LESS)  // less.
            // o = eqvar_t_cdf( tstat ); 
            o = 0.5 + 0.5 * y;
          else if (test_type == PVAL_GREATER)  // greater.  same as LESS for -t 
            // o = 1.0L - eqvar_t_cdf( tstat);
            o = 0.5 - 0.5 * y;
          else
            o = 1.0L;
        }
      }

      return o;
    }


template <typename T, typename DEGREE>
T t_distribution<T, DEGREE, false>::cdf(T const & t, DEGREE const & _v) {
      //  if _v is infinity, then use normal CDF, so erfc()?
              
      if (_v == std::numeric_limits<T>::infinity()) return 0.5 * erfc(-t * M_SQRT1_2);
      
      // this takes fewer iterations to converge (for large _v?), and seems to be numerically more stable (fewer inf)
      T out1; //, ib;
      double x1 = _v / (t * t + _v);

      int converge_flag = incbeta.get_convergence_type(_v * 0.5, 0.5, x1);

      if (converge_flag == 0) out1 = std::numeric_limits<T>::infinity();
      else if (converge_flag > 0) out1 =  1.0 - 0.5 * incbeta.incbeta(_v * 0.5, 0.5, x1, false);
      else out1 =  0.5 + 0.5 * incbeta.incbeta(_v * 0.5, 0.5, x1, true);

      // Rprintf("out %.17g. is this a double?\n", out1); 

      // this method here takes a lot more iterations to converge than the parameters above.
      // double out;
      // double sqttv = sqrt(t * t + _v);
      // double x = (t + sqttv) / (2.0 * sqttv);
      // out = incbeta(_v * 0.5, _v * 0.5, x);

      // Rprintf("cdf with :  t %.17g, dof  %.17g,  x1 %.17g,  x %.17g,  test %.17g, output %.17g \n", t, _v, x1, x, out1, out);
    
      return out1;
    }



    // note that dof/2.0 for a AND b is not what wikipedia has.
template <typename T, typename DEGREE>
t_distribution<T, DEGREE, true>::t_distribution(DEGREE const & dof) : _v(dof), incbeta(dof * 0.5, 0.5) {}

    // directly compute pval to avoid repeated 1-x that causes epsilon loss of precision.
template <typename T, typename DEGREE>
T t_distribution<T, DEGREE, true>::pval(T const & t, int const & test_type) {
      
      // first calcualte the cdf.
      // double _cdf;
      double o;
      if (_v == std::numeric_limits<T>::infinity()) {
        // _cdf = 0.5 * erfc(-t * M_SQRT1_2);

        // Rprintf("equal variance dof: %f \n", df);
        if (test_type == PVAL_TWO_SIDED)  // 2 sided
          // o = 2.0L * (1.0L - cdf( fabs(tstat)));   // this is correct.
          o = erfc(fabs(t) * M_SQRT1_2);   // the erfc part has range 1 to 2.
        else if (test_type == PVAL_LESS)  // less.
          o = 0.5 * erfc(-t * M_SQRT1_2); 
        else if (test_type == PVAL_GREATER)  // greater.  same as LESS for -t 
          // o = 1.0L - eqvar_t_cdf( tstat );
          o = 0.5 * erfc(t * M_SQRT1_2);
        else
          o = 1.0L;
      } else {
        double x1 = _v / (t * t + _v);
        int converge_flag = incbeta.get_convergence_type(x1);

        if (converge_flag == 0) return std::numeric_limits<T>::infinity();

        if (converge_flag > 0) {
          // out1 =  1.0 - 0.5 * incbeta.incbeta(x1, false);
          // Rprintf("equal variance dof: %f \n", df);
          double y = incbeta.incbeta(x1, false);
          if (test_type == PVAL_TWO_SIDED)  // 2 sided
            // o = 2.0L * (1.0L - eqvar_t_cdf( fabs(tstat)));   // this is correct.
            o = y;
          else if (test_type == PVAL_LESS)  // less.
            // o = eqvar_t_cdf( tstat ); 
            o = 1.0 - 0.5 * y;
          else if (test_type == PVAL_GREATER)  // greater.  same as LESS for -t 
            // o = 1.0L - eqvar_t_cdf( tstat);
            o = 0.5 * y;
          else
            o = 1.0L;
        } else {
          // out1 =  0.5 + 0.5 * incbeta.incbeta(x1, true);
          // Rprintf("equal variance dof: %f \n", df);
          double y = incbeta.incbeta(x1, true);
          if (test_type == PVAL_TWO_SIDED)  // 2 sided
            // o = 2.0L * (1.0L - eqvar_t_cdf( fabs(tstat)));   // this is correct.
            o = 1.0 - y;
          else if (test_type == PVAL_LESS)  // less.
            // o = eqvar_t_cdf( tstat ); 
            o = 0.5 + 0.5 * y;
          else if (test_type == PVAL_GREATER)  // greater.  same as LESS for -t 
            // o = 1.0L - eqvar_t_cdf( tstat);
            o = 0.5 - 0.5 * y;
          else
            o = 1.0L;
        }
      }
      
      return o;
    }

    // compute CDF directly
template <typename T, typename DEGREE>
T t_distribution<T, DEGREE, true>::cdf(T const & t) {
      //  if _v is infinity, then use normal CDF, so erfc()?
              
      if (_v == std::numeric_limits<T>::infinity()) return 0.5 * erfc(-t * M_SQRT1_2);
      
      // this takes fewer iterations to converge (for large _v?), and seems to be numerically more stable (fewer inf)
      T out1; //, ib;
      double x1 = _v / (t * t + _v);

      int converge_flag = incbeta.get_convergence_type(x1);

      if (converge_flag == 0) out1 = std::numeric_limits<T>::infinity();
      else if (converge_flag > 0) out1 =  1.0 - 0.5 * incbeta.incbeta(x1, false);
      else out1 =  0.5 + 0.5 * incbeta.incbeta(x1, true);

      // Rprintf(" out %.17g. is this a double?\n", out1); 

      // params below converges to same results but take more iterations.
      // double out;
      // double sqttv = sqrt(t * t + _v);
      // double x = (t + sqttv) / (2.0 * sqttv);
      // out = incbeta(_v * 0.5, _v * 0.5, x);

      // Rprintf("cdf with :  t %.17g, dof  %.17g,  x1 %.17g,  x %.17g,  test %.17g, output %.17g \n", t, _v, x1, x, out1, out);
    
      return out1;
    }


// types:  
template <typename LABEL, typename OT_ITER,
  typename OT>
void two_sample_ttest(
  std::unordered_map<LABEL, gaussian_stats<double> > const & gaussian_sums,
  std::vector<std::pair<LABEL, size_t> > const & clust_counts,
  OT_ITER out, 
  int const & test_type, 
  bool const equal_variance) {
  
  size_t count = 0;

  // sort the labels.
  double total_sum = 0, total_sum_squared = 0;
  for (auto item : gaussian_sums) {
    total_sum += item.second.sum;
    total_sum_squared += item.second.sum_of_square;
    count += (item.second.count + item.second.zeros);
  }

  if (count == 0) return;

  // iterate over the labels in sorted order.
  // compute 1 vs all else.
  double mean1, mean2;
  double sum1, sum2, sumsum2;
  double var1, var2;
  double inv_n1, inv_n2, inv_n1_1, inv_n2_1;
  size_t n1, n2;
  OT_ITER out_iter = out;
  double tstat;
  double df;
  double df_var1, df_var2;

  typename std::iterator_traits<OT_ITER>::value_type o;

  t_distribution<double, double, true> eqvar_t_cdf(count - 2);


  for (auto item : clust_counts) {
    LABEL l = item.first;
    auto stat = gaussian_sums.at(l);
    n1 = stat.count + stat.zeros;
    inv_n1 = 1.0 / static_cast<double>(n1);
    inv_n1_1 = 1.0 / static_cast<double>(n1 - 1);
    sum1 = stat.sum;
    mean1 = (n1 > 0) ? sum1 * inv_n1 : 0.0;
    var1 = (n1 > 1) ? (stat.sum_of_square - sum1 * sum1 * inv_n1) * inv_n1_1 : 0.0;
    // var1 = (n1 > 0) ? (stat.sum_of_square * inv_n1 - mean1 * mean1) : 0.0;

    n2 = count - n1;
    inv_n2 = 1.0 / static_cast<double>(n2);
    inv_n2_1 = 1.0 / static_cast<double>(n2 - 1);
    sum2 = total_sum - sum1;
    sumsum2 = total_sum_squared - stat.sum_of_square;
    mean2 = (n2 > 0) ? sum2 * inv_n2 : 0.0;
    var2 = (n2 > 1) ? (sumsum2 - sum2 * sum2 * inv_n2 ) * inv_n2_1 : 0.0;
    // var2 = (n2 > 0) ? (sumsum2 * inv_n2 - mean2 * mean2) : 0.0;

    // real computation:
    // if ((var1 > 2.0*var2) || (var2 > 2.0*var1)) {
    if (! equal_variance) {
      df_var1 = var1 * inv_n1;
      df_var2 = var2 * inv_n2;
      df = (df_var1 + df_var2) * (df_var1 + df_var2) / (df_var1 * df_var1 * inv_n1_1 + df_var2 * df_var2 * inv_n2_1);   // R's df?
      // different variance
      tstat = (mean1 - mean2) / std::sqrt(var1 * inv_n1 + var2 * inv_n2);
    } else {
      // similar, or equal variance.
      df = count - 2;
      if (n1 == n2) 
        tstat = (mean1 - mean2) / std::sqrt( ( var1 + var2 )  * inv_n1 );
      else
        tstat = (mean1 - mean2) / std::sqrt( ( static_cast<double>(n1-1) * var1 + static_cast<double>(n2-1) * var2 ) / static_cast<double>(n1 + n2 - 2) * ( inv_n1 + inv_n2 ));
    }
    if (test_type == TEST_VAL) {
      o = tstat;
    } else if (test_type == TEST_PARAM) {
      o = df;
    } else {
      // convert t statistics to pvalue.  approximation with Normal distribution for DOF >= 30 is not precise enough.  t-distribution table shows DOF=120 and empirical test with DOF = 1000 are not precise enough.  use exact compute instead.

      // doing this directly with incomplete beta instead of via t-distribution cdf to avoid repeated (1-x) reducing precision.

      // do it in the CDF class' pval funtion

      // if equal variance assumption, then we can reuse the cdf function.
      // else we'd need to set up new cdf for each dof.  Since dof in that case is dependent on variance of the 2 populations, unlikely to be able to reuse.
      if (equal_variance) {
        o = eqvar_t_cdf.pval(tstat, test_type);
      } else {
        // Rprintf("unequal variance label %d, count %d, dof: %f, tstat %f \n", l, item.second, df, tstat);
        t_distribution<double, double, false> cdf;
        o = cdf.pval(tstat, df, test_type);
      }
    }
    *out_iter = o;
    ++out_iter;
  }


}


// omp t-Test for dense matrix
template < typename XIT, typename LIT>
void omp_dense_ttest(
    XIT mat, size_t const & nsamples, size_t const & nfeatures,
    LIT lab, 
    int alternative, 
    bool var_equal, 
    std::vector<double> &pv,
    std::vector<std::pair<int, size_t> > &sorted_cluster_counts,
    int threads) {

  // std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double>> start = std::chrono::steady_clock::now();

  // get the number of unique labels.
  sorted_cluster_counts.clear();
  count_clusters(lab, nsamples, sorted_cluster_counts, threads);
  size_t label_count = sorted_cluster_counts.size();

  // ---- output pval matrix
  pv.clear(); pv.resize(nfeatures * label_count);

  // ------------------------ parameter
  // int alternative, threads;
  // bool as_dataframe, var_equal;
  // import_de_common_params(alternative, var_equal,
  //   alternative, var_equal);
  // import_r_common_params(as_dataframe, threads,
  //   as_dataframe, threads);

  // ------------------------ compute
  omp_set_num_threads(threads);
  // Rprintf("THREADING: using %d threads\n", threads);


#pragma omp parallel num_threads(threads)
  {
    int tid = omp_get_thread_num();
    int block = nfeatures / threads;
    int rem = nfeatures - threads * block;
    size_t offset = tid * block + (tid > rem ? rem : tid);
    int nid = tid + 1;
    size_t end = nid * block + (nid > rem ? rem : nid);

    std::unordered_map<int, gaussian_stats<double> > gaussian_sums;
    // printf("thread %d summarizing feature %ld\n", omp_get_thread_num(), offset);
    for(; offset < end; ++offset) {
      // directly compute matrix and res pointers.
      dense_ttest_summary(mat + offset * nsamples, 
        lab, nsamples, 0.0, sorted_cluster_counts, gaussian_sums);

      two_sample_ttest(gaussian_sums, sorted_cluster_counts, 
        pv.data() + offset * label_count, alternative, var_equal);
    }
  }
  // Rprintf("[TIME] T-test Elapsed(ms)= %f\n", since(start).count());

}




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
    int threads) {
  // Rprintf("here 1\n");

    std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double>> start;

  start = std::chrono::steady_clock::now();


  // get the number of unique labels.
  sorted_cluster_counts.clear();
  count_clusters(lab, nsamples, sorted_cluster_counts, threads);
  size_t label_count = sorted_cluster_counts.size();


  // ---- output pval matrix
  pv.clear(); pv.resize(nfeatures * label_count);

  // ------------------------ parameter
  // int alternative, threads;
  // bool as_dataframe, var_equal;
  // import_de_common_params(alternative, var_equal,
  //   alternative, var_equal);
  // import_r_common_params(as_dataframe, threads,
  //   as_dataframe, threads);

  // ------------------------- compute
  omp_set_num_threads(threads);
  // Rprintf("THREADING: using %d threads\n", threads);

  
#pragma omp parallel num_threads(threads)
  {
    int tid = omp_get_thread_num();
    int block = nfeatures / threads;
    int rem = nfeatures - threads * block;
    size_t offset = tid * block + (tid > rem ? rem : tid);
    int nid = tid + 1;
    size_t end = nid * block + (nid > rem ? rem : nid);
  
    long nz_offset, nz_count;
    std::unordered_map<int, gaussian_stats<double> > gaussian_sums;

    for(; offset < end; ++offset) {
      nz_offset = *(p + offset);
      nz_count = *(p + offset+1) - nz_offset;
      
      // directly compute matrix and res pointers.
      // Rprintf("thread %d processing feature %d, nonzeros %d\n", omp_get_thread_num(), offset, nz_count);
      sparse_ttest_summary(x +nz_offset, i + nz_offset, nz_count,
        lab, nsamples, 0.0, sorted_cluster_counts, gaussian_sums);

      two_sample_ttest(gaussian_sums, sorted_cluster_counts,
        pv.data() + offset * label_count, alternative, var_equal);
    }
  }
  // Rprintf("[TIME] T-test Elapsed(ms)= %f\n", since(start).count());

}

