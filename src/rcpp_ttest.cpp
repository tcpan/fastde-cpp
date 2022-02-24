#include <Rcpp.h>
#include <R.h>
#include <Rdefines.h>

// TRY WITH Rcpp classes (NumericMatrix, NumericVector), or use Rcpp::as<std::vector<double>>


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
#include <tuple>
#include <unordered_map>  // unordered_map
#include <unordered_set> 
#include <set>
#include <vector>
#include <algorithm>
#include <type_traits>

// Enable C++11
// [[Rcpp::plugins(cpp11)]]

// Enable OpenMP (excludes macOS)
// [[Rcpp::plugins(openmp)]]

// NOTE: do not use Rprintf inside thread.

#include <omp.h>

// #ifndef INV_SQRT2
// #define INV_SQRT2 0.70710678118
// #endif

#include "rcpp_cluster_utils.hpp"
#include "rcpp_data_utils.hpp"

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
  std::unordered_map<LABEL, gaussian_stats<double>> & gaussian_sums) {

  gaussian_sums.clear();

  if (count == 0) return;

  // =====  scan and create summary.

  // first walk through labels, treat as if all are zeros.
  LABEL l;
  LABEL_ITER lab_iter = labels;
  for (size_t i = 0; i < count; ++i, ++lab_iter) {
    l = *lab_iter;
    if (gaussian_sums.count(l) == 0 ) {
      // insert new one.
      gaussian_sums[l] = {0, 1, 0.0, 0.0};
    } else {
      ++gaussian_sums[l].zeros;
    }
  }

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
  }
}


// compute count, mean and variance for each subclass.
template <typename ITER, typename LABEL_ITER,
  typename IT = typename std::iterator_traits<ITER>::value_type,
  typename LABEL = typename std::iterator_traits<LABEL_ITER>::value_type>
void dense_ttest_summary(
  ITER in, LABEL_ITER labels, size_t const & count, IT const & zero_val,
  std::unordered_map<LABEL, gaussian_stats<double>> & gaussian_sums) {

  gaussian_sums.clear();

  if (count == 0) return;

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
    if (gaussian_sums.count(l) == 0 ) {
      // insert new one.
      if (val == zero_val)
        gaussian_sums[l] = {0, 1, 0.0, 0.0};
      else 
        gaussian_sums[l] = {1, 0, val, (val * val)};
    } else {
      if (val == zero_val)
        ++gaussian_sums[l].zeros;
      else {
        ++gaussian_sums[l].count;
        gaussian_sums[l].sum += val;
        gaussian_sums[l].sum_of_square += val * val;
      }
    }
  }

  if (zero_val != 0.0)
    for (auto it = gaussian_sums.begin(); it != gaussian_sums.end(); ++it) {
      it->second.sum += it->second.zeros * zero_val;
      it->second.sum_of_square += it->second.zeros * zero_val * zero_val;
    }

}

// TODO:
// [x]  test that t val is same as R's.  yes.
// [ ]  output not identical from tstat to pval.  could it be because of the (v/2, 1/2) vs (v/2, v/2) used in the codeplea method?
// [ ]  erfc-like function for t-distribution.

//F(t) = 1 - I_x(t) (dof/2, 1/2), x(t) = dof/ (t^2+dof), t > 0, dof is degree of freedom 
// I_x(t) (a, b) is the regularized incomplete beta function, 
// I_x(t) (a, b) = B(x; a, b)/ B(a, b), B(x; a,b) is the incomplete beta function, and B(a, b) is the beta function
// B(a,b) = \Gamma(a) * \Gamma(b) / \Gamma(a+b)
// <cmath> provides tgamma() and lgamma()
//  https://codeplea.com/incomplete-beta-function-c  for computing the regularized incomplete beta function and the cdf.
// the class below leverages the ideas and code from https://codeplea.com/incomplete-beta-function-c
// to produce the regularized incomplete beta function output given a degree of freedom (v), and an x value derived from t and v.
template <typename T, bool precompute>
class t_distribution_cdf;


template <typename T>
class t_distribution_cdf<T, false> {
  protected:
    static constexpr double STOP = 1.0e-8;
    static constexpr double TINY = 1.0e-30;
    static constexpr size_t ITERS = 200;

    double incbeta(double const & _a, double const & _b, double const & x) {
      if (x < 0.0 || x > 1.0) return INFINITY;

      /*The continued fraction converges nicely for x < (a+1)/(a+b+2)*/
      double _x_thresh = (_a + 1.0) / (_a + _b + 2.0);
      if (x >= _x_thresh) {
        return 1.0 - incbeta_ondemand(_b, _a, 1.0 - x);
      } else {
        return incbeta_ondemand(_a, _b, x); 
      }
    }

    // regularized incomplete beta function.  from https://codeplea.com/incomplete-beta-function-c
    double incbeta_ondemand(double const & _a, double const & _b, 
      double const & x) {

      /*Find the first part before the continued fraction.*/
      const double _lbeta_ab = lgamma(_a)+lgamma(_b)-lgamma(_a + _b);
      const double front = exp(log(x)*_a+log(1.0-x)*_b-_lbeta_ab) / _a;

      /*Use Lentz's algorithm to evaluate the continued fraction.*/
      double f = 1.0, c = 1.0, d = 0.0;

      int i, m;
      for (i = 0; i <= ITERS; ++i) {
          m = i/2;

          double numerator;
          if (i == 0) {
              numerator = 1.0; /*First numerator is 1.0.*/
          } else if (i % 2 == 0) {
              numerator = (m*(_b-m)*x)/((_a+2*m-1)*(_a+2*m)); /*Even term.*/
          } else {
              numerator = -((_a+m)*(_a+_b+m)*x)/((_a+2*m)*(_a+2*m+1)); /*Odd term.*/
          }

          /*Do an iteration of Lentz's algorithm.*/
          d = 1.0 + numerator * d;
          if (fabs(d) < TINY) d = TINY;
          d = 1.0 / d;

          c = 1.0 + numerator / c;
          if (fabs(c) < TINY) c = TINY;

          const double cd = c*d;
          f *= cd;

          /*Check for stop.*/
          if (fabs(1.0-cd) < STOP) {
              return front * (f-1.0);
          }
      }

      return INFINITY; /*Needed more loops, did not converge.*/
    }


  public:
    // note that dof/2.0 for a AND b is not what wikipedia has.
    t_distribution_cdf() {};

    T operator()(T const & t, double const & _v) {
      
      double x = _v / (t * t + _v);

      if (x < 0.0 || x > 1.0) return INFINITY;

      return 1.0 - 0.5 * incbeta(_v/2.0, 0.5, x);
    }

};
template <typename T>
class t_distribution_cdf<T, true> {
  protected:
    const double v;
    const double a;
    const double b;
    const double x_thresh;
    const double lbeta_ab;
    std::vector<double> numerator_ab;
    std::vector<double> numerator_ba;
    static constexpr double STOP = 1.0e-8;
    static constexpr double TINY = 1.0e-30;
    static constexpr size_t ITERS = 200;

    void precompute() {
      numerator_ab.resize(ITERS + 1);
      numerator_ab[0] = 1.0;
      numerator_ba.resize(ITERS + 1);
      numerator_ba[0] = 1.0;
      size_t m;
      for (size_t i = 1; i <= ITERS; ++i) {
        m = i >> 1;
        numerator_ab[i] = 
          (i & 0x1) ? 
            -((a+m)*(a+b+m))/((a+2*m)*(a+2*m+1))  // odd term
            : (m*(b-m))/((a+2*m-1)*(a+2*m));  // even term
        numerator_ba[i] = 
          (i & 0x1) ? 
            -((b+m)*(b+a+m))/((b+2*m)*(b+2*m+1))  // odd term
            : (m*(a-m))/((b+2*m-1)*(b+2*m));  // even term
      }
    }

    // regularized incomplete beta function.  from https://codeplea.com/incomplete-beta-function-c
    double incbeta(double const & x) {
      if (x < 0.0 || x > 1.0) return INFINITY;

      /*The continued fraction converges nicely for x < (a+1)/(a+b+2)*/
      if (x > x_thresh) {
        return 1.0 - incbeta_precomputed(b, a, 1.0 - x, numerator_ba);
      } else {
        return incbeta_precomputed(a, b, x, numerator_ab); 
      }
    }

    double incbeta_precomputed(double const & _a, double const & _b,
      double const & x,
      std::vector<double> const & numerator_front ) {

      /*Find the first part before the continued fraction.*/
      const double front = exp(log(x)*_a+log(1.0-x)*_b - lbeta_ab) / _a;

      /*Use Lentz's algorithm to evaluate the continued fraction.*/
      double f = 1.0, c = 1.0, d = 0.0;

      size_t i;
      for ( i = 0; i <= ITERS; ++i) {

          double numerator;
          if (i == 0) {
              numerator = 1.0; /*First numerator is 1.0.*/
          } else {
              numerator = numerator_front[i] * x;
          }

          /*Do an iteration of Lentz's algorithm.*/
          d = 1.0 + numerator * d;
          if (fabs(d) < TINY) d = TINY;
          d = 1.0 / d;

          c = 1.0 + numerator / c;
          if (fabs(c) < TINY) c = TINY;

          const double cd = c*d;
          f *= cd;

          /*Check for stop.*/
          if (fabs(1.0-cd) < STOP) {
              return front * (f-1.0);
          }
      }

      return INFINITY; /*Needed more loops, did not converge.*/
    }

  public:
    // note that dof/2.0 for a AND b is not what wikipedia has.
    t_distribution_cdf(double const & dof) : v(dof),
      a(dof / 2.0),
      b(0.5), 
      x_thresh((a + 1.0) / (a + b + 2.0)),
      lbeta_ab(lgamma(a)+lgamma(b)-lgamma(a+b)) {

      precompute();
    }

    T operator()(T const & t) {
      assert(precomputed && "class instance is not set up with precomputing");

      double x = v / (t * t + v);

      // double x = (t + sqrt(t * t + v)) / (2.0 * sqrt(t * t + v));
      if (x < 0.0 || x > 1.0) return INFINITY;

      return 1.0 - 0.5 * incbeta(x); // tdist = 1 - 0.5 incbeta(a, b, x)
    }
};


// types:  
template <typename LABEL, typename OT_ITER,
  typename OT = typename std::iterator_traits<OT_ITER>::value_type>
void two_sample_ttest(
  std::unordered_map<LABEL, gaussian_stats<double>> const & gaussian_sums,
  OT_ITER out, 
  int const & test_type = PVAL_TWO_SIDED, 
  bool const equal_variance = false) {
  
  size_t count = 0;

  // sort the labels.
  std::vector<LABEL> ordered_labels;
  double total_sum = 0, total_sum_squared = 0;
  for (auto item : gaussian_sums) {
    ordered_labels.emplace_back(item.first);
    total_sum += item.second.sum;
    total_sum_squared += item.second.sum_of_square;
    count += (item.second.count + item.second.zeros);
  }

  if (count == 0) return;

  std::sort(ordered_labels.begin(), ordered_labels.end());

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

  t_distribution_cdf<double, true> eqvar_t_cdf(count - 2);

  for (LABEL l : ordered_labels) {
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
    } else {
      // convert t statistics to pvalue.  approximation with Normal distribution for DOF >= 30 is not precise enough.  t-distribution table shows DOF=120 and empirical test with DOF = 1000 are not precise enough.  use exact compute instead.

      // if equal variance assumption, then we can reuse the cdf function.
      // else we'd need to set up new cdf for each dof.  Since dof in that case is dependent on variance of the 2 populations, unlikely to be able to reuse.
      if (equal_variance) {
        // Rprintf("equal variance dof: %f \n", df);
        if (test_type == PVAL_TWO_SIDED)  // 2 sided
          o = 2.0L * (1.0L - eqvar_t_cdf( fabs(tstat)));   // this is correct.
        else if (test_type == PVAL_LESS)  // less.
          o = eqvar_t_cdf( tstat ); 
        else if (test_type == PVAL_GREATER)  // greater.  same as LESS for -t 
          o = 1.0L - eqvar_t_cdf( tstat );
        else
          o = 1.0L;
      } else {
        // Rprintf("unequal variance dof: %f \n", df);
        t_distribution_cdf<double, false> cdf;
        if (test_type == PVAL_TWO_SIDED)  // 2 sided
          o = 2.0L * (1.0L - cdf( fabs(tstat), df));   // this is correct.
        else if (test_type == PVAL_LESS)  // less.
          o = cdf(tstat, df); 
        else if (test_type == PVAL_GREATER)  // greater.  same as LESS for -t 
          o = 1.0L - cdf( tstat, df );
        else
          o = 1.0L;
      }
    }
    *out_iter = o;
    ++out_iter;
  }


}


//' Fast t-Test for dense matrix
//'
//' This implementation uses normal approximation, which works reasonably well if sample size is large (say N>=20)
//' 
//' @rdname ttest_fast
//' @param matrix an expression matrix, COLUMN-MAJOR, each col is a feature, each row a sample
//' @param labels an integer vector, each element indicating the group to which a sample belongs.
//' @param alternative 
//' \itemize{
//' \item{0} : p(two.sided)
//' \item{1} : p(less)
//' \item{2} : p(greater)
//' }
//' @param var_equal TRUE/FALSE to indicate the variance is expected to be equal
//' @param as_dataframe TRUE/FALSE - TRUE returns a dataframe, FALSE returns a matrix
//' @param threads  number of concurrent threads.
//' @return array or dataframe.  for each gene/feature, the rows for the clusters are ordered by id.
//' @name ttest_fast
//' @export
// [[Rcpp::export]]
extern SEXP ttest_fast(
    SEXP matrix, SEXP labels, 
    SEXP alternative, 
    SEXP var_equal, 
    SEXP as_dataframe,
    SEXP threads) {

  // ----------- copy to local
  // ---- input matrix
  std::vector<double> mat;
  size_t nsamples, nfeatures, nelem;
  Rcpp::StringVector features = 
    rmatrix_to_vector(matrix, mat, nsamples, nfeatures, nelem);

  // ---- label vector
  std::vector<int> lab;
  rvector_to_vector(labels, lab, nsamples);

  // get the number of unique labels.
  std::vector<int> sorted_labels;
  get_unique_clusters(lab.cbegin(), nsamples, sorted_labels);
  size_t label_count = sorted_labels.size();

  // ---- output pval matrix
  std::vector<double> pv(nfeatures * label_count);

  // ------------------------ parameter
  int type, nthreads;
  bool _as_dataframe, var_eq;
  import_r_common_params(alternative, var_equal, as_dataframe, threads,
    type, var_eq, _as_dataframe, nthreads);

  // ------------------------ compute
  omp_set_num_threads(nthreads);
  Rprintf("THREADING: using %d threads\n", nthreads);

#pragma omp parallel num_threads(nthreads)
  {
    int tid = omp_get_thread_num();
    size_t block = nfeatures / nthreads;
    size_t rem = nfeatures - nthreads * block;
    size_t offset = tid * block + (tid > rem ? rem : tid);
    int nid = tid + 1;
    size_t end = nid * block + (nid > rem ? rem : nid);

    std::unordered_map<int, gaussian_stats<double>> gaussian_sums;
    // printf("thread %d summarizing feature %ld\n", omp_get_thread_num(), offset);
    for(offset; offset < end; ++offset) {
      // directly compute matrix and res pointers.
      dense_ttest_summary(&(mat[offset * nsamples]), lab.data(), nsamples, 0.0, gaussian_sums);
      two_sample_ttest(gaussian_sums, &(pv[offset * label_count]), type, var_eq);
    }
  }


  // ------------------------ generate output
  // GET features.
  Rcpp::StringVector new_features = populate_feature_names(features, nfeatures);


  if (_as_dataframe) {
    return(Rcpp::wrap(export_de_to_r_dataframe(pv, sorted_labels, new_features)));
  } else {
    // use clust for column names.
    return (Rcpp::wrap(export_de_to_r_matrix(pv, sorted_labels, new_features)));
  }
}




// =================================

//' Fast T-Test for sparse matrix.  2 sample t-test.
//'
//' This implementation uses normal approximation, which works reasonably well if sample size is large (say N>=20)
//' 
//' @rdname sparse_ttest_fast
//' @param matrix an expression matrix, COLUMN-MAJOR, each col is a feature, each row a sample
//' @param labels an integer vector, each element indicating the group to which a sample belongs.
//' @param alternative 
//' \itemize{
//' \item{0} : p(two.sided)
//' \item{1} : p(less)
//' \item{2} : p(greater)
//' }
//' @param var_equal TRUE/FALSE to indicate the variance is expected to be equal
//' @param as_dataframe TRUE/FALSE - TRUE returns a dataframe, FALSE returns a matrix
//' @param threads  number of concurrent threads.
//' @return array or dataframe.  for each gene/feature, the rows for the clusters are ordered by id.
//' @name sparse_ttest_fast
//' @export
// [[Rcpp::export]]
extern SEXP sparse_ttest_fast(
    SEXP matrix, SEXP labels, 
    SEXP alternative, 
    SEXP var_equal, 
    SEXP as_dataframe,
    SEXP threads) {
  // Rprintf("here 1\n");

  // ----------- copy to local
  // ---- input matrix
  std::vector<double> x;
  std::vector<int> i, p;
  size_t nsamples, nfeatures, nelem;
  Rcpp::StringVector features = 
    rsparsematrix_to_vectors(matrix, x, i, p, nsamples, nfeatures, nelem);

  Rprintf("Sparse DIM: samples %lu x features %lu, non-zeros %lu\n", nsamples, nfeatures, nelem); 

  // ---- label vector
  std::vector<int> lab;
  rvector_to_vector(labels, lab, nsamples);

  // get the number of unique labels.
  std::vector<int> sorted_labels;
  get_unique_clusters(lab.cbegin(), nsamples, sorted_labels);
  size_t label_count = sorted_labels.size();

  // ---- output pval matrix
  std::vector<double> pv(nfeatures * label_count);

  // ------------------------ parameter
  int type, nthreads;
  bool _as_dataframe, var_eq;
  import_r_common_params(alternative, var_equal, as_dataframe, threads,
    type, var_eq, _as_dataframe, nthreads);

  // ------------------------- compute
  omp_set_num_threads(nthreads);
  Rprintf("THREADING: using %d threads\n", nthreads);

#pragma omp parallel num_threads(nthreads)
  {
    int tid = omp_get_thread_num();
    size_t block = nfeatures / nthreads;
    size_t rem = nfeatures - nthreads * block;
    size_t offset = tid * block + (tid > rem ? rem : tid);
    int nid = tid + 1;
    size_t end = nid * block + (nid > rem ? rem : nid);
  
    int nz_offset, nz_count;
    std::unordered_map<int, gaussian_stats<double>> gaussian_sums;

    for(offset; offset < end; ++offset) {
      nz_offset = p[offset];
      nz_count = p[offset+1] - nz_offset;
      
      // directly compute matrix and res pointers.
      // Rprintf("thread %d processing feature %d\n", omp_get_thread_num(), i);
      sparse_ttest_summary(&(x[nz_offset]), &(i[nz_offset]), nz_count,
        lab.data(), nsamples, 0.0, gaussian_sums);
      two_sample_ttest(gaussian_sums, &(pv[offset * label_count]), type, var_eq);
    }
  }


  // ----------------------- make output
  Rcpp::StringVector new_features = populate_feature_names(features, nfeatures);
  
  if (_as_dataframe) {
    return(Rcpp::wrap(export_de_to_r_dataframe(pv, sorted_labels, new_features)));
  } else {
    // use clust for column names.
    return (Rcpp::wrap(export_de_to_r_matrix(pv, sorted_labels, new_features)));
  }
}

