#include <Rcpp.h>
#include <R.h>
#include <Rdefines.h>


/*
 *  copied and modified from wave/include/transform/wilcox.hpp
 *
 *  Created on: Jun 8, 2021
 *  Author: Tony Pan
 *  Affiliation: Institute for Data Engineering and Science
 *  			Georgia Institute of Technology, Atlanta, GA 30332
 */

#include <cmath>  // sqrt
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <type_traits>

// Enable C++11
// [[Rcpp::plugins(cpp11)]]

// Enable OpenMP (excludes macOS)
// [[Rcpp::plugins(openmp)]]

#include <omp.h>

// #ifndef INV_SQRT2
// #define INV_SQRT2 0.70710678118
// #endif

#include "rcpp_cluster_utils.hpp"
#include "rcpp_data_utils.hpp"
#include "rcpp_benchmark_utils.hpp"


template <typename IT, typename LABEL>
void sparse_sum_rank(std::vector<std::pair<IT, LABEL> > const & temp, 
  size_t const & count, IT const & zero_val,
  std::unordered_map<LABEL, size_t> const & z_cl_counts,
  std::unordered_map<LABEL, size_t> & rank_sums,
  double & tie_sum
) {
  rank_sums.clear();
  tie_sum = 0.0;

  size_t nzc = temp.size();
  double tiesum_scale = 1.0 / (static_cast<double>(count) * static_cast<double>(count - 1));

  // initialize
  LABEL key;
  for (auto item : z_cl_counts) {
    key = item.first;
    rank_sums[key] = 0;
  }

  double tie;

  // ---- below zero, forward walk
  std::pair<IT, LABEL> zero = {zero_val, 0};
  auto lower = std::lower_bound(temp.begin(), temp.end(), zero,
    [](std::pair<IT, LABEL> const & a, std::pair<IT, LABEL> const & b){
      return a.first < b.first;
    });
  // init loop.
  auto start = temp.begin();
  size_t rank = 0;
  
  // insert first and save value for last.
  IT last_val = std::numeric_limits<IT>::lowest();

  // forward walk loop.
  size_t i = 1;
  for (; start != lower; ++start, ++i) {
    // rank update
    if (last_val < start->first) {  // ties are dealt with only when exiting ties.
        tie = i - rank;  // rank is of the first element in the tie.
        if (i > rank + 1) {
          tie_sum += tiesum_scale * (tie * tie - 1.0) * tie;
        }
        rank = i;  // not tie, so advance to new rank
    } // ELSE tie, so don't change the rank  (note :  this is the lower of the range for median)
    // add rank to sum
    rank_sums[start->second] += rank;
    last_val = start->first;
  }

  // ----- at zeros.  here start is at y.  
  //                       (s)--------------------->s
  // in:    w   w   x   x   0   0   0   0   0   0   0   y    z    z   v
  // label: c   a   b   b   d   a   f   a   d   c   e   e    a    d   f
  // r:     1   2   3   4   5   6   7   8   9  10  11  12   13   14  15
  // rank:1.5 1.5 3.5 3.5   8   8   8   8   8   8   8  12 13.5 13.5  15
  //                        i--------------------->(i)
  // compute last tie sum before zero, and new rank.
  // this part is only called on transition to and from zero ties.
  tie = i - rank;
  if (i > rank + 1) {
    tie_sum += tiesum_scale * (tie * tie - 1.0) * tie;
  }
  rank = i;  // pointing at first zero.
  // add to rank_sum - all the "missing" elements are zeros.
  // this is called many times.   last_val during this stretch will all be 0.
  for (auto item : z_cl_counts) {
    rank_sums[item.first] += rank * item.second;  // basically, adding rank z_cl_counts[key] times.
  }
  // about to re-enter the second half of the loop.
  // set i to new location at (i)
  i += count - nzc;  // count - non-zero == zero count.
  last_val = zero_val;

  // ----- above zero
  for (; start != temp.end(); ++start, ++i) {
    // rank update
    if (last_val < start->first) {  // ties are dealt with only when exiting ties.
        tie = i - rank;  // rank is of the first element in the tie.
        if (i > rank + 1) {
          tie_sum += tiesum_scale * (tie * tie - 1.0) * tie;
        }
        rank = i;  // not tie, so advance to new rank
    } // ELSE tie, so don't change the rank  (note :  this is the lower of the range for median)
    // add rank to sum
    rank_sums[start->second] += rank;
    last_val = start->first;
  }
  // last part.
  tie = i - rank;
  if (i > rank + 1) {
    tie_sum += tiesum_scale * (tie * tie - 1.0) * tie;
  }
    
  // =====  reverse walk
  auto upper = std::lower_bound(temp.rbegin(), temp.rend(), zero, 
      [](std::pair<IT, LABEL> const & a, std::pair<IT, LABEL> const & b){
      return a.first > b.first;
    });
  auto rstart = temp.rbegin();
  rank = count;

  // insert first and save value for last.
  last_val = std::numeric_limits<IT>::max();

  // reverse walking 
  // ------ greater than zero.
  i = count;
  for (; rstart != upper; ++rstart, --i) {
    // rank update
    if (last_val > rstart->first) rank = i;  // not tie, so advance to new rank
    // ELSE tie, so don't change the rank  (note :  this is the higher of the range for median)

    // add rank to sum
    rank_sums[rstart->second] += rank;
    last_val = rstart->first;
  }

  // ----- at zeros.  here start is at y.  
  //              s<-----------(s)
  // in:    x x x x 0 0 0 0 0 0 0 y y y y
  // label: c a b b d a f a d c e e a d f
  // rank:  1 1 2 2 3 3 3 3 3 3 3 4 5 5 7
  //             (i)<-----------i
  // compute last tie sum before zero, and new rank.
  // this part is only called on transition to and from zero ties.
  rank = i;
  // add to rank_sum - all the "missing" elements are zeros.
  // this is called many times.   last_val during this stretch will all be 0.
  for (auto item : z_cl_counts) {
    rank_sums[item.first] += rank * item.second;  // basically, adding rank z_cl_counts[key] times.
  }
  // about to re-enter the second half of the loop.
  // set i to new location at (i)
  i -= count - nzc;  // count - non-zero == zero count.
  last_val = zero_val;

  for (; rstart != temp.rend(); ++rstart, --i) {
    // rank update
    if (last_val > rstart->first) rank = i;  // not tie, so advance to new rank
    // ELSE tie, so don't change the rank  (note :  this is the higher of the range for median)

    // add rank to sum
    rank_sums[rstart->second] += rank;
    last_val = rstart->first;
  }

}

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
  std::unordered_map<LABEL, size_t> & rank_sums, double & tie_sum) {
  
  rank_sums.clear();
  tie_sum = 0.0;

  // tuplize in and labels.
  if (count == 0) return;
  assert((count < 20) && "ERROR: count is too small (< 20) for a normal approximation\n");

  // initialize ranksums, also need per cluster zero counts.
  std::unordered_map<LABEL, size_t> z_cl_counts;
  LABEL key;
  for (auto item : cl_counts) {
    z_cl_counts[item.first] = item.second;
  }

  // comparator.
  // auto first_less = 
  //   [](std::pair<IT, LABEL> const & a, std::pair<IT, LABEL> const & b){
  //     return a.first < b.first;
  //     return (a.first == b.first) ? (a.second < b.second) : (a.first < b.first);

  //   };

  // ============== populate the sort structure and Sort ONLY NON-ZERO ENTRIES.  if any input is actually 0.0, skip
  std::vector<std::pair<IT, LABEL> > temp;
  temp.reserve(nz_count);
  for (size_t i = 0; i < nz_count; ++i) {
    // EXCLUDE values with 0, as they would and should tie with the elements not present.
    if (in[i] != zero_val) {
      key = labels[ids[i]];
      temp.emplace_back(in[i], key);
      --z_cl_counts[key];   // get the zero counts by subtracting non-zero entries
    }  
  }
  // sort by label...
  std::sort(temp.begin(), temp.end(), 
    [](std::pair<IT, LABEL> const & a, std::pair<IT, LABEL> const & b){
      return (a.first == b.first) ? (a.second < b.second) : (a.first < b.first);
    });

  // to handle ties, we need 0.5.  here is an approach:  
  // walk once from front to add rank, and once from back and add rank.
  //  when ties, forward walk uses the starting rank, reverse walk uses the end rank
  //  each sum /2 would give correct rank
  //    at end of tie, reset to use array index as rank.
  // alternative is likely not better:  
  //    any tie detection would need 1x iteration.  
  //    each tie would need another pass to get median so may break prefetching.  ties common?
  //    should minimize writes in the big array as well.


  // ** NEED TO SEPARATELY HANDLE ZEROS.  
  //    note that ALL ZEROS are tied for all elements.

  // == local zero crossing in the sorted list (temp), 
  //   then process 2 halves separately.
  sparse_sum_rank(temp, count, zero_val, z_cl_counts, rank_sums, tie_sum);
}


// rank_sums output:  map cluster to rank_sum.
// counts zeros and skip sorting of those.  this is faster than sorting all.
template <typename IT_ITER, typename LABEL_ITER,
  typename IT = typename std::iterator_traits<IT_ITER>::value_type,
  typename LABEL = typename std::iterator_traits<LABEL_ITER>::value_type>
void pseudosparse_wmw_summary(
  IT_ITER in, LABEL_ITER labels, size_t const & count, IT const & zero_val,
  std::vector<std::pair<LABEL, size_t> > const & cl_counts,
  std::unordered_map<LABEL, size_t> & rank_sums, double & tie_sum) {

  // tuplize in and labels.
  rank_sums.clear();
  tie_sum = 0.0;

  if (count == 0) return;
  assert((count < 20) && "ERROR: count is too small (< 20) for a normal approximation\n");

  std::unordered_map<LABEL, size_t> z_cl_counts;
  LABEL key;
  for (auto item : cl_counts) {
    z_cl_counts[item.first] = item.second;
  }


  // ============== populate the sort structure and Sort
  std::vector<std::pair<IT, LABEL> > temp;
  temp.reserve(count);
  for (size_t i = 0; i < count; ++i) {
    // EXCLUDE values with 0, as they would and should tie with the elements not present.
    if (in[i] != zero_val) {
      key = labels[i];
      temp.emplace_back(in[i], key);
      --z_cl_counts[key];   // get the zero counts by subtracting non-zero entries
    }  
  }

  // sort by value...
  std::sort(temp.begin(), temp.end(), 
    [](std::pair<IT, LABEL> const & a, std::pair<IT, LABEL> const & b){
      return (a.first == b.first) ? (a.second < b.second) : (a.first < b.first);
  });

  // assign rank and accumulate count and rank

  // to handle ties, we need 0.5.  To make comparison easier and not float based, double the rank value here and divide later.
  // actually, we may need to have median of the tie ranks.
  // here is an approach:  keep undoubled rank, walk once from front, and once from back.  if ties, do not change rank until value changes.  each sum /2 would give correct rank
  //    at end of tie, reset to use array index as rank.
  // any tie detection would need 1x iteration.   each tie would need another pass to get median so may break prefetching, but this may not be common?
  // should minimize writes in the big array as well.

  sparse_sum_rank(temp, count, zero_val, z_cl_counts, rank_sums, tie_sum);

}


// rank_sums output:  map cluster to rank_sum.
// old, do not use.
template <typename IT_ITER, typename LABEL_ITER,
  typename IT = typename std::iterator_traits<IT_ITER>::value_type,
  typename LABEL = typename std::iterator_traits<LABEL_ITER>::value_type>
void dense_wmw_summary(
  IT_ITER in, LABEL_ITER labels, size_t const & count,
  std::unordered_map<LABEL, size_t> & rank_sums, double & tie_sum) {

  // tuplize in and labels.
  rank_sums.clear();
  double tiesum_scale = 1.0 / (static_cast<double>(count) * static_cast<double>(count - 1));

  tie_sum = 0.0;
  if (count == 0) return;
  assert((count < 20) && "ERROR: count is too small (< 20) for a normal approximation\n");

  // ============== populate the sort structure and Sort
  std::vector<std::pair<IT, LABEL> > temp;
  temp.reserve(count);
  for (size_t i = 0; i < count; ++i) {
      temp.emplace_back(in[i], labels[i]);
  }
  // sort by value...
  std::sort(temp.begin(), temp.end(), 
    [](std::pair<IT, LABEL> const & a, std::pair<IT, LABEL> const & b){
      return (a.first == b.first) ? (a.second < b.second) : (a.first < b.first);
  });

  // assign rank and accumulate count and rank

  // to handle ties, we need 0.5.  To make comparison easier and not float based, double the rank value here and divide later.
  // actually, we may need to have median of the tie ranks.
  // here is an approach:  keep undoubled rank, walk once from front, and once from back.  if ties, do not change rank until value changes.  each sum /2 would give correct rank
  //    at end of tie, reset to use array index as rank.
  // any tie detection would need 1x iteration.   each tie would need another pass to get median so may break prefetching, but this may not be common?
  // should minimize writes in the big array as well.

  size_t rank = 1;
  double tie;
  LABEL key = temp[0].second;
  // insert first
  if (rank_sums.count(key) == 0) rank_sums[key] = rank;
  else {
      rank_sums[key] += rank;
  }
  // forward walk
  size_t i = 1;
  for (size_t i = 1; i < count; ++i) {
      // rank update
      if (temp[i-1].first < temp[i].first) {  // ties are dealt with only when exiting ties.
          tie = i + 1 - rank;  // rank is of the first element in the tie.
          tie_sum += tiesum_scale * (tie * tie - 1.0) * tie;
          rank = i + 1;  // not tie, so advance to new rank
      } // ELSE tie, so don't change the rank  (note :  this is the lower of the range for median)
      key = temp[i].second;
      // add rank to sum
      if (rank_sums.count(key) == 0) rank_sums[key] = rank;
      else {
          rank_sums[key] += rank;
      }
  }
  if (i > rank + 1) {
    tie = i + 1 - rank;
    tie_sum += tiesum_scale * (tie * tie - 1.0) * tie;
  }
    
  // reverse walk
  rank = count;
  i = count - 1;
  for (; i > 0; --i) {
      // add rank to sum
      rank_sums[temp[i].second] += rank;

      // rank update
      if (temp[i].first > temp[i-1].first) rank = i;  // not tie, so advance to new rank
      // ELSE tie, so don't change the rank  (note :  this is the higher of the range for median)
  }
  // need to count ties...
  rank_sums[temp[0].second] += rank;

}


// types:  0 = less, 1 = greater, 2 = twosided (default), 3 = U2
template <typename LABEL, typename OT_ITER,
  typename OT = typename std::iterator_traits<OT_ITER>::value_type>
void wmw(
  std::vector<std::pair<LABEL, size_t> > const & cl_counts, 
  std::unordered_map<LABEL, size_t> const & rank_sums, 
  double const & tie_sum,
  size_t const & count,
  OT_ITER out,
  int const & test_type, bool const & continuity) {
  
  // tuplize in and labels.
  if (count == 0) return;
  assert((count < 20) && "ERROR: count is too small (< 20) for a normal approximation\n");

  if (cl_counts.size() == 0) return;

  // dealing with ties amongst multiple classes?
  // since final differential comparisons are one vs all others, and the ties are split so 0.5 is added for all.

  // compute the 1 to other U stats, storing the minimum U
  double R, U1, U2, U, prod, mu, sigma, z;
  size_t n1;
  double c1 = static_cast<double>(count + 1);
  double tie_mean = tie_sum; // / static_cast<double>(count * (count - 1));
  constexpr double inv12 = 1.0 / 12.0;
  size_t i = 0;
  LABEL key;
  size_t val;
  for (auto item : cl_counts) {
      key = item.first;
      n1 = item.second;
      val = rank_sums.at(key);
      R = static_cast<double>(val) * 0.5;  // since val was computed from both forward and reverse walk, we need to half it.

      // compute U stats
      prod = static_cast<double>(n1 * (count - n1));
      U1 = R - static_cast<double>((n1 * (n1 + 1)) >> 1);  // same as STATISTIC in 
      //  https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/wilcox.test.R
      U2 = prod - U1; 
      if (test_type == PVAL_GREATER) { // greater
          U = U1;
      } else if (test_type == PVAL_LESS) {  // less
          U = U2;
      } else if (test_type == PVAL_TWO_SIDED) {  // two sided
          U = std::max(U1, U2);
      }
      // normal approximation
      mu = prod * 0.5;
      // sigma = sqrt((prod / 12.0) * ((count + 1) - tie_sum / (count * (count - 1))));
      sigma = sqrt(prod * inv12 * (c1 - tie_mean));
      z = U1 - mu;
      if (continuity) {
        if (test_type == PVAL_GREATER) { // greater
            z -= 0.5;
        } else if (test_type == PVAL_LESS) {  // less
            z += 0.5;
        } else if (test_type == PVAL_TWO_SIDED) {  // two sided
            z -= (z > 0) ? 0.5 : ((z < 0) ? -0.5 : 0.0);
        }
      } 
      z /= sigma;

      // convert to p-value
      // https://stackoverflow.com/questions/2328258/cumulative-normal-distribution-function-in-c-c
      // use erfc function  This is 2-tailed.
      if (test_type == PVAL_TWO_SIDED)  // 2 sided
        out[i] = erfc( fabs(z) * M_SQRT1_2 );   // this is correct.  erf is has stdev of 1/sqrt(2), and return prob(-x < X < x).   erfc = 1- erf(x)
      else if (test_type == PVAL_LESS)
        out[i] = 0.5 * erfc( -z * M_SQRT1_2 );
      else if (test_type == PVAL_GREATER)  // greater or less - U is changed (either U1 or U2), thus z is changed.  still calculating 1-CDF = survival
        out[i] = 1.0 - 0.5 * erfc( -z * M_SQRT1_2 );
      else if (test_type == TEST_VAL)
        out[i] = U1;
      else
        out[i] = z;
      ++i;

  }

}


//' Fast Wilcoxon-Mann-Whitney Test for dense matrix
//'
//' This implementation uses normal approximation, which works reasonably well if sample size is large (say N>=20)
//' 
//' @rdname wmw_fast
//' @param matrix an expression matrix, COLUMN-MAJOR, each col is a feature, each row a sample
//' @param labels an integer vector, each element indicating the group to which a sample belongs.
//' @param rtype 
//' \itemize{
//' \item{0} : p(less)
//' \item{1} : p(greater)
//' \item{2} : p(twoSided)
//' \item{3} : U
//' }
//' @param continuity_correction TRUE/FALSE for continuity correction
//' @param as_dataframe TRUE/FALSE - TRUE returns a dataframe, FALSE returns a matrix
//' @param threads  number of concurrent threads.
//' @return array or dataframe.  for each gene/feature, the rows for the clusters are ordered by id.
//' @name wmw_fast
//' @export
// [[Rcpp::export]]
extern SEXP wmw_fast(
    Rcpp::NumericMatrix const & matrix, Rcpp::IntegerVector const & labels, 
    int rtype, 
    bool continuity_correction, 
    bool as_dataframe,
    int threads) {

  // ----------- copy to local
  // ---- input matrix
  std::vector<double> mat;
  size_t nsamples, nfeatures, nelem;
  Rcpp::StringVector features = 
    copy_rmatrix_to_cppvector(matrix, mat, nsamples, nfeatures, nelem);

  // ---- label vector
  std::vector<int> lab;
  copy_rvector_to_cppvector(labels, lab, nsamples);

  // get the number of unique labels.
  std::vector<std::pair<int, size_t> > sorted_cluster_counts;
  count_clusters(lab, sorted_cluster_counts);
  size_t label_count = sorted_cluster_counts.size();

  // ---- output pval matrix
  std::vector<double> pv(nfeatures * label_count);

  // ------------------------ parameter
  // int rtype, threads;
  // bool as_dataframe, continuity_correction;
  // import_de_common_params(rtype, continuity_correction, rtype, continuity_correction);
  // import_r_common_params(as_dataframe, threads,
  //   as_dataframe, threads);

  // ------------------------ compute
  omp_set_num_threads(threads);
  Rprintf("THREADING: using %d threads\n", threads);

  std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double>> start = std::chrono::steady_clock::now();

#pragma omp parallel num_threads(threads)
  {
    int tid = omp_get_thread_num();
    size_t block = nfeatures / threads;
    size_t rem = nfeatures - threads * block;
    size_t offset = tid * block + (tid > rem ? rem : tid);
    int nid = tid + 1;
    size_t end = nid * block + (nid > rem ? rem : nid);

    std::unordered_map<int, size_t> rank_sums;
    double tie_sum;

    for(; offset < end; ++offset) {
      // directly compute matrix and res pointers.
      // Rprintf("thread %d processing feature %d\n", omp_get_thread_num(), i);
      pseudosparse_wmw_summary(&(mat[offset * nsamples]), 
        lab.data(), nsamples, 
        static_cast<double>(0),
        sorted_cluster_counts, rank_sums, tie_sum);
      wmw(sorted_cluster_counts, rank_sums, tie_sum, nsamples, &(pv[offset * label_count]), rtype, continuity_correction);
    }
  }
  Rprintf("[TIME] WMW Elapsed(ms)= %f\n", since(start).count());


  // ------------------------ generate output
  // GET features.
  Rcpp::StringVector new_features = populate_feature_names(features, nfeatures);

  if (as_dataframe) {
    return(Rcpp::wrap(export_de_to_r_dataframe(pv, "p_val", sorted_cluster_counts, new_features)));
  } else {
    // use clust for column names.
    return (Rcpp::wrap(export_de_to_r_matrix(pv, sorted_cluster_counts, new_features)));
  }
}


// =================================

//' Fast Wilcoxon-Mann-Whitney Test for dense matrix
//'
//' This implementation uses normal approximation, which works reasonably well if sample size is large (say N>=20)
//' 
//' @rdname sparse_wmw_fast
//' @param matrix an expression matrix, COLUMN-MAJOR, each col is a feature, each row a sample
//' @param labels an integer vector, each element indicating the group to which a sample belongs.
//' @param rtype 
//' \itemize{
//' \item{0} : p(less)
//' \item{1} : p(greater)
//' \item{2} : p(twoSided)
//' \item{3} : U
//' }
//' @param continuity_correction TRUE/FALSE for continuity_correction correction
//' @param as_dataframe TRUE/FALSE - TRUE returns a dataframe, FALSE returns a matrix
//' @param threads  number of concurrent threads.
//' @return array or dataframe.  for each gene/feature, the rows for the clusters are ordered by id.
//' @name sparse_wmw_fast
//' @export
// [[Rcpp::export]]
extern SEXP sparse_wmw_fast(
    Rcpp::dgCMatrix const & matrix, Rcpp::IntegerVector const & labels, 
    int rtype, 
    bool continuity_correction, 
    bool as_dataframe,
    int threads) {
  // Rprintf("here 1\n");

  // ----------- copy to local
  // ---- input matrix
  std::vector<double> x;
  std::vector<int> i, p;
  size_t nsamples, nfeatures, nelem;
  Rcpp::StringVector features = 
    copy_rsparsematrix_to_cppvectors(matrix, x, i, p, nsamples, nfeatures, nelem);

  // Rprintf("Sparse DIM: samples %lu x features %lu, non-zeros %lu\n", nsamples, nfeatures, nelem); 

  // ---- label vector
  std::vector<int> lab;
  copy_rvector_to_cppvector(labels, lab, nsamples);

  // get the number of unique labels.
  std::vector<std::pair<int, size_t> > sorted_cluster_counts;
  count_clusters(lab, sorted_cluster_counts);
  size_t label_count = sorted_cluster_counts.size();

  // ---- output pval matrix
  std::vector<double> pv(nfeatures * label_count);

  // ------------------------ parameter
  // int rtype, threads;
  // bool as_dataframe, continuity_correction;
  // import_de_common_params(rtype, continuity_correction, rtype, continuity_correction);
  // import_r_common_params(as_dataframe, threads,
  //   as_dataframe, threads);


  // ------------------------- compute
  omp_set_num_threads(threads);
  Rprintf("THREADING: using %d threads\n", threads);

  std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double>> start = std::chrono::steady_clock::now();

#pragma omp parallel num_threads(threads)
  {
    int tid = omp_get_thread_num();
    size_t block = nfeatures / threads;
    size_t rem = nfeatures - threads * block;
    size_t offset = tid * block + (tid > rem ? rem : tid);
    int nid = tid + 1;
    size_t end = nid * block + (nid > rem ? rem : nid);

    int nz_offset, nz_count;
    std::unordered_map<int, size_t> rank_sums;
    double tie_sum;

    for(; offset < end; ++offset) {
      nz_offset = p[offset];
      nz_count = p[offset+1] - nz_offset;

      // directly compute matrix and res pointers.
      // Rprintf("thread %d processing feature %d\n", omp_get_thread_num(), i);
      sparse_wmw_summary(&(x[nz_offset]), &(i[nz_offset]), nz_count,
        lab.data(), nsamples, 
        static_cast<double>(0),
        sorted_cluster_counts, rank_sums, tie_sum);
      wmw(sorted_cluster_counts, rank_sums, tie_sum, nsamples, &(pv[offset * label_count]), rtype, continuity_correction);
    }
  }
  Rprintf("[TIME] WMW Elapsed(ms)= %f\n", since(start).count());

  // ----------------------- make output
  Rcpp::StringVector new_features = populate_feature_names(features, nfeatures);

  if (as_dataframe) {
    return(Rcpp::wrap(export_de_to_r_dataframe(pv, "p_val", sorted_cluster_counts, new_features)));
  } else {
    // use clust for column names.
    return (Rcpp::wrap(export_de_to_r_matrix(pv, sorted_cluster_counts, new_features)));
  }
}



// =================================



// =================================

//' Fast Wilcoxon-Mann-Whitney Test for dense matrix
//'
//' This implementation uses normal approximation, which works reasonably well if sample size is large (say N>=20)
//' 
//' @rdname sparse64_wmw_fast
//' @param matrix an expression matrix, COLUMN-MAJOR, each col is a feature, each row a sample
//' @param labels an integer vector, each element indicating the group to which a sample belongs.
//' @param rtype 
//' \itemize{
//' \item{0} : p(less)
//' \item{1} : p(greater)
//' \item{2} : p(twoSided)
//' \item{3} : U
//' }
//' @param continuity_correction TRUE/FALSE for continuity_correction correction
//' @param as_dataframe TRUE/FALSE - TRUE returns a dataframe, FALSE returns a matrix
//' @param threads  number of concurrent threads.
//' @return array or dataframe.  for each gene/feature, the rows for the clusters are ordered by id.
//' @name sparse64_wmw_fast
//' @export
// [[Rcpp::export]]
extern SEXP sparse64_wmw_fast(
    Rcpp::dgCMatrix64 const & matrix, Rcpp::IntegerVector const & labels, 
    int rtype, 
    bool continuity_correction, 
    bool as_dataframe,
    int threads) {
  // Rprintf("here 1\n");

  // ----------- copy to local
  // ---- input matrix
  std::vector<double> x;
  std::vector<long> i, p;
  size_t nsamples, nfeatures, nelem;
  Rcpp::StringVector features = 
    copy_rsparsematrix_to_cppvectors(matrix, x, i, p, nsamples, nfeatures, nelem);

  // Rprintf("Sparse DIM: samples %lu x features %lu, non-zeros %lu\n", nsamples, nfeatures, nelem); 

  // ---- label vector
  std::vector<int> lab;
  copy_rvector_to_cppvector(labels, lab, nsamples);

  // get the number of unique labels.
  std::vector<std::pair<int, size_t> > sorted_cluster_counts;
  count_clusters(lab, sorted_cluster_counts);
  size_t label_count = sorted_cluster_counts.size();

  // ---- output pval matrix
  std::vector<double> pv(nfeatures * label_count);

  // ------------------------ parameter
  // int rtype, threads;
  // bool as_dataframe, continuity_correction;
  // import_de_common_params(rtype, continuity_correction, rtype, continuity_correction);
  // import_r_common_params(as_dataframe, threads,
  //   as_dataframe, threads);


  // ------------------------- compute
  omp_set_num_threads(threads);
  Rprintf("THREADING: using %d threads\n", threads);

  std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double>> start = std::chrono::steady_clock::now();

#pragma omp parallel num_threads(threads)
  {
    int tid = omp_get_thread_num();
    size_t block = nfeatures / threads;
    size_t rem = nfeatures - threads * block;
    size_t offset = tid * block + (tid > rem ? rem : tid);
    int nid = tid + 1;
    size_t end = nid * block + (nid > rem ? rem : nid);

    long nz_offset, nz_count;
    std::unordered_map<int, size_t> rank_sums;
    double tie_sum;

    for(; offset < end; ++offset) {
      nz_offset = p[offset];
      nz_count = p[offset+1] - nz_offset;

      // directly compute matrix and res pointers.
      // Rprintf("thread %d processing feature %d\n", omp_get_thread_num(), i);
      sparse_wmw_summary(&(x[nz_offset]), &(i[nz_offset]), nz_count,
        lab.data(), nsamples, 
        static_cast<double>(0),
        sorted_cluster_counts, rank_sums, tie_sum);
      wmw(sorted_cluster_counts, rank_sums, tie_sum, nsamples, &(pv[offset * label_count]), rtype, continuity_correction);
    }
  }
  Rprintf("[TIME] WMW Elapsed(ms)= %f\n", since(start).count());

  // ----------------------- make output
  Rcpp::StringVector new_features = populate_feature_names(features, nfeatures);

  if (as_dataframe) {
    return(Rcpp::wrap(export_de_to_r_dataframe(pv, "p_val", sorted_cluster_counts, new_features)));
  } else {
    // use clust for column names.
    return (Rcpp::wrap(export_de_to_r_matrix(pv, sorted_cluster_counts, new_features)));
  }
}
