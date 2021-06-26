#include <Rcpp.h>
#include <R.h>
#include <Rdefines.h>
using namespace Rcpp;


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

// Enable C++11
// [[Rcpp::plugins(cpp11)]]

// Enable OpenMP (excludes macOS)
// [[Rcpp::plugins(openmp)]]

#include <omp.h>

// #ifndef INV_SQRT2
// #define INV_SQRT2 0.70710678118
// #endif

#include "common.hpp"

// ids points to start of the column's positive element row ids
// x  points to the start fo the columns positive element values
// count is hte number of positive elements in the column.
template <typename IT, typename IDX, typename LABEL>
void sparse_wmw_summary(IT const * in, IDX const * ids, size_t const & nz_count, 
  LABEL const * labels, size_t const & count, 
  std::map<LABEL, size_t> const & cl_counts,
  std::unordered_map<LABEL, size_t> & rank_sums, double & tie_sum) {
  
  rank_sums.clear();
  tie_sum = 0.0;

  // tuplize in and labels.
  if (count == 0) return;

  assert((count < 20) && "ERROR: count is too small (< 20) for a normal approximation\n");
  // if (count < 100) {
  //     printf("WARNING: count is small for a normal approximation: %ld\n", count);
  // }

  // initialize ranksums, also need per cluster zero counts.
  std::unordered_map<LABEL, size_t> z_cl_counts;
  LABEL key;
  for (auto item : cl_counts) {
    key = item.first;
    rank_sums[key] = 0;
    z_cl_counts[key] = item.second;
  }

  // comparator.
  auto first_less = 
    [](std::pair<IT, LABEL> const & a, std::pair<IT, LABEL> const & b){
      return a.first < b.first;
    };

  // ============== populate the sort structure and Sort ONLY NON-ZERO ENTRIES.  if any input is actually 0.0, skip
  std::vector<std::pair<IT, LABEL>> temp;
  temp.reserve(nz_count);
  for (size_t i = 0; i < nz_count; ++i) {
    // EXCLUDE values with 0, as they would and should tie with the elements not present.
    if (in[i] != static_cast<IT>(0)) {
      key = labels[ids[i]];
      temp.emplace_back(in[i], key);
      --z_cl_counts[key];   // get the zero counts by subtracting non-zero entries
    }  
  }
  size_t nzc = temp.size();
  // sort by label...
  std::sort(temp.begin(), temp.end(), first_less);

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

  double tie;

  // ---- below zero, forward walk
  std::pair<IT, LABEL> zero = {0, 0};
  auto lower = std::lower_bound(temp.begin(), temp.end(), zero, first_less);
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
        if (tie > 1) tie_sum += tie * tie * tie - tie;
        rank = i;  // not tie, so advance to new rank
    } // ELSE tie, so don't change the rank  (note :  this is the lower of the range for median)
    // add rank to sum
    rank_sums[start->second] += rank;
    last_val = start->first;
  }

  // ----- at zeros.  here start is at y.  
  //               (s)----------->s
  // in:    x x x x 0 0 0 0 0 0 0 y y y y
  // label: c a b b d a f a d c e e a d f
  // rank:  1 1 2 2 3 3 3 3 3 3 3 4 5 5 7
  //                i----------->(i)
  // compute last tie sum before zero, and new rank.
  // this part is only called on transition to and from zero ties.
  tie = i - rank;
  if (tie > 1) tie_sum += tie * tie * tie - tie;
  rank = i;
  // add to rank_sum - all the "missing" elements are zeros.
  // this is called many times.   last_val during this stretch will all be 0.
  for (auto item : z_cl_counts) {
    rank_sums[item.first] += rank * item.second;  // basically, adding rank z_cl_counts[key] times.
  }
  // about to re-enter the second half of the loop.
  // set i to new location at (i)
  i += count - nzc;  // count - non-zero == zero count.
  last_val = 0;

  // ----- above zero
  for (; start != temp.end(); ++start, ++i) {
    // rank update
    if (last_val < start->first) {  // ties are dealt with only when exiting ties.
        tie = i - rank;  // rank is of the first element in the tie.
        if (tie > 1) tie_sum += tie * tie * tie - tie;
        rank = i;  // not tie, so advance to new rank
    } // ELSE tie, so don't change the rank  (note :  this is the lower of the range for median)
    // add rank to sum
    rank_sums[start->second] += rank;
    last_val = start->first;
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
  last_val = 0;

  for (; rstart != temp.rend(); ++rstart, --i) {
    // rank update
    if (last_val > rstart->first) rank = i;  // not tie, so advance to new rank
    // ELSE tie, so don't change the rank  (note :  this is the higher of the range for median)

    // add rank to sum
    rank_sums[rstart->second] += rank;
    last_val = rstart->first;
  }

}

// rank_sums output:  map cluster to rank_sum.
template <typename IT, typename LABEL>
void dense_wmw_summary(
  IT const * in, LABEL const * labels, size_t const & count, 
  std::unordered_map<LABEL, size_t> & rank_sums, double & tie_sum) {

  // tuplize in and labels.
  rank_sums.clear();
  tie_sum = 0.0;
  if (count == 0) return;
  assert((count < 20) && "ERROR: count is too small (< 20) for a normal approximation\n");

  // ============== populate the sort structure and Sort
  std::vector<std::pair<IT, LABEL>> temp;
  temp.reserve(count);
  for (size_t i = 0; i < count; ++i) {
      temp.emplace_back(in[i], labels[i]);
  }
  // sort by value...
  std::sort(temp.begin(), temp.end(), 
    [](std::pair<IT, LABEL> const & a, std::pair<IT, LABEL> const & b){
      return a.first < b.first;
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
  tie_sum = 0.0;
  LABEL key = temp[0].second;
  // insert first
  if (rank_sums.count(key) == 0) rank_sums[key] = rank;
  else {
      rank_sums[key] += rank;
  }
  // forward walk
  for (size_t i = 1; i < count; ++i) {
      // rank update
      if (temp[i-1].first < temp[i].first) {  // ties are dealt with only when exiting ties.
          tie = i + 1 - rank;  // rank is of the first element in the tie.
          tie_sum += tie * tie * tie - tie;
          rank = i + 1;  // not tie, so advance to new rank
      } // ELSE tie, so don't change the rank  (note :  this is the lower of the range for median)
      key = temp[i].second;
      // add rank to sum
      if (rank_sums.count(key) == 0) rank_sums[key] = rank;
      else {
          rank_sums[key] += rank;
      }
  }
    
  // reverse walk
  rank = count;
  for (size_t i = count - 1; i > 0; --i) {
      // add rank to sum
      rank_sums[temp[i].second] += rank;

      // rank update
      if (temp[i].first > temp[i-1].first) rank = i;  // not tie, so advance to new rank
      // ELSE tie, so don't change the rank  (note :  this is the higher of the range for median)
  }
  // need to count ties...
  rank_sums[temp[0].second] += rank;

}


// types:  0 = greater, 1 = less, 2 = twosided (default), 3 = U
template <typename LABEL, typename OT>
void wmw(
  std::map<LABEL, size_t> const & cl_counts, 
  std::unordered_map<LABEL, size_t> const & rank_sums, double const & tie_sum,
  size_t const & count,
  OT * out,
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
  double tie_mean = tie_sum / static_cast<double>(count * (count - 1));
  constexpr double inv12 = 1.0 / 12.0;
  size_t i = 0;
  LABEL key;
  size_t val;
  for (auto item : cl_counts) {
      key = item.first;
      n1 = item.second;
      val = rank_sums.at(key);
      R = static_cast<double>(val) * 0.5;

      // compute U stats
      prod = static_cast<double>(n1 * (count - n1));
      U1 = R - static_cast<double>((n1 * (n1 + 1)) >> 1);
      U2 = prod - U1; 
      if (test_type == 0)  // greater
          U = U1;
      else if (test_type == 1)   // less
          U = U2;
      else   // two sided
          U = std::max(U1, U2);

      // normal approximation
      mu = prod * 0.5;
      // sigma = sqrt((prod / 12.0) * ((count + 1) - tie_sum / (count * (count - 1))));
      sigma = sqrt(prod * inv12 * (c1 - tie_mean));
      z = U - mu - (continuity ? 0.5 : 0.0);
      z /= sigma;

      // convert to p-value
      // https://stackoverflow.com/questions/2328258/cumulative-normal-distribution-function-in-c-c
      // use erfc function  This is 2-tailed.
      if (test_type == 2)  // 2 sided
          out[i] = erfc( z * M_SQRT1_2 );   // this is correct.
      else if (test_type == 1)   // less
          out[i] = 1.0 - 0.5 * erfc( -z * M_SQRT1_2 );
      else if (test_type == 0)   // greater
          out[i] = 1.0 - 0.5 * erfc( -z * M_SQRT1_2 );  // this is 1-CDF = survival function
      else
          out[i] = std::min(U1, U2);
      ++i;

  }

}


//' Fast Wilcoxon-Mann-Whitney Test for dense matrix
//'
//' This implementation uses normal approximation, which works reasonably well if sample size is large (say N>=20)
//' 
//' @rdname wmwfast
//' @param matrix an expression matrix, COLUMN-MAJOR, each col is a feature, each row a sample
//' @param labels an integer vector, each element indicating the group to which a sample belongs.
//' @param rtype 
//' \itemize{
//' \item{0} : p(greater)
//' \item{1} : p(less)
//' \item{2} : p(twoSided)
//' \item{3} : U
//' }
//' @param continuity_correction TRUE/FALSE for continuity correction
//' @param as_dataframe TRUE/FALSE - TRUE returns a dataframe, FALSE returns a matrix
//' @param threads  number of concurrent threads.
//' @return array or dataframe.  for each gene/feature, the rows for the clusters are ordered by id.
//' @name wmwfast
//' @export
// [[Rcpp::export]]
extern SEXP wmwfast(
    SEXP matrix, SEXP labels, SEXP rtype, 
    SEXP continuity_correction, 
    SEXP as_dataframe,
    SEXP threads) {
  // Rprintf("here 1\n");

  const int type=Rf_asInteger(rtype);
  if (type > 3) Rprintf("ERROR: unsupported type: %d. Supports only greater, less, twosided, and U\n", type);
  int nthreads=Rf_asInteger(threads);
  if (nthreads < 1) nthreads = 1;
  bool continuity = Rf_asLogical(continuity_correction);
  const size_t nsamples=NROW(matrix);  // n
  const size_t nfeatures=NCOL(matrix); // m
  bool _as_dataframe = Rf_asLogical(as_dataframe);
  // Rprintf("here 2\n");
  
  // get the number of unique labels.
  int *label_ptr=INTEGER(labels);
  std::map<int, size_t> cl_counts;
  count_clusters(label_ptr, nsamples, cl_counts);
  size_t label_count = cl_counts.size();
  label_ptr = INTEGER(labels);
  // Rprintf("here 3\n");

  int intlen = snprintf(NULL, 0, "%lu", std::numeric_limits<size_t>::max());
  char * str = reinterpret_cast<char *>(malloc(intlen + 1));
  int proc_count = 0;

  // https://stackoverflow.com/questions/23547625/returning-a-data-frame-from-c-to-r

  // GET features.
  SEXP features = Rf_getAttrib(matrix, R_NamesSymbol);
  // check if features is null.  if so, make a new one.
  // https://stackoverflow.com/questions/25172419/how-can-i-get-the-sexptype-of-an-sexp-value
  if (TYPEOF(features) == NILSXP) {
    PROTECT(features = Rf_allocVector(STRSXP, nfeatures));
    ++proc_count;

    for (size_t j = 0; j < nfeatures; ++j) {
      // create string and set in clust.
      sprintf(str, "%lu", j);
      SET_STRING_ELT(features, j, Rf_mkChar(str));
      memset(str, 0, intlen + 1);
    }

  }
  // Rprintf("here 4\n");

  // double *matColPtr; // pointer to the current column of the matrix
  SEXP res, pv, clust, genenames, names, rownames, cls;
  int ncols = 1 + (_as_dataframe ? 2 : 0);
  int col_id = 0;

  if (_as_dataframe) {
    PROTECT(res = Rf_allocVector(VECSXP, ncols));   // cluster id and gene names are repeated.
    proc_count += 1;

    PROTECT(clust = Rf_allocVector(INTSXP, label_count * nfeatures));
    PROTECT(genenames = Rf_allocVector(STRSXP, label_count * nfeatures));
    PROTECT(pv = Rf_allocVector(REALSXP, label_count * nfeatures));
    proc_count += 3;

    PROTECT(cls = Rf_allocVector(STRSXP, 1)); // class attribute
    SET_STRING_ELT(cls, 0, Rf_mkChar("data.frame"));
    Rf_classgets(res, cls);
    ++proc_count;

  } else {
    // use clust for column names.
    PROTECT(clust = Rf_allocVector(STRSXP, label_count));
    ++proc_count;
    size_t j = 0;
    for (auto item : cl_counts) {
      sprintf(str, "%d", item.first);
      SET_STRING_ELT(clust, j, Rf_mkChar(str));
      memset(str, 0, intlen + 1);
      ++j;
    }

    PROTECT(pv = Rf_allocMatrix(REALSXP, label_count, nfeatures));
    ++proc_count;

  }
  // Rprintf("here 5\n");

  // alloc output res (list?)
  if (_as_dataframe) {
    PROTECT(rownames = Rf_allocVector(STRSXP, label_count * nfeatures));  // dataframe column names.
    ++proc_count;

    // make the clusters vector.
    int * clust_ptr = INTEGER(clust);
    SEXP * features_ptr = STRING_PTR(features);
    size_t j = 0;
    for (size_t i = 0; i < nfeatures; ++i) {
      for (auto item : cl_counts) {
        // rotate through cluster labels for this feature.        
        *clust_ptr = item.first;
        ++clust_ptr;

        // same feature name for the set of cluster
        SET_STRING_ELT(genenames, j, Rf_duplicate(*features_ptr));

        sprintf(str, "%lu", j);
        SET_STRING_ELT(rownames, j, Rf_mkChar(str));
        memset(str, 0, intlen + 1);
  
        ++j;
      }
      ++features_ptr;

    }

    PROTECT(names = Rf_allocVector(STRSXP, ncols));  // dataframe column names.
    ++proc_count;

    SET_VECTOR_ELT(res, col_id, clust);  
    SET_STRING_ELT(names, col_id, Rf_mkChar("cluster"));
    ++col_id;
    SET_VECTOR_ELT(res, col_id, genenames);
    SET_STRING_ELT(names, col_id, Rf_mkChar("gene"));
    ++col_id;
    SET_VECTOR_ELT(res, col_id, pv);
    SET_STRING_ELT(names, col_id, Rf_mkChar("p_val"));   
    ++col_id;
    Rf_namesgets(res, names);  // colnames

    // set row names - NEEDED to print the dataframe!!!
    Rf_setAttrib(res, R_RowNamesSymbol, rownames);


  } else {
    // set col and row names for pv.
    // https://stackoverflow.com/questions/5709940/r-extension-in-c-setting-matrix-row-column-names
    PROTECT(names = Rf_allocVector(VECSXP, 2));
    ++proc_count;
    SET_VECTOR_ELT(names, 0, clust);  // rows = clusters
    SET_VECTOR_ELT(names, 1, features);  // columns  = features (genes)
    Rf_setAttrib(pv, R_DimNamesSymbol, names);
  }
  // Rprintf("here 6\n");

  // Rprintf("inputsize  r %lu c %lu , output r %lu c %lu \n", nsamples, nfeatures, NROW(res), NCOL(res));

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

    double *matptr = REAL(matrix) + offset * nsamples;
    double *pvptr = REAL(pv) + offset * label_count;
    std::unordered_map<int, size_t> rank_sums;
    double tie_sum;
    for(size_t i=offset; i < end; ++i) {
      // directly compute matrix and res pointers.
      // Rprintf("thread %d processing feature %d\n", omp_get_thread_num(), i);
      dense_wmw_summary(matptr, label_ptr, nsamples, rank_sums, tie_sum);
      wmw(cl_counts, rank_sums, tie_sum, nsamples, pvptr, type, continuity);
      matptr += nsamples;
      pvptr += label_count;
    }
  }
  // Rprintf("here 7\n");

  UNPROTECT(proc_count);
  free(str);
  if (_as_dataframe) return(res);
  else return(pv);
}


// =================================

//' Fast Wilcoxon-Mann-Whitney Test for dense matrix
//'
//' This implementation uses normal approximation, which works reasonably well if sample size is large (say N>=20)
//' 
//' @rdname sparsewmwfast
//' @param matrix an expression matrix, COLUMN-MAJOR, each col is a feature, each row a sample
//' @param labels an integer vector, each element indicating the group to which a sample belongs.
//' @param rtype 
//' \itemize{
//' \item{0} : p(greater)
//' \item{1} : p(less)
//' \item{2} : p(twoSided)
//' \item{3} : U
//' }
//' @param continuity_correction TRUE/FALSE for continuity correction
//' @param as_dataframe TRUE/FALSE - TRUE returns a dataframe, FALSE returns a matrix
//' @param threads  number of concurrent threads.
//' @return array or dataframe.  for each gene/feature, the rows for the clusters are ordered by id.
//' @name sparsewmwfast
//' @export
// [[Rcpp::export]]
extern SEXP sparsewmwfast(
    SEXP matrix, SEXP labels, SEXP rtype, 
    SEXP continuity_correction, 
    SEXP as_dataframe,
    SEXP threads) {
  // Rprintf("here 1\n");

  const int type=Rf_asInteger(rtype);
  if (type > 3) Rprintf("ERROR: unsupported type: %d. Supports only greater, less, twosided, and U\n", type);
  int nthreads=Rf_asInteger(threads);
  if (nthreads < 1) nthreads = 1;

  bool continuity = Rf_asLogical(continuity_correction);
  bool _as_dataframe = Rf_asLogical(as_dataframe);
  // Rprintf("here 2\n");
  
  S4 obj(matrix);
  SEXP i = obj.slot("i");
  SEXP p = obj.slot("p");  // ncol + 1
  SEXP x = obj.slot("x");
  SEXP dim = obj.slot("Dim");
  SEXP dimnms = obj.slot("Dimnames");
  // SEXP rownms = VECTOR_ELT(dimnms, 0);    // samples
  SEXP features = VECTOR_ELT(dimnms, 1);   // features_names = columnames

  size_t nsamples = INTEGER(dim)[0];   // sample count
  size_t nfeatures = INTEGER(dim)[1];   // feature/gene count
  size_t nelem = INTEGER(p)[nfeatures];   // since p is offsets, the ncol+1 entry has the total count.
  

  // get the number of unique labels.
  int *label_ptr=INTEGER(labels);
  std::map<int, size_t> cl_counts;
  count_clusters(label_ptr, nsamples, cl_counts);
  size_t label_count = cl_counts.size();
  label_ptr = INTEGER(labels);

  int intlen = snprintf(NULL, 0, "%lu", std::numeric_limits<size_t>::max());
  char * str = reinterpret_cast<char *>(malloc(intlen + 1));
  int proc_count = 0;

  // https://stackoverflow.com/questions/23547625/returning-a-data-frame-from-c-to-r

  // GET features.
  // check if features is null.  if so, make a new one.
  // https://stackoverflow.com/questions/25172419/how-can-i-get-the-sexptype-of-an-sexp-value
  if (TYPEOF(features) == NILSXP) {
    PROTECT(features = Rf_allocVector(STRSXP, nfeatures));
    ++proc_count;

    for (size_t j = 0; j < nfeatures; ++j) {
      // create string and set in clust.
      sprintf(str, "%lu", j);
      SET_STRING_ELT(features, j, Rf_mkChar(str));
      memset(str, 0, intlen + 1);
    }

  }
  // Rprintf("here 4\n");

  // double *matColPtr; // pointer to the current column of the matrix
  SEXP res, pv, clust, genenames, names, rownames, cls;
  int ncols = 1 + (_as_dataframe ? 2 : 0);
  int col_id = 0;

  if (_as_dataframe) {
    PROTECT(res = Rf_allocVector(VECSXP, ncols));   // cluster id and gene names are repeated.
    proc_count += 1;

    PROTECT(clust = Rf_allocVector(INTSXP, label_count * nfeatures));
    PROTECT(genenames = Rf_allocVector(STRSXP, label_count * nfeatures));
    PROTECT(pv = Rf_allocVector(REALSXP, label_count * nfeatures));
    proc_count += 3;

    PROTECT(cls = Rf_allocVector(STRSXP, 1)); // class attribute
    SET_STRING_ELT(cls, 0, Rf_mkChar("data.frame"));
    Rf_classgets(res, cls);
    ++proc_count;

  } else {
    // use clust for column names.
    PROTECT(clust = Rf_allocVector(STRSXP, label_count));
    ++proc_count;
    size_t j = 0;
    for (auto item : cl_counts) {
      sprintf(str, "%d", item.first);
      SET_STRING_ELT(clust, j, Rf_mkChar(str));
      memset(str, 0, intlen + 1);
      ++j;
    }

    PROTECT(pv = Rf_allocMatrix(REALSXP, label_count, nfeatures));
    ++proc_count;

  }
  // Rprintf("here 5\n");

  // alloc output res (list?)
  if (_as_dataframe) {
    PROTECT(rownames = Rf_allocVector(STRSXP, label_count * nfeatures));  // dataframe column names.
    ++proc_count;

    // make the clusters vector.
    int * clust_ptr = INTEGER(clust);
    SEXP * features_ptr = STRING_PTR(features);
    size_t j = 0;
    for (size_t i = 0; i < nfeatures; ++i) {
      for (auto item : cl_counts) {
        // rotate through cluster labels for this feature.        
        *clust_ptr = item.first;
        ++clust_ptr;

        // same feature name for the set of cluster
        SET_STRING_ELT(genenames, j, Rf_duplicate(*features_ptr));

        sprintf(str, "%lu", j);
        SET_STRING_ELT(rownames, j, Rf_mkChar(str));
        memset(str, 0, intlen + 1);
  
        ++j;
      }
      ++features_ptr;

    }

    PROTECT(names = Rf_allocVector(STRSXP, ncols));  // dataframe column names.
    ++proc_count;

    SET_VECTOR_ELT(res, col_id, clust);  
    SET_STRING_ELT(names, col_id, Rf_mkChar("cluster"));
    ++col_id;
    SET_VECTOR_ELT(res, col_id, genenames);
    SET_STRING_ELT(names, col_id, Rf_mkChar("gene"));
    ++col_id;
    SET_VECTOR_ELT(res, col_id, pv);
    SET_STRING_ELT(names, col_id, Rf_mkChar("p_val"));   
    ++col_id;
    Rf_namesgets(res, names);  // colnames

    // set row names - NEEDED to print the dataframe!!!
    Rf_setAttrib(res, R_RowNamesSymbol, rownames);


  } else {
    // set col and row names for pv.
    // https://stackoverflow.com/questions/5709940/r-extension-in-c-setting-matrix-row-column-names
    PROTECT(names = Rf_allocVector(VECSXP, 2));
    ++proc_count;
    SET_VECTOR_ELT(names, 0, clust);  // rows = clusters
    SET_VECTOR_ELT(names, 1, features);  // columns  = features (genes)
    Rf_setAttrib(pv, R_DimNamesSymbol, names);
  }
  // Rprintf("here 6\n");

  // Rprintf("inputsize  r %lu c %lu , output r %lu c %lu \n", nsamples, nfeatures, NROW(res), NCOL(res));

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

    double *mat_ptr = REAL(x);
    int * rowid_ptr = INTEGER(i);
    int * offset_ptr = INTEGER(p);
    int nz_offset, nz_count;

    double *pvptr = REAL(pv) + offset * label_count;
    std::unordered_map<int, size_t> rank_sums;
    double tie_sum;

    for(size_t i=offset; i < end; ++i) {
      nz_offset = offset_ptr[i];
      nz_count = offset_ptr[i+1] - nz_offset;

      // directly compute matrix and res pointers.
      // Rprintf("thread %d processing feature %d\n", omp_get_thread_num(), i);
      sparse_wmw_summary(mat_ptr + nz_offset, rowid_ptr + nz_offset, nz_count,
      label_ptr, nsamples, cl_counts, rank_sums, tie_sum);
      wmw(cl_counts, rank_sums, tie_sum, nsamples, pvptr, type, continuity);
      pvptr += label_count;
    }
  }
  // Rprintf("here 7\n");

  UNPROTECT(proc_count);
  free(str);
  if (_as_dataframe) return(res);
  else return(pv);
}
