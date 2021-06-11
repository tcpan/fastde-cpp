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

#ifndef __APPLE__
#include <omp.h>
#endif

#ifndef INV_SQRT2
#define INV_SQRT2 0.70710678118
#endif
#define NROW(x) INTEGER(GET_DIM((x)))[0]
#define NCOL(x) INTEGER(GET_DIM((x)))[1]

// types:  0 = greater, 1 = less, 2 = twosided (default), 3 = U
template <typename IT, typename LABEL, typename OT>
void wmw(IT const * in, LABEL const * labels, size_t const & count, OT * out, size_t const & label_count, int const & type) {
  
  // tuplize in and labels.
  if (count == 0) return;
  assert((count < 20) && "ERROR: count is too small (< 20) for a normal approximation\n");
  // if (count < 100) {
  //     printf("WARNING: count is small for a normal approximation: %ld\n", count);
  // }

  // populate the sort structure.
  std::vector<std::pair<IT, LABEL>> temp;
  temp.reserve(count);
  for (size_t i = 0; i < count; ++i) {
      temp.emplace_back(in[i], labels[i]);
  }

  // sort 
  std::sort(temp.begin(), temp.end(), [](std::pair<IT, LABEL> const & a, std::pair<IT, LABEL> const & b){
      return a.first < b.first;
  });

  // assign rank and accumulate.
  std::unordered_map<LABEL, std::pair<size_t, size_t>> rank_sums;  // random access of this...

  // to handle ties, we need 0.5.  To make comparison easier and not float based, double the rank value here and divide later.
  // actually, we may need to have median of the tie ranks.
  // here is an approach:  keep undoubled rank, walk once from front, and once from back.  if ties, do not change rank until value changes.  each sum /2 would give correct rank
  //    at end of tie, reset to use array index as rank.
  // any tie detection would need 1x iteration.   each tie would need another pass to get median so may break prefetching, but this may not be common?
  // should minimize writes in the big array as well.

  size_t rank = 1;
  double tie, tie_sum = 0.0;
  LABEL key = temp[0].second;
  // insert first
  if (rank_sums.count(key) == 0) rank_sums[key] = std::make_pair(1, rank);
  else {
      ++rank_sums[key].first;
      rank_sums[key].second += rank;
  }
  // forward walk
  for (size_t i = 1; i < count; ++i) {
      // rank update
      if (temp[i-1].first < temp[i].first) {
          tie = i + 1 - rank;
          tie_sum += tie * tie * tie + tie;
          rank = i + 1;  // not tie, so advance to new rank
      } // ELSE tie, so don't change the rank  (note :  this is the lower of the range for median)
      key = temp[i].second;
      // add rank to sum
      if (rank_sums.count(key) == 0) rank_sums[key] = std::make_pair(1, rank);
      else {
          ++rank_sums[key].first;
          rank_sums[key].second += rank;
      }
  }
    
  // reverse walk
  rank = count;
  for (size_t i = count - 1; i > 0; --i) {
      // add rank to sum
      rank_sums[temp[i].second].second += rank;

      // rank update
      if (temp[i].first > temp[i-1].first) rank = i;  // not tie, so advance to new rank
      // ELSE tie, so don't change the rank  (note :  this is the higher of the range for median)
  }
  // need to count ties...
  rank_sums[temp[0].second].second += rank;

    
  // dealing with ties amongst multiple classes?
  // since final differential comparisons are one vs all others, and the ties are split so 0.5 is added for all.

  // compute the 1 to other U stats, storing the minimum U
  double R, U, prod, mu, sigma, z;
  size_t c;
  double c1 = static_cast<double>(count + 1);
  double cc_1 = static_cast<double>(count * (count - 1));
  constexpr double inv12 = 1.0 / 12.0;
  for (size_t i = 0; i < label_count; ++i) {
      key = i+1;
      R = static_cast<double>(rank_sums[key].second) * 0.5;
      c = rank_sums[key].first;

      // compute U stats
      U = R - static_cast<double>((c * (c + 1)) >> 1); 
      prod = static_cast<double>(c * (count - c));
      U = std::min(U, prod - U);

      // normal approximation
      mu = prod * 0.5;
      // sigma = sqrt((prod / 12.0) * ((count + 1) - tie_sum / (count * (count - 1))));
      sigma = sqrt(prod * inv12 * (c1 - tie_sum / cc_1));
      z = (U - mu) / sigma;

      // convert to p-value
      // https://stackoverflow.com/questions/2328258/cumulative-normal-distribution-function-in-c-c
      // use erfc function  This is 2-tailed.
      if (type == 2) 
        out[i] = erfc(fabs(z) * INV_SQRT2 );
      else if (type == 0)
        out[i] = erfc(z);
      else if (type == 1)
        out[i] = erf(z);
      else
        out[i] = U;
  }

}




/*! \brief Wilcoxon-Mann-Whitney Test
 *
 * \param matrix: an expression matrix, each row is a feature, each column corresponds to a samples
 * \param labels: an integer vector, each element indicating the group to which a sample belongs.
 * \param rtype:
 * \parblock
 * Define f(x)=abs(log10(x))
 * 0=p(greater), 1=p(less), 2=p(twoSided), 3=U,
 * 4=f(p(greater)),5=p(less),6=f(p(twoSided)), 7=p(greater)<p(less) ? f(p(twoSided)) : -f(p(twoSided))
 * \endparblock
 * 
 * This implementation uses normal approximation, which works reasonably well if sample size is large (say N>=20)
 */
// [[Rcpp::export]]
extern SEXP wmwtest(SEXP matrix, SEXP labels, SEXP rtype, SEXP threads) {
  const int type=INTEGER(rtype)[0];
  if (type > 3) Rprintf("ERROR: unsupported type: %d. Supports only greater, less, twosided, and U\n", type);
  const int nthreads=INTEGER(threads)[0];
  if (nthreads < 1) nthreads = 1;
  const size_t nsamples=NROW(matrix);  // n
  const size_t nfeatures=NCOL(matrix); // m
  
  // get the number of unique labels.
  size_t label_count = 0;
  int *label_ptr=INTEGER(labels);
  {
    std::unordered_set<int> unique_labels;
    for (size_t l = 0; l < nsamples; ++l, ++label_ptr) {
      unique_labels.insert(*label_ptr);
    }
    label_count = unique_labels.size();   // k
  }
  label_ptr = INTEGER(labels);

  double *matColPtr; // pointer to the current column of the matrix
  SEXP res;
  double *resPtr;  // pointer to the current column of the output 
   

  res=PROTECT(Rf_allocMatrix(REALSXP, label_count, nfeatures));
  // Rprintf("inputsize  r %lu c %lu , output r %lu c %lu \n", nsamples, nfeatures, NROW(res), NCOL(res));

  resPtr=REAL(res);
  matColPtr=REAL(matrix);

#pragma omp parallel for num_threads(nthreads)
  for(size_t i=0; i < nfeatures; ++i) {
    wmw(matColPtr, label_ptr, nsamples, resPtr, label_count, type);
    // if (i == 0) {
    //   for (size_t j = 0; j < nsamples; ++j) {
    //     Rprintf("%d, ", label_ptr[j]);
    //   }
    //   Rprintf("\n");
    //   for (size_t j = 0; j < nsamples; ++j) {
    //     Rprintf("%f, ", matColPtr[j]);
    //   }
    //   Rprintf("\n");
    // }
    resPtr += label_count;
    matColPtr += nsamples;
  }

  UNPROTECT(1);
  return(res);
}


// =================================
