#include <Rcpp.h>
#include <R.h>
#include <Rdefines.h>

using namespace Rcpp;

// https://stackoverflow.com/questions/44795190/how-can-i-return-a-list-of-matrices-from-rcpp-to-r

/*
 *  similar structure as wilcox.hpp
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

#ifndef NROW
#define NROW(x) INTEGER(GET_DIM((x)))[0]
#endif

#ifndef NCOL
#define NCOL(x) INTEGER(GET_DIM((x)))[1]
#endif

template <typename IT, typename LABEL, typename OT, typename PT>
void foldchange(IT const * in, LABEL const * labels, size_t const & count, 
  OT * out, PT * percent1, PT * percent2, size_t const & label_count, 
  IT const & min_thresh = static_cast<IT>(0)) {
  
  if (count == 0) return;

  struct clust_info {
    double sum;
    size_t thresh_count;
    size_t count;
  };

  // assign rank and accumulate.
  std::unordered_map<LABEL, clust_info> sums;  // first is per label sum, second is count above thresh min, third is per label count
  size_t thresh_count = 0;
  double total_sum = 0.0;

  LABEL key;
  IT val;
  // iterate over the input and labels.
  for (size_t i = 0; i < count; ++i) {
    key = labels[i];
    val = in[i];

    if (sums.count(key) == 0) sums[key] = { .sum=val, .thresh_count = (val > min_thresh), .count =1 };
    else {
        sums[key].sum += val;
        sums[key].thresh_count += (val > min_thresh);
        ++sums[key].count;
    }  
    thresh_count += (val > min_thresh);
    total_sum += val;
  }

  // for each gene/cluster combo, report fold change, and percent above threshold
  double mean1 = 0;
  double mean2 = 0;
  for (size_t i = 0; i < label_count; ++i) {
    key = i+1;
    mean1 = sums[key].sum / static_cast<double>(sums[key].count);
    mean2 = (total_sum - sums[key].sum) / static_cast<double>(count - sums[key].count);
    out[i] = mean1 - mean2;

    percent1[i] = static_cast<double>(sums[key].thresh_count) / static_cast<double>(sums[key].count);
    percent2[i] = static_cast<double>(thresh_count - sums[key].thresh_count) / static_cast<double>(count - sums[key].count);
  }
}

template <typename IT, typename LABEL, typename OT>
void foldchange(IT const * in, LABEL const * labels, size_t const & count, 
  OT * out, size_t const & label_count) {
  
  if (count == 0) return;

  struct clust_info {
    double sum;
    size_t count;
  };

  // assign rank and accumulate.
  std::unordered_map<LABEL, clust_info> sums;  // first is per label sum, second is count above thresh min, third is per label count
  double total_sum = 0.0;

  LABEL key;
  IT val;
  // iterate over the input and labels.
  for (size_t i = 0; i < count; ++i) {
    key = labels[i];
    val = in[i];

    if (sums.count(key) == 0) sums[key] = { .sum = val, .count = 1 };
    else {
        sums[key].sum += val;
        ++sums[key].count;
    }  
    total_sum += val;
  }

  // for each gene/cluster combo, report fold change, and percent above threshold
  double mean1 = 0;
  double mean2 = 0;
  for (size_t i = 0; i < label_count; ++i) {
    key = i+1;
    mean1 = sums[key].sum / static_cast<double>(sums[key].count);
    mean2 = (total_sum - sums[key].sum) / static_cast<double>(count - sums[key].count);
    out[i] = mean1 - mean2;
  }
}



/*! \brief Wilcoxon-Mann-Whitney Test
 *
 * \param matrix: an expression matrix, each row is a feature, each column corresponds to a samples
 * \param labels: an integer vector, each element indicating the group to which a sample belongs.
 * \param calc_percents:  a boolean to indicate whether to compute percents or not.
 * \param threads: number of threads to use
 * https://stackoverflow.com/questions/38338270/how-to-return-a-named-vecsxp-when-writing-r-extensions
 */
// [[Rcpp::export]]
extern SEXP foldChange(SEXP matrix, SEXP labels, SEXP calc_percents, SEXP threads, SEXP fc_name, SEXP min_threshold) {
  
  bool perc = Rf_asLogical(calc_percents);

  int nthreads=INTEGER(threads)[0];
  if (nthreads < 1) nthreads = 1;
  const size_t nsamples=NROW(matrix);  // n
  const size_t nfeatures=NCOL(matrix); // m
  double min_thresh = REAL(min_threshold)[0];
  
  // ========= get the number of unique labels.
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

  // alloc output
  SEXP fc = PROTECT(Rf_allocMatrix(REALSXP, label_count, nfeatures));
  SEXP p1;
  SEXP p2;

  if (perc) {
    p1 = PROTECT(Rf_allocMatrix(REALSXP, label_count, nfeatures));
    p2 = PROTECT(Rf_allocMatrix(REALSXP, label_count, nfeatures)); 
  }

  // alloc output res.
  SEXP res = PROTECT(Rf_allocVector(VECSXP, (perc ? 3 : 1)));
  SET_VECTOR_ELT(res, 0, fc);
  if (perc) {
    SET_VECTOR_ELT(res, 1, p1);
    SET_VECTOR_ELT(res, 2, p2);
  }

  SEXP names = PROTECT(Rf_allocVector(STRSXP, (perc ? 3 :  1)));
  SET_STRING_ELT(names, 0, fc_name);
  if (perc) {
    SET_STRING_ELT(names, 1, Rf_mkChar("pct.1"));
    SET_STRING_ELT(names, 2, Rf_mkChar("pct.2"));
  }
  Rf_setAttrib(res, R_NamesSymbol, names);

  // Rprintf("inputsize  r %lu c %lu , output r %lu c %lu \n", nsamples, nfeatures, NROW(res), NCOL(res));

  // double *matColPtr = REAL(matrix); // pointer to the current column of the matrix
  
  omp_set_num_threads(nthreads);


  if (perc) {
#pragma omp parallel for num_threads(nthreads)
    for(size_t i=0; i < nfeatures; ++i) {
      size_t out_offset = i * label_count;
      foldchange(REAL(matrix) + i * nsamples, label_ptr, nsamples, 
        REAL(fc) + out_offset, REAL(p1) + out_offset, REAL(p2) + out_offset, label_count, min_thresh);
    }

  } else {
#pragma omp parallel for num_threads(nthreads)
    for(size_t i=0; i < nfeatures; ++i) {
      foldchange(REAL(matrix) + i * nsamples, label_ptr, nsamples, REAL(fc) + i * label_count, label_count);
    }
  }

  UNPROTECT(5);
  return(res);
}


// =================================
