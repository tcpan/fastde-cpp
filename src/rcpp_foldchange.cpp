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

struct clust_info {
  double sum;
  size_t thresh_count;
  size_t count;
};

// compute base info for foldchange.  
// if expm1 is used, then it's suitable to replace seurat "data" fold change
// else it is used for "scaled.data" or default action.
template <typename IT, typename LABEL>
std::unordered_map<LABEL, clust_info> foldchange_stats(
  IT const * in, LABEL const * labels, size_t const & count,
  IT const & min_thresh = static_cast<IT>(0), 
  bool const & exponential = false) {
  std::unordered_map<LABEL, clust_info> sums;  // first is per label sum, second is count above thresh min, third is per label count
  
  if (count == 0) return sums;

  // assign rank and accumulate.
  clust_info total = { .sum = 0.0, .thresh_count = 0, .count = 0};

  LABEL key;
  IT val, vv;
  // iterate over the input and labels.
  if (exponential) {
    for (size_t i = 0; i < count; ++i) {
      key = labels[i];
      val = in[i];
      vv = std::expm1(val);

      if (sums.count(key) == 0) sums[key] = { .sum=vv, .thresh_count = (val > min_thresh), .count =1 };
      else {
          sums[key].sum += vv;
          sums[key].thresh_count += (val > min_thresh);
          ++sums[key].count;
      }  
      total.thresh_count += (val > min_thresh);
      total.sum += vv;
    }
  } else {
    for (size_t i = 0; i < count; ++i) {
      key = labels[i];
      val = in[i];

      if (sums.count(key) == 0) sums[key] = { .sum=val, .thresh_count = (val > min_thresh), .count =1 };
      else {
          sums[key].sum += val;
          sums[key].thresh_count += (val > min_thresh);
          ++sums[key].count;
      }  
      total.thresh_count += (val > min_thresh);
      total.sum += val;
    }
  }
  total.count = count;
  sums[std::numeric_limits<LABEL>::max()] = total;
  // Rprintf("%lu %lu %lu %f\n", sums.size(), total.count, total.thresh_count, total.sum);

  // for each gene/cluster combo, report fold change, and percent above threshold
  return sums;
}

// compute percentages.  This is common amongst all variants
// this is used in seurat foldchange.default 
template <typename LABEL, typename PT>
void foldchange_percents(std::unordered_map<LABEL, clust_info> const & sums, 
  PT * percent1, PT * percent2) {
  
  if (sums.size() == 0) return;

  size_t label_count = sums.size() - 1;
  size_t count = sums.at(std::numeric_limits<LABEL>::max()).count;
  size_t thresh_count = sums.at(std::numeric_limits<LABEL>::max()).thresh_count;
  // Rprintf("%lu %lu %lu %f\n", label_count, count, thresh_count, sums.at(std::numeric_limits<LABEL>::max()).sum);
  
  // for each gene/cluster combo, report fold change, and percent above threshold
  LABEL key;
  double mean1 = 0;
  double mean2 = 0;
  for (size_t i = 0; i < label_count; ++i) {
    key = i+1;
    percent1[i] = static_cast<double>(sums.at(key).thresh_count) / static_cast<double>(sums.at(key).count);
    percent2[i] = static_cast<double>(thresh_count - sums.at(key).thresh_count) / static_cast<double>(count - sums.at(key).count);
  }
}

// compute mean .  this is used for seurat "scaled.data".  no log or pseudocount applied.
//  Should not be used with expm1
template <typename LABEL, typename OT>
void foldchange_mean(std::unordered_map<LABEL, clust_info> const & sums, 
  OT * out) {

  if (sums.size() == 0) {
    return;
  }
  size_t label_count = sums.size() - 1;
  double total_sum = sums.at(std::numeric_limits<LABEL>::max()).sum;
  size_t count = sums.at(std::numeric_limits<LABEL>::max()).count;
  // Rprintf("MEAN: %lu %lu %lu %f\n", sums.size(), count, sums.at(std::numeric_limits<LABEL>::max()).thresh_count, total_sum);

  // for each gene/cluster combo, report fold change, and percent above threshold
  LABEL key;
  double mean1 = 0;
  double mean2 = 0;
  for (size_t i = 0; i < label_count; ++i) {
    key = i+1;
    mean1 = sums.at(key).sum / static_cast<double>(sums.at(key).count);
    mean2 = (total_sum - sums.at(key).sum) / static_cast<double>(count - sums.at(key).count);
    out[i] = mean1 - mean2;
  }

}


// compute log of mean with pseudocount.  this is used for seurat "data" and default fold change.
template <typename LABEL, typename OT>
void foldchange_logmean(std::unordered_map<LABEL, clust_info> const & sums, 
  OT * out, bool const & pseudocount = true, double const & base = 2.0) {
  
  if (sums.size() == 0) return;

  size_t label_count = sums.size() - 1;
  double total_sum = sums.at(std::numeric_limits<LABEL>::max()).sum;
  size_t count = sums.at(std::numeric_limits<LABEL>::max()).count;
  
  // for each gene/cluster combo, report fold change, and percent above threshold
  LABEL key;
  double mean1 = 0;
  double mean2 = 0;
  double inv_log = 1.0 / log(base);
  for (size_t i = 0; i < label_count; ++i) {
    key = i+1;
    mean1 = sums.at(key).sum / static_cast<double>(sums.at(key).count);
    mean2 = (total_sum - sums.at(key).sum) / static_cast<double>(count - sums.at(key).count);
    if (base == 2.0) {
      mean1 = log2(mean1 + pseudocount);
      mean2 = log2(mean2 + pseudocount);
    } else if (base == 10.0) {
      mean1 = log10(mean1 + pseudocount);
      mean2 = log10(mean2 + pseudocount);
    } else if (pseudocount) {
      mean1 = log1p(mean1) * inv_log;
      mean2 = log1p(mean2) * inv_log;
    } else {
      mean1 = log(mean1) * inv_log;
      mean2 = log(mean2) * inv_log;
    }
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
extern SEXP FoldChangeBatch(SEXP matrix, SEXP labels, SEXP calc_percents, SEXP fc_name, 
  SEXP use_expm1, SEXP min_threshold, 
  SEXP use_log, SEXP log_base, SEXP use_pseudocount, SEXP threads) {
  
  bool perc = Rf_asLogical(calc_percents);
  double min_thresh = REAL(min_threshold)[0];
  bool _use_expm1 = Rf_asLogical(use_expm1);
  bool _use_log = Rf_asLogical(use_log);
  double _log_base = REAL(log_base)[0];
  bool _use_pseudocount = Rf_asLogical(use_pseudocount);

  int nthreads=INTEGER(threads)[0];
  if (nthreads < 1) nthreads = 1;
  const size_t nsamples=NROW(matrix);  // n
  const size_t nfeatures=NCOL(matrix); // m
  
  // ========= count number of labels so we can allocate output
  int *label_ptr=INTEGER(labels);
  size_t label_count = 0;
  {
    std::unordered_set<int> unique_labels;
    for (size_t l = 0; l < nsamples; ++l, ++label_ptr) {
      unique_labels.insert(*label_ptr);
    }
    label_count = unique_labels.size();   // k
  }
  label_ptr = INTEGER(labels);   // reset to start.

  // ============= alloc output
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
  SET_STRING_ELT(names, 0, Rf_mkChar("fc"));   // temporary
  if (perc) {
    SET_STRING_ELT(names, 1, Rf_mkChar("pct.1"));
    SET_STRING_ELT(names, 2, Rf_mkChar("pct.2"));
  }
  Rf_setAttrib(res, R_NamesSymbol, names);

  // Rprintf("inputsize  r %lu c %lu , output r %lu c %lu \n", nsamples, nfeatures, NROW(res), NCOL(res));

  // double *matColPtr = REAL(matrix); // pointer to the current column of the matrix
  
  omp_set_num_threads(nthreads);

  // ======= compute.
#pragma omp parallel for num_threads(nthreads)
  for(size_t i=0; i < nfeatures; ++i) {
    std::unordered_map<int, clust_info> info =
      foldchange_stats(REAL(matrix) + i * nsamples, label_ptr, nsamples, min_thresh, _use_expm1);

    // if percent is required, calc
    size_t out_offset = i * label_count;
    if (perc) {
      foldchange_percents(info, REAL(p1) + out_offset, REAL(p2) + out_offset);
    }

    if (_use_log) {
      foldchange_logmean(info, REAL(fc) + out_offset, _use_pseudocount, _log_base);
    } else {
      foldchange_mean(info, REAL(fc) + out_offset);
    }
  }


  UNPROTECT(5);
  return(res);
}


// =================================
