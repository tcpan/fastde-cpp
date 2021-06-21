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

  size_t count = sums.at(std::numeric_limits<LABEL>::max()).count;
  size_t thresh_count = sums.at(std::numeric_limits<LABEL>::max()).thresh_count;
  // Rprintf("%lu %lu %lu %f\n", label_count, count, thresh_count, sums.at(std::numeric_limits<LABEL>::max()).sum);
  
  // for each gene/cluster combo, report fold change, and percent above threshold
  double mean1 = 0;
  double mean2 = 0;
  size_t i = 0;
  for (auto item : sums) {
    if (item.first == std::numeric_limits<LABEL>::max()) continue;  // skip the max count

    percent1[i] = static_cast<double>(item.second.thresh_count) / static_cast<double>(item.second.count);
    percent2[i] = static_cast<double>(thresh_count - item.second.thresh_count) / static_cast<double>(count - item.second.count);
    ++i;
    // order of entry is same as order of sums traversal.
  }
  
}

// compute mean .  this is used for seurat "scaled.data".  no log or pseudocount applied.
//  Should not be used with expm1
template <typename LABEL, typename OT>
void foldchange_mean(std::unordered_map<LABEL, clust_info> const & sums, 
  OT * out) {

  if (sums.size() == 0) return;
 
  double total_sum = sums.at(std::numeric_limits<LABEL>::max()).sum;
  size_t count = sums.at(std::numeric_limits<LABEL>::max()).count;
  // Rprintf("MEAN: %lu %lu %lu %f\n", sums.size(), count, sums.at(std::numeric_limits<LABEL>::max()).thresh_count, total_sum);

  // for each gene/cluster combo, report fold change, and percent above threshold
  double mean1 = 0;
  double mean2 = 0;
  size_t i = 0;
  for (auto item : sums) {    
    if (item.first == std::numeric_limits<LABEL>::max()) continue;  // skip the max count

    mean1 = item.second.sum / static_cast<double>(item.second.count);
    mean2 = (total_sum - item.second.sum) / static_cast<double>(count - item.second.count);
    out[i] = mean1 - mean2;
    ++i;
  }

}


// compute log of mean with pseudocount.  this is used for seurat "data" and default fold change.
template <typename LABEL, typename OT>
void foldchange_logmean(std::unordered_map<LABEL, clust_info> const & sums, 
  OT * out, bool const & pseudocount = true, double const & base = 2.0) {
  
  if (sums.size() == 0) return;

  double total_sum = sums.at(std::numeric_limits<LABEL>::max()).sum;
  size_t count = sums.at(std::numeric_limits<LABEL>::max()).count;
  
  // for each gene/cluster combo, report fold change, and percent above threshold
  double mean1 = 0;
  double mean2 = 0;
  double inv_log = 1.0 / log(base);
  size_t i = 0;
  for (auto item : sums) {    
    if (item.first == std::numeric_limits<LABEL>::max()) continue;  // skip the max count

    mean1 = item.second.sum / static_cast<double>(item.second.count);
    mean2 = (total_sum - item.second.sum) / static_cast<double>(count - item.second.count);
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
    ++i;
  }
}




/*! \brief Fold Change
 *
 * \param matrix: an expression matrix, each row is a feature, each column corresponds to a samples
 * \param labels: an integer vector, each element indicating the group to which a sample belongs.
 * \param calc_percents:  a boolean to indicate whether to compute percents or not.
 * \param fc_name: column name to use for the fold change results 
 * \param use_expm1: for "data", use expm1
 * \param min_threshold: minimum threshold to count towards pct.1 and pct.2 percentages.
 * \param use_log: for "data" and default log type, indicate log of the sum is to be used.
 * \param log_base: base for the log
 * \param use_pseudocount: for "data" and default log type, add pseudocount after log.
 * \param as_dataframe: TRUE/FALSE.  TRUE = return a linearized dataframe.  FALSE = return matrices.
 * \param threads: number of threads to use
 * 
 * https://stackoverflow.com/questions/38338270/how-to-return-a-named-vecsxp-when-writing-r-extensions
 */
// [[Rcpp::export]]
extern SEXP ComputeFoldChange(SEXP matrix, SEXP labels, SEXP calc_percents, SEXP fc_name, 
  SEXP use_expm1, SEXP min_threshold, 
  SEXP use_log, SEXP log_base, SEXP use_pseudocount, 
  SEXP as_dataframe,
  SEXP threads) {
  
  bool perc = Rf_asLogical(calc_percents);
  double min_thresh = REAL(min_threshold)[0];
  bool _use_expm1 = Rf_asLogical(use_expm1);
  bool _use_log = Rf_asLogical(use_log);
  double _log_base = REAL(log_base)[0];
  bool _use_pseudocount = Rf_asLogical(use_pseudocount);
  bool _as_dataframe = Rf_asLogical(as_dataframe);

  int nthreads=INTEGER(threads)[0];
  if (nthreads < 1) nthreads = 1;
  const size_t nsamples=NROW(matrix);  // n
  const size_t nfeatures=NCOL(matrix); // m
  
  // ========= count number of labels so we can allocate output  run the first one.
  int *label_ptr=INTEGER(labels);
  // run 1, so we can get the cluster order.
  std::unordered_map<int, clust_info> info =
      foldchange_stats(REAL(matrix), label_ptr, nsamples, min_thresh, false);
  size_t label_count = info.size() - 1;   // k
  label_ptr = INTEGER(labels);   // reset to start.
  
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

  // ============= alloc output
  SEXP fc, p1, p2, clust, genenames, res, names, rownames, cls, dimnames;
  int ncols = 1 + (perc ? 2 : 0) + (_as_dataframe ? 2 : 0);
  int col_id = 0;

  PROTECT(res = Rf_allocVector(VECSXP, ncols));   // cluster id and gene names are repeated.
  PROTECT(names = Rf_allocVector(STRSXP, ncols));
  proc_count += 2;

  if (_as_dataframe) {
    PROTECT(clust = Rf_allocVector(INTSXP, label_count * nfeatures));
    PROTECT(genenames = Rf_allocVector(STRSXP, label_count * nfeatures));
    PROTECT(fc = Rf_allocVector(REALSXP, label_count * nfeatures));
    proc_count += 3;
    if (perc) {
      PROTECT(p1 = Rf_allocVector(REALSXP, label_count * nfeatures));
      PROTECT(p2 = Rf_allocVector(REALSXP, label_count * nfeatures)); 
      proc_count += 2;
    }

    PROTECT(cls = Rf_allocVector(STRSXP, 1)); // class attribute
    SET_STRING_ELT(cls, 0, Rf_mkChar("data.frame"));
    Rf_classgets(res, cls);

    ++proc_count;

  } else {
    // use clust for column names.
    PROTECT(clust = Rf_allocVector(STRSXP, label_count));
    ++proc_count;
    int key;
    size_t j = 0;
    for (auto item : info) {
      key = item.first;
      if (key == std::numeric_limits<int>::max()) continue;

      sprintf(str, "%d", key);
      SET_STRING_ELT(clust, j, Rf_mkChar(str));
      memset(str, 0, intlen + 1);
      ++j;
    }
    
    PROTECT(fc = Rf_allocMatrix(REALSXP, label_count, nfeatures));
    ++proc_count;
    if (perc) {
      PROTECT(p1 = Rf_allocMatrix(REALSXP, label_count, nfeatures));
      PROTECT(p2 = Rf_allocMatrix(REALSXP, label_count, nfeatures)); 
      proc_count += 2;
    }
  }

  // set up the output structure.
  if (_as_dataframe) {
    PROTECT(rownames = Rf_allocVector(STRSXP, label_count * nfeatures));  // dataframe column names.
    ++proc_count;

    // make the clusters vector.
    int * clust_ptr = INTEGER(clust);
    SEXP * features_ptr = STRING_PTR(features);
    size_t j = 0;
    for (size_t i = 0; i < nfeatures; ++i) {
      for (auto item : info) {
        if (item.first == std::numeric_limits<int>::max()) continue;

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

    SET_VECTOR_ELT(res, col_id, clust);  
    SET_STRING_ELT(names, col_id, Rf_mkChar("cluster"));
    ++col_id;
    SET_VECTOR_ELT(res, col_id, genenames);
    SET_STRING_ELT(names, col_id, Rf_mkChar("gene"));
    ++col_id;
  }


  SET_VECTOR_ELT(res, col_id, fc);
  // convert from STRSXP to CHARSXP
  SET_STRING_ELT(names, col_id, Rf_mkChar(CHAR(STRING_ELT(fc_name, 0)))); 
  ++col_id;

  if (perc) {
    SET_VECTOR_ELT(res, col_id, p1);
    SET_STRING_ELT(names, col_id, Rf_mkChar("pct.1"));
    ++col_id;
    SET_VECTOR_ELT(res, col_id, p2);
    SET_STRING_ELT(names, col_id, Rf_mkChar("pct.2"));
    ++col_id;
  }

  // column names.
  if (_as_dataframe) {
    Rf_namesgets(res, names);

    // set row names - NEEDED to print the dataframe!!!
    Rf_setAttrib(res, R_RowNamesSymbol, rownames);
  } else {
    // set list element names.
    Rf_setAttrib(res, R_NamesSymbol, names);

    PROTECT(dimnames = Rf_allocVector(VECSXP, 2));
    ++proc_count;
    SET_VECTOR_ELT(dimnames, 0, clust);  // rows = clusters
    SET_VECTOR_ELT(dimnames, 1, features);  // columns  = features (genes)

    // https://stackoverflow.com/questions/5709940/r-extension-in-c-setting-matrix-row-column-names
    Rf_setAttrib(fc, R_DimNamesSymbol, dimnames);
    if (perc) {
      Rf_setAttrib(p1, R_DimNamesSymbol, Rf_duplicate(dimnames));
      Rf_setAttrib(p2, R_DimNamesSymbol, Rf_duplicate(dimnames));
    }

  }

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

  UNPROTECT(proc_count);
  free(str);
  return(res);
}


/*! \brief Filter based on FoldChange
 *
 * \param fc: foldchange values, either as a vector or a matrix
 * \param pct1: percent greater than threshold (0) in class 1.
 * \param pct2: percent greater than threshold (0) in class 2.
 * \param init_mask:  initial mask, based on an external "features" vector.
 * \param min_pct: minimum threshold for max pct1, pct2
 * \param min_diff_pct: minimum threshold for difference between max and min {pct1, pct2}
 * \param logfc_threshold: if not scaled.data, then compare to the logfc.
 * \param only_pos: keep only positive fc value, and not use abs when thresholding.
 * \param not_count:  not scaled.data
 * \param threads: number of threads to use
 * 
 *  https://stackoverflow.com/questions/38338270/how-to-return-a-named-vecsxp-when-writing-r-extensions
 */
// [[Rcpp::export]]
extern SEXP FilterFoldChange(SEXP fc, SEXP pct1, SEXP pct2,
  SEXP init_mask,
  SEXP min_pct, SEXP min_diff_pct, SEXP logfc_threshold, 
  SEXP only_pos, SEXP not_count,
  SEXP threads) {
  

  // feature_mask <- ifelse(
  //   test = is.null(features),
  //   yes = ! logical(nrow(data)),   # all rows
  //   no = rownames(data) %in% features   # just the one specified.
  // )
  double _min_pct = REAL(min_pct)[0];
  double _min_diff_pct = REAL(min_diff_pct)[0];
  double _logfc_threshold = REAL(logfc_threshold)[0];

  bool _only_pos = Rf_asLogical(only_pos);  // operate on 
  bool _not_count = Rf_asLogical(not_count);

  int nthreads=INTEGER(threads)[0];
  if (nthreads < 1) nthreads = 1;

  bool has_init = (TYPEOF(init_mask) != NILSXP);

  // ========= count number of labels so we can allocate output  run the first one.  
  const size_t nclusters=NROW(fc);  // n
  const size_t nfeatures=NCOL(fc); // m

  int proc_count = 0;

  omp_set_num_threads(nthreads);

  SEXP mask;
  int *init_ptr = LOGICAL(init_mask);
  double *fc_ptr = REAL(fc);
  double *pct1_ptr = REAL(pct1);
  double *pct2_ptr = REAL(pct2);
  int *mask_ptr = LOGICAL(mask);
  double mx, mn;
  bool out;
  size_t offset;
  if (nfeatures > 1) {
    PROTECT(mask = Rf_allocMatrix(LGLSXP, nclusters, nfeatures));
    ++proc_count;

    // ======= compute.
#pragma omp parallel for num_threads(nthreads) private(init_ptr, fc_ptr, pct1_ptr, pct2_ptr, mask_ptr, mx, mn, out, offset)
    for(size_t i=0; i < nfeatures; ++i) {
      offset = i * nclusters;
      init_ptr = LOGICAL(init_mask) + offset;
      fc_ptr = REAL(fc) + offset;
      pct1_ptr = REAL(pct1) + offset;
      pct2_ptr = REAL(pct2) + offset;
      mask_ptr = LOGICAL(mask) + offset;

      for(size_t j=0; j < nclusters; ++j) {
        out = (has_init ? *init_ptr : true);

        // alpha.min <- pmax(fc.results$pct.1, fc.results$pct.2)
        // features <- names(x = which(x = alpha.min >= min.pct))
        mx = std::max(*pct1_ptr, *pct2_ptr);
        out &= (mx >= _min_pct);
        // alpha.diff <- alpha.min - pmin(fc.results$pct.1, fc.results$pct.2)
        // x = which(x = alpha.min >= min.pct & alpha.diff >= min.diff.pct)
        mn = std::min(*pct1_ptr, *pct2_ptr);
        out &= ((mx - mn) >= _min_diff_pct);


        // if (slot != "scale.data") {
        //   total.diff <- fc.results[, 1] #first column is logFC
        //   features.diff <- if (only.pos) {
        //     names(x = which(x = total.diff >= logfc.threshold))
        //   } else {
        //     names(x = which(x = abs(x = total.diff) >= logfc.threshold))
        //   }
        mn = _only_pos ? *fc_ptr : fabs(*fc_ptr);
        //   features <- intersect(x = features, y = features.diff)
        // }
        if (_not_count) {
          out &= (mn >= _logfc_threshold);
        }

        // if (only.pos)
        // de.results <- de.results[de.results[, fc.name] > 0, , drop = FALSE]
        //              | _only_pos   | !_only_pos 
        //  fc_ptr > 0  | y           |   y
        //  fc_ptr <= 0 | n           |   y
        if (_only_pos) out &= (*fc_ptr > 0);

        *mask_ptr = out;

        if (has_init) ++init_ptr;
        ++fc_ptr;
        ++pct1_ptr;
        ++pct2_ptr;
        ++mask_ptr;
      }
    }

  } else {
    PROTECT(mask = Rf_allocVector(LGLSXP, nclusters));
    ++proc_count;

    // ======= compute.
    // ======= compute.
#pragma omp parallel for num_threads(nthreads) schedule(static) private(init_ptr, fc_ptr, pct1_ptr, pct2_ptr, mask_ptr, mx, mn, out)
    for(size_t j=0; j < nclusters; ++j) {
        out = (has_init ? init_ptr[j] : true);
        // alpha.min <- pmax(fc.results$pct.1, fc.results$pct.2)
        // features <- names(x = which(x = alpha.min >= min.pct))
        mx = std::max(pct1_ptr[j], pct2_ptr[j]);
        out &= (mx >= _min_pct);
        // alpha.diff <- alpha.min - pmin(fc.results$pct.1, fc.results$pct.2)
        // x = which(x = alpha.min >= min.pct & alpha.diff >= min.diff.pct)
        mn = std::min(pct1_ptr[j], pct2_ptr[j]);
        out &= ((mx - mn) >= _min_diff_pct);

        // if (slot != "scale.data") {
        //   total.diff <- fc.results[, 1] #first column is logFC
        //   features.diff <- if (only.pos) {
        //     names(x = which(x = total.diff >= logfc.threshold))
        //   } else {
        //     names(x = which(x = abs(x = total.diff) >= logfc.threshold))
        //   }
        mx = fc_ptr[j];
        mn = _only_pos ? mx : fabs(mx);
        //   features <- intersect(x = features, y = features.diff)
        // }
        if (_not_count) {
          out &= (mn >= _logfc_threshold);
        }

        // if (only.pos)
        // de.results <- de.results[de.results[, fc.name] > 0, , drop = FALSE]
        //              | _only_pos   | !_only_pos 
        //  fc_ptr > 0  | y           |   y
        //  fc_ptr <= 0 | n           |   y
        if (_only_pos) out &= (mx > 0);

        mask_ptr[j] = out;
    }

  }
  // set the row and col names.
  // https://stackoverflow.com/questions/5709940/r-extension-in-c-setting-matrix-row-column-names
  Rf_setAttrib(mask, R_DimNamesSymbol, Rf_getAttrib(fc, R_DimNamesSymbol));
  
  UNPROTECT(proc_count);
  return(mask);
}



// =================================
