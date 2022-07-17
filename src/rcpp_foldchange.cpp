#include <Rcpp.h>
#include <R.h>
#include <Rdefines.h>


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
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <chrono>  // clock, duration

// Enable C++11
// [[Rcpp::plugins(cpp11)]]

// Enable OpenMP (excludes macOS)
// [[Rcpp::plugins(openmp)]]

#include <omp.h>

#include "rcpp_cluster_utils.hpp"
#include "rcpp_data_utils.hpp"
#include "rcpp_benchmark_utils.hpp"

struct clust_info {
  double sum;
  size_t thresh_count;
};

// compute base info for foldchange.  
// if expm1 is used, then it's suitable to replace seurat "data" fold change
// else it is used for "scaled.data" or default action.
// NOTE this has an extra entry for the total for the whole "in" array.
template <typename IT_ITER, typename IDX_ITER, typename LABEL_ITER,
  typename IT = typename std::iterator_traits<IT_ITER>::value_type,
  typename IDX = typename std::iterator_traits<IDX_ITER>::value_type,
  typename LABEL = typename std::iterator_traits<LABEL_ITER>::value_type>
void sparse_foldchange_summary(
  IT_ITER in, IDX_ITER row_ids, size_t const & nz_count,
  LABEL_ITER labels, size_t const & count, IT const & zero_val,
  std::vector<std::pair<LABEL, size_t> > const & cl_counts,
  std::unordered_map<LABEL, clust_info> & sums,
  IT const & min_thresh = static_cast<IT>(0), 
  bool const & exponential = false) {
  sums.clear();  // first is per label sum, second is count above thresh min, third is per label count
  

  // ======= first handle the non-zero counts.  also locally track the local cl counts.

  // initialize.   z_cl_count will have the total number of implicit entries per cluster.
  std::unordered_map<LABEL, size_t> z_cl_counts;
  LABEL key;
  for (auto item : cl_counts) {
    key = item.first;
    sums[key] = { .sum = 0.0, .thresh_count = 0 };
    z_cl_counts[key] = item.second;
  }

  // iterate over the input and labels.
  IDX rowid;
  IT val, vv;
  for (size_t i = 0; i < nz_count; ++i) {
    val = in[i];
    vv = exponential ? std::expm1(val) : val;  // exp(val) - 1.  if val==0, then this return 0.
    rowid = row_ids[i];
    key = labels[rowid];

    sums[key].sum += vv;
    sums[key].thresh_count += (val > min_thresh);

    // if min_thresh is smaller than zero value, then include the "zero" entries of the sparse matrix.
    // if here, non-zero entry.  if min_thresh < zero_value, then zeros are included in the result, so we'd need to count zero, by counting non-zeros.
    if (val != zero_val)
      --z_cl_counts[key];
  }

  // ========= next deal with zero counts.  this affects sum when zero-val is not zero.
  if (zero_val > min_thresh) {   // if min_thresh >= 0.0, then in[i] does not contribute.
    vv = exponential ? std::expm1(zero_val) : zero_val;
    for (auto item : cl_counts) {
      key = item.first;
      sums[key].thresh_count += z_cl_counts[key];  // add all the implicit (zero) entries for the cluster.
      sums[key].sum += z_cl_counts[key] * vv;
    }

  }

  clust_info total = { .sum = 0.0, .thresh_count = 0 };
  for (auto item : sums) {
    total.thresh_count += item.second.thresh_count;
    total.sum += item.second.sum;
  }

  // now finally set the total.
  sums[std::numeric_limits<LABEL>::max()] = total;
  // Rprintf("%lu %lu %f\n", sums.size(), total.thresh_count, total.sum);

}

// compute base info for foldchange.  
// if expm1 is used, then it's suitable to replace seurat "data" fold change
// else it is used for "scaled.data" or default action.
// NOTE this has an extra entry for the total for the whole "in" array.
template <typename IT_ITER, typename LABEL_ITER,
  typename IT = typename std::iterator_traits<IT_ITER>::value_type,
  typename LABEL = typename std::iterator_traits<LABEL_ITER>::value_type>
void pseudosparse_foldchange_summary(
  IT_ITER in, LABEL_ITER labels, size_t const & count, IT const & zero_val,
  std::vector<std::pair<LABEL, size_t> > const & cl_counts,
  std::unordered_map<LABEL, clust_info> & sums,
  IT const & min_thresh = static_cast<IT>(0), 
  bool const & exponential = false) {
  sums.clear();  // first is per label sum, second is count above thresh min, third is per label count
  
  if (count == 0) return;

  // assign rank and accumulate.
  std::unordered_map<LABEL, size_t> z_cl_counts;
  LABEL key;
  for (auto item : cl_counts) {
    key = item.first;
    sums[key] = { .sum = 0.0, .thresh_count = 0 };
    z_cl_counts[key] = 0;
  }

  IT val, vv;
  // iterate over the input and labels.
  for (size_t i = 0; i < count; ++i) {
    key = labels[i];
    val = in[i];

    if (val != zero_val) {
      vv = exponential ? std::expm1(val) : val;

      sums[key].sum += vv;
      sums[key].thresh_count += (val > min_thresh);
    } else {
      ++z_cl_counts[key];
    }
  }

  if (zero_val > min_thresh) {   // if min_thresh >= 0.0, then in[i] does not contribute.
    vv = exponential ? std::expm1(zero_val) : zero_val;
    for (auto item : cl_counts) {
      key = item.first;
      sums[key].thresh_count += z_cl_counts[key];  // add all the implicit (zero) entries for the cluster.
      sums[key].sum += z_cl_counts[key] * vv;
    }

  }

  clust_info total = { .sum = 0.0, .thresh_count = 0 };
  for (auto item : sums) {
    total.thresh_count += item.second.thresh_count;
    total.sum += item.second.sum;
  }
  sums[std::numeric_limits<LABEL>::max()] = total;
  // Rprintf("%lu %lu %f\n", sums.size(), total.thresh_count, total.sum);

}

// compute percentages.  This is common amongst all variants
// this is used in seurat foldchange.default 
template <typename LABEL, typename PT_ITER, typename PT =
  typename std::iterator_traits<PT_ITER>::value_type >
void foldchange_percents(
  std::vector<std::pair<LABEL, size_t> > const & cl_counts, 
  std::unordered_map<LABEL, clust_info> const & sums, 
  size_t const & count,
  PT_ITER percent1, PT_ITER percent2) {
  
  if (count == 0) return;
  if (cl_counts.size() == 0) return;

  size_t thresh_count = sums.at(std::numeric_limits<LABEL>::max()).thresh_count;
  // Rprintf("%lu %lu %lu %f\n", label_count, count, thresh_count, sums.at(std::numeric_limits<LABEL>::max()).sum);

  // for each gene/cluster combo, report fold change, and percent above threshold
  size_t i = 0;
  size_t cl_count;
  LABEL key;
  clust_info val;
  for (auto item : cl_counts) {
    key = item.first;
    cl_count = item.second;
    val = sums.at(key);
    // percent1[i] = round(static_cast<double>(val.thresh_count) / static_cast<double>(cl_count) * 1000.0) * 0.001;
    // percent2[i] = round(static_cast<double>(thresh_count - val.thresh_count) / static_cast<double>(count - cl_count) * 1000.0) * 0.001;
    percent1[i] = static_cast<double>(val.thresh_count) / static_cast<double>(cl_count);
    percent2[i] = static_cast<double>(thresh_count - val.thresh_count) / static_cast<double>(count - cl_count);
    ++i;
    // order of entry is same as order of sums traversal.
  }
  
}

// compute mean .  this is used for seurat "scaled.data".  no log or pseudocount applied.
//  Should not be used with expm1
template <typename LABEL, typename OT_ITER, typename OT =
  typename std::iterator_traits<OT_ITER>::value_type>
void foldchange_mean(
  std::vector<std::pair<LABEL, size_t> > const & cl_counts, 
  std::unordered_map<LABEL, clust_info> const & sums, 
  size_t const & count,
  OT_ITER out) {

  if (count == 0) return;
  if (cl_counts.size() == 0) return;
  
  double total_sum = sums.at(std::numeric_limits<LABEL>::max()).sum;
  // Rprintf("MEAN: %lu %lu %lu %f\n", sums.size(), count, sums.at(std::numeric_limits<LABEL>::max()).thresh_count, total_sum);

  // for each gene/cluster combo, report fold change, and percent above threshold
  double mean1 = 0;
  double mean2 = 0;
  size_t i = 0;
  size_t cl_count;
  LABEL key;
  clust_info val;
  for (auto item : cl_counts) {
    key = item.first;
    cl_count = item.second;
    val = sums.at(key);
    mean1 = val.sum / static_cast<double>(cl_count);
    mean2 = (total_sum - val.sum) / static_cast<double>(count - cl_count);
    out[i] = mean1 - mean2;
    ++i;
  }

}


// compute log of mean with pseudocount.  this is used for seurat "data" and default fold change.
template <typename LABEL, typename OT_ITER, typename OT = 
  typename std::iterator_traits<OT_ITER>::value_type>
void foldchange_logmean(
  std::vector<std::pair<LABEL, size_t> > const & cl_counts, 
  std::unordered_map<LABEL, clust_info> const & sums, 
  size_t const & count,
  OT_ITER out, 
  bool const & pseudocount = true, double const & base = 2.0) {
  
  if (count == 0) return;
  if (cl_counts.size() == 0) return;

  double total_sum = sums.at(std::numeric_limits<LABEL>::max()).sum;
  
  // for each gene/cluster combo, report fold change, and percent above threshold
  double mean1 = 0;
  double mean2 = 0;
  double inv_log = 1.0 / log(base);
  size_t i = 0;
  size_t cl_count;
  LABEL key;
  clust_info val;
  for (auto item : cl_counts) {   
    key = item.first;
    cl_count = item.second;
    val = sums.at(key);
    mean1 = val.sum / static_cast<double>(cl_count);
    mean2 = (total_sum - val.sum) / static_cast<double>(count - cl_count);
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



//' Fold Change
//' 
//' https://stackoverflow.com/questions/38338270/how-to-return-a-named-vecsxp-when-writing-r-extensions
//' 
//' @rdname ComputeFoldChange
//' @param matrix an expression matrix, COLUMN-MAJOR, each row is a sample, each column a gene
//' @param labels an integer vector, each element indicating the group to which a sample belongs.
//' @param calc_percents  a boolean to indicate whether to compute percents or not.
//' @param fc_name column name to use for the fold change results 
//' @param use_expm1 for "data", use expm1
//' @param min_threshold minimum threshold to count towards pct.1 and pct.2 percentages.
//' @param use_log for "data" and default log type, indicate log of the sum is to be used.
//' @param log_base base for the log
//' @param use_pseudocount for "data" and default log type, add pseudocount after log.
//' @param as_dataframe TRUE/FALSE.  TRUE = return a linearized dataframe.  FALSE = return matrices.
//' @param threads number of threads to use
//' @return array or dataframe
//' @name ComputeFoldChange
//' @export
// [[Rcpp::export]]
extern SEXP ComputeFoldChange(
  Rcpp::NumericMatrix const & matrix, Rcpp::IntegerVector const & labels, 
  bool calc_percents, std::string fc_name, 
  bool use_expm1, double min_threshold, 
  bool use_log, double log_base, bool use_pseudocount, 
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
  std::vector<double> fc(nfeatures * label_count);
  std::vector<double> p1(nfeatures * label_count);
  std::vector<double> p2(nfeatures * label_count);
  
  // ------------------------ parameter
  // int threads;
  // double min_threshold, log_base;
  // bool as_dataframe, calc_percents, use_expm1, use_log, use_pseudocount;
  // import_fc_common_params(calc_percents, min_threshold, use_expm1,
  //   use_log, log_base, use_pseudocount,
  //   calc_percents, min_threshold, use_expm1, use_log, log_base, use_pseudocount);
  // import_r_common_params(as_dataframe, threads,
  //   as_dataframe, threads);
  // std::string fc_name = Rcpp::as<std::string>(fc_name);

  // ------------------------ compute
  omp_set_num_threads(threads);
  Rprintf("THREADING: using %d threads\n", threads);


  std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double>> start = std::chrono::steady_clock::now();
  // ======= compute.
#pragma omp parallel num_threads(threads)
  {
    
    int tid = omp_get_thread_num();
    size_t block = nfeatures / threads;
    size_t rem = nfeatures - threads * block;
    size_t offset = tid * block + (tid > rem ? rem : tid);
    int nid = tid + 1;
    size_t end = nid * block + (nid > rem ? rem : nid);

    std::unordered_map<int, clust_info> cl_sums;

    for(; offset < end; ++offset) {
      pseudosparse_foldchange_summary(&(mat[offset * nsamples]), 
        lab.data(), nsamples, static_cast<double>(0),
        sorted_cluster_counts, cl_sums, min_threshold, use_expm1);

      // if percent is required, calc
      if (calc_percents) {
        foldchange_percents(sorted_cluster_counts, cl_sums,
          nsamples, &(p1[offset * label_count]), &(p2[offset * label_count]));
      }

      if (use_log) {
        foldchange_logmean(sorted_cluster_counts, cl_sums,
          nsamples, &(fc[offset * label_count]), use_pseudocount, log_base);
      } else {
        foldchange_mean(sorted_cluster_counts, cl_sums, 
          nsamples, &(fc[offset * label_count]));
      }
    }
  }

  Rprintf("[TIME] FC DN Elapsed(ms)= %f\n", since(start).count());

  // ------------------------ generate output
  // GET features.
  Rcpp::StringVector new_features = populate_feature_names(features, nfeatures);

  if (as_dataframe) {
    if (calc_percents)
      return(Rcpp::wrap(export_fc_to_r_dataframe(
        fc, fc_name, 
        p1, "pct.1", 
        p2, "pct.2",
        sorted_cluster_counts, new_features)));
    else
      return(Rcpp::wrap(export_de_to_r_dataframe(fc, fc_name,
      sorted_cluster_counts, new_features)));
  } else {
    if (calc_percents)
      return (Rcpp::wrap(export_fc_to_r_matrix(
        fc, fc_name, 
        p1, "pct.1", 
        p2, "pct.2",
        sorted_cluster_counts, new_features)));
    else
    // use clust for column names.
      return (Rcpp::wrap(export_fc_to_r_matrix(fc, fc_name,
        sorted_cluster_counts, new_features)));
  }
  
}


//' Fold Change
//' 
//' https://stackoverflow.com/questions/38338270/how-to-return-a-named-vecsxp-when-writing-r-extensions
//' 
//' @rdname ComputeFoldChangeSparse
//' @param matrix an expression matrix, COLUMN-MAJOR, each row is a sample, each column a gene
//' @param labels an integer vector, each element indicating the group to which a sample belongs.
//' @param calc_percents  a boolean to indicate whether to compute percents or not.
//' @param fc_name column name to use for the fold change results 
//' @param use_expm1 for "data", use expm1
//' @param min_threshold minimum threshold to count towards pct.1 and pct.2 percentages.
//' @param use_log for "data" and default log type, indicate log of the sum is to be used.
//' @param log_base base for the log
//' @param use_pseudocount for "data" and default log type, add pseudocount after log.
//' @param as_dataframe TRUE/FALSE.  TRUE = return a linearized dataframe.  FALSE = return matrices.
//' @param threads number of threads to use
//' @return dense array or dataframe of size features*clusters
//' @name ComputeFoldChangeSparse
//' @export
// [[Rcpp::export]]
extern SEXP ComputeFoldChangeSparse(
  Rcpp::dgCMatrix const & matrix, Rcpp::IntegerVector const & labels,
  bool calc_percents, std::string fc_name, 
  bool use_expm1, double min_threshold, 
  bool use_log, double log_base, bool use_pseudocount, 
  bool as_dataframe,
  int threads) {

  // ----------- copy to local
  // ---- input matrix
  std::vector<double> x;
  std::vector<int> i, p;
  size_t nsamples, nfeatures, nelem;
  Rcpp::StringVector features = 
    copy_rsparsematrix_to_cppvectors(matrix, x, i, p, nsamples, nfeatures, nelem);

  Rprintf("Sparse DIM: samples %lu x features %lu, non-zeros %lu\n", nsamples, nfeatures, nelem); 

  // ---- label vector
  std::vector<int> lab;
  copy_rvector_to_cppvector(labels, lab, nsamples);

  // get the number of unique labels.
  std::vector<std::pair<int, size_t> > sorted_cluster_counts;
  count_clusters(lab, sorted_cluster_counts);
  size_t label_count = sorted_cluster_counts.size();

  // ---- output pval matrix
  std::vector<double> fc(nfeatures * label_count);
  std::vector<double> p1(nfeatures * label_count);
  std::vector<double> p2(nfeatures * label_count);
  
  // ------------------------ parameter
  // int threads;
  // double min_threshold, log_base;
  // bool as_dataframe, calc_percents, use_expm1, use_log, use_pseudocount;
  // import_fc_common_params(calc_percents, min_threshold, use_expm1,
  //   use_log, log_base, use_pseudocount,
  //   calc_percents, min_threshold, use_expm1, use_log, log_base, use_pseudocount);
  // import_r_common_params(as_dataframe, threads,
  //   as_dataframe, threads);
  // std::string fc_name = Rcpp::as<std::string>(fc_name);


  // ======= compute.
  
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
    std::unordered_map<int, clust_info> cl_sums;

    for(; offset < end; ++offset) {
      nz_offset = p[offset];
      nz_count = p[offset+1] - nz_offset;

      sparse_foldchange_summary(&(x[nz_offset]), &(i[nz_offset]), nz_count,
        lab.data(), nsamples,
        static_cast<double>(0),
        sorted_cluster_counts, cl_sums, min_threshold, use_expm1);

      // if percent is required, calc
      if (calc_percents) {
        foldchange_percents(sorted_cluster_counts, cl_sums, nsamples, 
        &(p1[offset * label_count]), &(p2[offset * label_count]));
      }

      if (use_log) {
        foldchange_logmean(sorted_cluster_counts, cl_sums, nsamples, &(fc[offset * label_count]), use_pseudocount, log_base);
      } else {
        foldchange_mean(sorted_cluster_counts, cl_sums, nsamples, &(fc[offset * label_count]));
      }
    }
  }
  Rprintf("[TIME] FC SP Elapsed(ms)= %f\n", since(start).count());

  // ------------------------ generate output
  // GET features.
  Rcpp::StringVector new_features = populate_feature_names(features, nfeatures);

  if (as_dataframe) {
    if (calc_percents)
      return(Rcpp::wrap(export_fc_to_r_dataframe(
        fc, fc_name, 
        p1, "pct.1", 
        p2, "pct.2",
        sorted_cluster_counts, new_features)));
    else
      return(Rcpp::wrap(export_de_to_r_dataframe(fc, fc_name,
      sorted_cluster_counts, new_features)));
  } else {
    if (calc_percents)
      return (Rcpp::wrap(export_fc_to_r_matrix(
        fc, fc_name, 
        p1, "pct.1", 
        p2, "pct.2",
        sorted_cluster_counts, new_features)));
    else
    // use clust for column names.
      return (Rcpp::wrap(export_fc_to_r_matrix(fc, fc_name,
        sorted_cluster_counts, new_features)));
  }

}


//' Fold Change
//' 
//' https://stackoverflow.com/questions/38338270/how-to-return-a-named-vecsxp-when-writing-r-extensions
//' 
//' @rdname ComputeFoldChangeSpamx32
//' @param matrix an expression matrix, COLUMN-MAJOR, each row is a sample, each column a gene
//' @param labels an integer vector, each element indicating the group to which a sample belongs.
//' @param calc_percents  a boolean to indicate whether to compute percents or not.
//' @param fc_name column name to use for the fold change results 
//' @param use_expm1 for "data", use expm1
//' @param min_threshold minimum threshold to count towards pct.1 and pct.2 percentages.
//' @param use_log for "data" and default log type, indicate log of the sum is to be used.
//' @param log_base base for the log
//' @param use_pseudocount for "data" and default log type, add pseudocount after log.
//' @param as_dataframe TRUE/FALSE.  TRUE = return a linearized dataframe.  FALSE = return matrices.
//' @param threads number of threads to use
//' @return dense array or dataframe of size features*clusters
//' @name ComputeFoldChangeSpamx32
//' @export
// [[Rcpp::export]]
extern SEXP ComputeFoldChangeSpamx32(
  Rcpp::spamx32 const & matrix, Rcpp::IntegerVector const & labels,
  bool calc_percents, std::string fc_name, 
  bool use_expm1, double min_threshold, 
  bool use_log, double log_base, bool use_pseudocount, 
  bool as_dataframe,
  int threads) {

  // ----------- copy to local
  // ---- input matrix
  std::vector<double> x;
  std::vector<int> i, p;
  size_t nsamples, nfeatures, nelem;  
  Rcpp::StringVector features = 
    copy_rsparsematrix_to_cppvectors(matrix, x, i, p, nsamples, nfeatures, nelem);

  Rprintf("Spamx32 DIM: samples %lu x features %lu, non-zeros %lu\n", nsamples, nfeatures, nelem); 
  // Rprintf("Spamx32 DIM: x %lu, p %lu, i %lu\n", x.size(), p.size(), i.size()); 

  // ---- label vector
  std::vector<int> lab;
  copy_rvector_to_cppvector(labels, lab, nsamples);

  // get the number of unique labels.
  std::vector<std::pair<int, size_t> > sorted_cluster_counts;
  count_clusters(lab, sorted_cluster_counts);
  size_t label_count = sorted_cluster_counts.size();

  // ---- output pval matrix
  std::vector<double> fc(nfeatures * label_count);
  std::vector<double> p1(nfeatures * label_count);
  std::vector<double> p2(nfeatures * label_count);
  
  // ------------------------ parameter
  // int threads;
  // double min_threshold, log_base;
  // bool as_dataframe, calc_percents, use_expm1, use_log, use_pseudocount;
  // import_fc_common_params(calc_percents, min_threshold, use_expm1,
  //   use_log, log_base, use_pseudocount,
  //   calc_percents, min_threshold, use_expm1, use_log, log_base, use_pseudocount);
  // import_r_common_params(as_dataframe, threads,
  //   as_dataframe, threads);
  // std::string fc_name = Rcpp::as<std::string>(fc_name);


  // ======= compute.
  
  omp_set_num_threads(threads);
  Rprintf("THREADING: using %d threads\n", threads);
  
  std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double>> start = std::chrono::steady_clock::now();

#pragma omp parallel num_threads(threads)
  {
    int tid = omp_get_thread_num();
    int block = nfeatures / threads;
    int rem = nfeatures - threads * block;
    int offset = tid * block + (tid > rem ? rem : tid);
    int nid = tid + 1;
    int end = nid * block + (nid > rem ? rem : nid);

    int nz_offset, nz_count;
    std::unordered_map<int, clust_info> cl_sums;

    for(; offset < end; ++offset) {
      nz_offset = p[offset];
      nz_count = p[offset+1] - nz_offset;

      sparse_foldchange_summary(&(x[nz_offset]), &(i[nz_offset]), nz_count,
        lab.data(), nsamples,
        static_cast<double>(0),
        sorted_cluster_counts, cl_sums, min_threshold, use_expm1);

      // if percent is required, calc
      if (calc_percents) {
        foldchange_percents(sorted_cluster_counts, cl_sums, nsamples, 
        &(p1[offset * label_count]), &(p2[offset * label_count]));
      }

      if (use_log) {
        foldchange_logmean(sorted_cluster_counts, cl_sums, nsamples, &(fc[offset * label_count]), use_pseudocount, log_base);
      } else {
        foldchange_mean(sorted_cluster_counts, cl_sums, nsamples, &(fc[offset * label_count]));
      }
    }
  }
  Rprintf("[TIME] FC 32 Elapsed(ms)= %f\n", since(start).count());

  // ------------------------ generate output
  // GET features.
  Rcpp::StringVector new_features = populate_feature_names(features, nfeatures);

  if (as_dataframe) {
    if (calc_percents)
      return(Rcpp::wrap(export_fc_to_r_dataframe(
        fc, fc_name, 
        p1, "pct.1", 
        p2, "pct.2",
        sorted_cluster_counts, new_features)));
    else
      return(Rcpp::wrap(export_de_to_r_dataframe(fc, fc_name,
      sorted_cluster_counts, new_features)));
  } else {
    if (calc_percents)
      return (Rcpp::wrap(export_fc_to_r_matrix(
        fc, fc_name, 
        p1, "pct.1", 
        p2, "pct.2",
        sorted_cluster_counts, new_features)));
    else
    // use clust for column names.
      return (Rcpp::wrap(export_fc_to_r_matrix(fc, fc_name,
        sorted_cluster_counts, new_features)));
  }

}


//' Fold Change
//' 
//' https://stackoverflow.com/questions/38338270/how-to-return-a-named-vecsxp-when-writing-r-extensions
//' 
//' @rdname ComputeFoldChangeSpamx64
//' @param matrix an expression matrix, COLUMN-MAJOR, each row is a sample, each column a gene
//' @param labels an integer vector, each element indicating the group to which a sample belongs.
//' @param calc_percents  a boolean to indicate whether to compute percents or not.
//' @param fc_name column name to use for the fold change results 
//' @param use_expm1 for "data", use expm1
//' @param min_threshold minimum threshold to count towards pct.1 and pct.2 percentages.
//' @param use_log for "data" and default log type, indicate log of the sum is to be used.
//' @param log_base base for the log
//' @param use_pseudocount for "data" and default log type, add pseudocount after log.
//' @param as_dataframe TRUE/FALSE.  TRUE = return a linearized dataframe.  FALSE = return matrices.
//' @param threads number of threads to use
//' @return dense array or dataframe of size features*clusters
//' @name ComputeFoldChangeSpamx64
//' @export
// [[Rcpp::export]]
extern SEXP ComputeFoldChangeSpamx64(
  Rcpp::spamx64 const & matrix, Rcpp::IntegerVector const & labels,
  bool calc_percents, std::string fc_name, 
  bool use_expm1, double min_threshold, 
  bool use_log, double log_base, bool use_pseudocount, 
  bool as_dataframe,
  int threads) {

  // ----------- copy to local
  // ---- input matrix
  std::vector<double> x;
  std::vector<long> i, p;
  size_t nsamples, nfeatures, nelem;
  Rcpp::StringVector features = 
    copy_rsparsematrix_to_cppvectors(matrix, x, i, p, nsamples, nfeatures, nelem);

  Rprintf("Sparse DIM: samples %lu x features %lu, non-zeros %lu\n", nsamples, nfeatures, nelem); 

  // ---- label vector
  std::vector<int> lab;
  copy_rvector_to_cppvector(labels, lab, nsamples);

  // get the number of unique labels.
  std::vector<std::pair<int, size_t> > sorted_cluster_counts;
  count_clusters(lab, sorted_cluster_counts);
  size_t label_count = sorted_cluster_counts.size();

  // ---- output pval matrix
  std::vector<double> fc(nfeatures * label_count);
  std::vector<double> p1(nfeatures * label_count);
  std::vector<double> p2(nfeatures * label_count);
  
  // ------------------------ parameter
  // int threads;
  // double min_threshold, log_base;
  // bool as_dataframe, calc_percents, use_expm1, use_log, use_pseudocount;
  // import_fc_common_params(calc_percents, min_threshold, use_expm1,
  //   use_log, log_base, use_pseudocount,
  //   calc_percents, min_threshold, use_expm1, use_log, log_base, use_pseudocount);
  // import_r_common_params(as_dataframe, threads,
  //   as_dataframe, threads);
  // std::string fc_name = Rcpp::as<std::string>(fc_name);


  // ======= compute.
  
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
    std::unordered_map<int, clust_info> cl_sums;

    for(; offset < end; ++offset) {
      nz_offset = p[offset];
      nz_count = p[offset+1] - nz_offset;

      sparse_foldchange_summary(&(x[nz_offset]), &(i[nz_offset]), nz_count,
        lab.data(), nsamples,
        static_cast<double>(0),
        sorted_cluster_counts, cl_sums, min_threshold, use_expm1);

      // if percent is required, calc
      if (calc_percents) {
        foldchange_percents(sorted_cluster_counts, cl_sums, nsamples, 
        &(p1[offset * label_count]), &(p2[offset * label_count]));
      }

      if (use_log) {
        foldchange_logmean(sorted_cluster_counts, cl_sums, nsamples, &(fc[offset * label_count]), use_pseudocount, log_base);
      } else {
        foldchange_mean(sorted_cluster_counts, cl_sums, nsamples, &(fc[offset * label_count]));
      }
    }
  }
  Rprintf("[TIME] FC 64 Elapsed(ms)= %f\n", since(start).count());

  // ------------------------ generate output
  // GET features.
  Rcpp::StringVector new_features = populate_feature_names(features, nfeatures);

  if (as_dataframe) {
    if (calc_percents)
      return(Rcpp::wrap(export_fc_to_r_dataframe(
        fc, fc_name, 
        p1, "pct.1", 
        p2, "pct.2",
        sorted_cluster_counts, new_features)));
    else
      return(Rcpp::wrap(export_de_to_r_dataframe(fc, fc_name,
      sorted_cluster_counts, new_features)));
  } else {
    if (calc_percents)
      return (Rcpp::wrap(export_fc_to_r_matrix(
        fc, fc_name, 
        p1, "pct.1", 
        p2, "pct.2",
        sorted_cluster_counts, new_features)));
    else
    // use clust for column names.
      return (Rcpp::wrap(export_fc_to_r_matrix(fc, fc_name,
        sorted_cluster_counts, new_features)));
  }

}



//' Filter based on FoldChange
//' 
//'  https://stackoverflow.com/questions/38338270/how-to-return-a-named-vecsxp-when-writing-r-extensions
//' 
//' @rdname FilterFoldChange
//' @param fc foldchange values, either as a vector or a matrix.  if matrix: 1 gene per col, 1 clust per row
//' @param pct1 percent greater than threshold (0) in class 1.
//' @param pct2 percent greater than threshold (0) in class 2.
//' @param init_mask  initial mask, based on an external "features" vector.
//' @param min_pct minimum threshold for max pct1, pct2
//' @param min_diff_pct minimum threshold for difference between max and min {pct1, pct2}
//' @param logfc_threshold if not scaled.data, then compare to the logfc.
//' @param only_pos keep only positive fc value, and not use abs when thresholding.
//' @param not_count  not scaled.data
//' @param threads number of threads to use
//' @return array of same shape as fc (dense)
//' @name FilterFoldChange
//' @export
// [[Rcpp::export]]
extern SEXP FilterFoldChange(SEXP fc, SEXP pct1, SEXP pct2,
  SEXP init_mask,
  double min_pct, double min_diff_pct, double logfc_threshold, 
  bool only_pos, bool not_count,
  int threads) {
  
    // ----------- copy to local
    // ---- input matrix
    bool has_init = (TYPEOF(init_mask) != NILSXP);

    bool is_matrix = Rf_isMatrix(fc);
    std::vector<double> _fc, _pct1, _pct2;
    std::vector<bool> init;
    size_t nclusters, nfeatures, nelem;
    Rcpp::StringVector features;
    if (is_matrix) {
      size_t dummyr, dummyc, dummyel;
    
      features = copy_rmatrix_to_cppvector(Rcpp::NumericMatrix(fc), _fc, nclusters, nfeatures, nelem);
      copy_rmatrix_to_cppvector(pct1, _pct1, dummyr, dummyc, dummyel);
      copy_rmatrix_to_cppvector(pct2, _pct2, dummyr, dummyc, dummyel);
      // Rprintf("FilterFoldChange matrix: nclust %ld  nfeatures %ld  nelem %ld\n", nclusters, nfeatures, nelem);
    } else {
      // Rprintf("FilterFoldChange is NOT matrix\n");
      copy_rvector_to_cppvector(Rcpp::NumericVector(fc), _fc);
      copy_rvector_to_cppvector(pct1, _pct1);
      copy_rvector_to_cppvector(pct2, _pct2);
      nelem = _fc.size();
      // Rprintf("FilterFoldChange vector nelem %ld\n", nelem);
    }

    // ---- output matrix
    std::vector<bool> mask;
    if (has_init) {
      size_t dummyr, dummyc, dummyel;
      if (is_matrix) copy_rmatrix_to_cppvector(Rcpp::LogicalMatrix(init_mask), mask, dummyr, dummyc, dummyel);
      else copy_rvector_to_cppvector(Rcpp::LogicalVector(init_mask), mask);
    } 
    else mask.assign(nelem, true);


    // ------------------------ parameter
    // int threads;
    // double min_pct, min_diff_pct, logfc_threshold;
    // bool not_count, only_pos;
    // import_filterfc_common_params(min_pct, min_diff_pct, logfc_threshold,
    //   only_pos, min_pct, min_diff_pct, logfc_threshold, only_pos);
    // import_r_common_params(not_count, threads,
    //   not_count, threads);

      // Rprintf("C++ pct %f, diff pct %f, logfc_thresh %f, pos %d\n", min_pct, min_diff_pct, logfc_threshold, (only_pos ? 1 : 0));
      // Rprintf("C++ not count %d, threads %d\n", not_count, threads);

    // ----- compute
    omp_set_num_threads(threads);
    Rprintf("THREADING: using %d threads\n", threads);

  std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double>> start = std::chrono::steady_clock::now();

    // ======= compute.
    // size_t init_count = 0, fc_count = 0, scale_count = 0, pos_count = 0, final_count = 0;
#pragma omp parallel num_threads(threads) 
// reduction( + : init_count, scale_count, fc_count, pos_count, final_count)
    {
      int tid = omp_get_thread_num();
      size_t block = nelem / threads;
      size_t rem = nelem - threads * block;
      size_t offset = tid * block + (tid > rem ? rem : tid);
      int nid = tid + 1;
      size_t end = nid * block + (nid > rem ? rem : nid);

      double mx, mn;
      bool out;
      
      for(; offset < end; ++offset) {
          out = mask[offset];
          if (!out) continue;

          // ++init_count;
          // alpha.min <- pmax(fc.results$pct.1, fc.results$pct.2)
          // features <- names(x = which(x = alpha.min >= min.pct))
          mx = std::max(_pct1[offset], _pct2[offset]);
          out &= (mx >= min_pct);
          // alpha.diff <- alpha.min - pmin(fc.results$pct.1, fc.results$pct.2)
          // x = which(x = alpha.min >= min.pct & alpha.diff >= min.diff.pct)
          mn = std::min(_pct1[offset], _pct2[offset]);
          out &= ((mx - mn) >= min_diff_pct);

          // fc_count += (out == true);

          // if (slot != "scale.data") {
          //   total.diff <- fc.results[, 1] #first column is logFC
          //   features.diff <- if (only.pos) {
          //     names(x = which(x = total.diff >= logfc.threshold))
          //   } else {
          //     names(x = which(x = abs(x = total.diff) >= logfc.threshold))
          //   }
          mx = _fc[offset];
          //   features <- intersect(x = features, y = features.diff)
          // }
          if (not_count) {
            mn = only_pos ? mx : fabs(mx);
            out &= (mn >= logfc_threshold);
            // scale_count += (out == true);
          }

          // if (only.pos)
          // de.results <- de.results[de.results[, fc.name] > 0, , drop = FALSE]
          //              | only_pos   | !only_pos 
          //  fc_ptr > 0  | y           |   y
          //  fc_ptr <= 0 | n           |   y
          if (only_pos) {
            out &= (mx > 0);
            // pos_count += (out == true);
          }
          mask[offset] = out;
          // final_count += (out == true);
      }
    }
    Rprintf("[TIME] FC Filter Elapsed(ms)= %f\n", since(start).count());
    // Rprintf("FC filter counts:  init  %ld, fc %ld, non_scale %ld, pos %ld, final %ld\n", init_count, fc_count, scale_count, pos_count, final_count);
    

  if (is_matrix) {
    Rcpp::LogicalMatrix rmask(nclusters, nfeatures, mask.begin());
    Rcpp::rownames(rmask) = Rcpp::rownames(fc);
    Rcpp::colnames(rmask) = Rcpp::colnames(fc);
    return (Rcpp::wrap(rmask));
  } else {
    Rcpp::LogicalVector rmask(mask.size());
    rmask.assign(mask.begin(), mask.end());
    return (Rcpp::wrap(rmask));
  }

}


