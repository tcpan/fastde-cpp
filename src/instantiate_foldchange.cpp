
// https://stackoverflow.com/questions/44795190/how-can-i-return-a-list-of-matrices-from-rcpp-to-r

/*
 *  similar structure as wilcox.hpp
 *
 *  Created on: Jun 8, 2021
 *  Author: Tony Pan
 *  Affiliation: Institute for Data Engineering and Science
 *  			Georgia Institute of Technology, Atlanta, GA 30332
 */

#include "fastde/foldchange.tpp"

#include <vector>
#include <unordered_map>
#include <utility>


struct clust_info;

template void sparse_foldchange_summary(
  double * in, int * row_ids, size_t const & nz_count,
  int * labels, size_t const & count, double const & zero_val,
  std::vector<std::pair<int, size_t> > const & cl_counts,
  std::unordered_map<int, clust_info> & sums,
  double const & min_thresh, 
  bool const & exponential);

// compute base info for foldchange.  
// if expm1 is used, then it's suitable to replace seurat "data" fold change
// else it is used for "scaled.data" or default action.
// NOTE this has an extra entry for the total for the whole "in" array.
template void pseudosparse_foldchange_summary(
  double * in, int * labels, size_t const & count, double const & zero_val,
  std::vector<std::pair<int, size_t> > const & cl_counts,
  std::unordered_map<int, clust_info> & sums,
  double const & min_thresh, 
  bool const & exponential);

// compute percentages.  This is common amongst all variants
// this is used in seurat foldchange.default 
template void foldchange_percents(
  std::vector<std::pair<int, size_t> > const & cl_counts, 
  std::unordered_map<int, clust_info> const & sums, 
  size_t const & count,
  double * percent1, double * percent2);

// compute mean .  this is used for seurat "scaled.data".  no log or pseudocount applied.
//  Should not be used with expm1
template void foldchange_mean(
  std::vector<std::pair<int, size_t> > const & cl_counts, 
  std::unordered_map<int, clust_info> const & sums, 
  size_t const & count,
  double * out);


// compute log of mean with pseudocount.  this is used for seurat "data" and default fold change.
template void foldchange_logmean(
  std::vector<std::pair<int, size_t> > const & cl_counts, 
  std::unordered_map<int, clust_info> const & sums, 
  size_t const & count,
  double * out, 
  bool const & pseudocount, double const & base);


template void omp_dense_foldchange(
    double * mat, size_t nsamples, size_t nfeatures,
    int * lab,
    bool calc_percents, std::string fc_name, 
    bool use_expm1, double min_threshold, 
    bool use_log, double log_base, bool use_pseudocount, 
    std::vector<double> &fc,
    std::vector<double> &p1,
    std::vector<double> &p2,
    std::vector<std::pair<int, size_t> > &sorted_cluster_counts,
    int const & threads
);

template void omp_sparse_foldchange(
    double * x, int * i, int * p, size_t nsamples, size_t nfeatures,
    int * lab,
    bool calc_percents, std::string fc_name, 
    bool use_expm1, double min_threshold, 
    bool use_log, double log_base, bool use_pseudocount, 
    std::vector<double> &fc,
    std::vector<double> &p1,
    std::vector<double> &p2,
    std::vector<std::pair<int, size_t> > &sorted_cluster_counts,
    int const & threads
);


template void omp_sparse_foldchange(
    double * x, int * i, long * p, size_t nsamples, size_t nfeatures,
    int * lab,
    bool calc_percents, std::string fc_name, 
    bool use_expm1, double min_threshold, 
    bool use_log, double log_base, bool use_pseudocount, 
    std::vector<double> &fc,
    std::vector<double> &p1,
    std::vector<double> &p2,
    std::vector<std::pair<int, size_t> > &sorted_cluster_counts,
    int const & threads
);

template void omp_filter_foldchange(
  double * fc,
  double * pct1, 
  double * pct2,
  bool * mask,
  size_t nelem,
  double min_pct, double min_diff_pct, 
  double logfc_threshold, 
  bool only_pos, bool not_count,
  int threads);



