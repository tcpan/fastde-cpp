
#pragma once

// https://stackoverflow.com/questions/44795190/how-can-i-return-a-list-of-matrices-from-rcpp-to-r

/*
 *  similar structure as wilcox.hpp
 *
 *  Created on: Jun 8, 2021
 *  Author: Tony Pan
 *  Affiliation: Institute for Data Engineering and Science
 *  			Georgia Institute of Technology, Atlanta, GA 30332
 */

#include <iterator>
#include <vector>
#include <unordered_map>
#include <utility>

struct clust_info;

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
  bool const & exponential = false);

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
  bool const & exponential = false);

// compute percentages.  This is common amongst all variants
// this is used in seurat foldchange.default 
template <typename LABEL, typename PT_ITER, typename PT =
  typename std::iterator_traits<PT_ITER>::value_type >
void foldchange_percents(
  std::vector<std::pair<LABEL, size_t> > const & cl_counts, 
  std::unordered_map<LABEL, clust_info> const & sums, 
  size_t const & count,
  PT_ITER percent1, PT_ITER percent2);

// compute mean .  this is used for seurat "scaled.data".  no log or pseudocount applied.
//  Should not be used with expm1
template <typename LABEL, typename OT_ITER, typename OT =
  typename std::iterator_traits<OT_ITER>::value_type>
void foldchange_mean(
  std::vector<std::pair<LABEL, size_t> > const & cl_counts, 
  std::unordered_map<LABEL, clust_info> const & sums, 
  size_t const & count,
  OT_ITER out);


// compute log of mean with pseudocount.  this is used for seurat "data" and default fold change.
template <typename LABEL, typename OT_ITER, typename OT = 
  typename std::iterator_traits<OT_ITER>::value_type>
void foldchange_logmean(
  std::vector<std::pair<LABEL, size_t> > const & cl_counts, 
  std::unordered_map<LABEL, clust_info> const & sums, 
  size_t const & count,
  OT_ITER out, 
  bool const & pseudocount = true, double const & base = 2.0);


template <typename XIT, typename LIT>
void omp_dense_foldchange(
    XIT mat, size_t nsamples, size_t nfeatures,
    LIT lab,
    bool calc_percents, std::string fc_name, 
    bool use_expm1, double min_threshold, 
    bool use_log, double log_base, bool use_pseudocount, 
    std::vector<double> &fc,
    std::vector<double> &p1,
    std::vector<double> &p2,
    std::vector<std::pair<int, size_t> > &sorted_cluster_counts,
    int const & threads
);

template <typename XIT, typename IIT, typename PIT, 
    typename LIT>
void omp_sparse_foldchange(
    XIT x, IIT i, PIT p, size_t nsamples, size_t nfeatures,
    LIT lab,
    bool calc_percents, std::string fc_name, 
    bool use_expm1, double min_threshold, 
    bool use_log, double log_base, bool use_pseudocount, 
    std::vector<double> &fc,
    std::vector<double> &p1,
    std::vector<double> &p2,
    std::vector<std::pair<int, size_t> > &sorted_cluster_counts,
    int const & threads
);

template <typename VITER, typename MITER>
void omp_filter_foldchange(
  VITER fc,
  VITER pct1, 
  VITER pct2,
  MITER mask,
  size_t nelem,
  double min_pct, double min_diff_pct, 
  double logfc_threshold, 
  bool only_pos, bool not_count,
  int threads);