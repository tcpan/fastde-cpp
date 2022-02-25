#pragma once

// ==== function definition
#include <algorithm>
#include <unordered_set>
#include "rcpp_cluster_utils.hpp"


// count the number of clusters.  does not assume the labels to be in consecutive range
// produce map UNORDERED by key, so further sort is needed to get ordered key
template <typename LABEL_ITER,
  typename LABEL = typename std::iterator_traits<LABEL_ITER>::value_type >
void count_clusters(LABEL_ITER labels, size_t const & count,
  std::unordered_map<LABEL, size_t> & clust_counts) {

  clust_counts.clear();
  if (count == 0) return;

  LABEL key;
  // iterate over the input and labels.
  for (size_t i = 0; i < count; ++i) {
    key = labels[i];

    if (clust_counts.count(key) == 0) clust_counts[key] = 1;
    else {
        ++clust_counts[key];
    }
  }
}
template <typename LABEL>
void count_clusters(std::vector<LABEL> const & labels,
  std::unordered_map<LABEL, size_t> & clust_counts) {

  clust_counts.clear();

  LABEL key;
  // iterate over the input and labels.
  for (auto key : labels) {
    if (clust_counts.count(key) == 0) 
      clust_counts[key] = 1;
    else {
      ++clust_counts[key];
    }
  }
}


// count the number of clusters.  does not assume the labels to be in consecutive range
// produce map ORDERED by key, so no further sort is needed to get ordered key
template <typename LABEL_ITER,
  typename LABEL = typename std::iterator_traits<LABEL_ITER>::value_type >
void count_clusters(LABEL_ITER labels, size_t const & count,
  std::map<LABEL, size_t> & clust_counts) {

  clust_counts.clear();
  if (count == 0) return;

  // temporary, to speed up updates and avoid the O(log(n)) time for [] operator
  std::unordered_map<LABEL, size_t> temp;
  count_clusters(labels, count, temp);

  // now copy temp to final
  for (auto item : temp) {
      clust_counts[item.first] = item.second;
  }
}

// count the number of clusters.  does not assume the labels to be in consecutive range
// produce map ORDERED by key, so no further sort is needed to get ordered key
template <typename LABEL_ITER,
  typename LABEL = typename std::iterator_traits<LABEL_ITER>::value_type >
void count_clusters(LABEL_ITER labels, size_t const & count,
  std::vector<std::pair<LABEL, size_t>> & clust_counts) {

  clust_counts.clear();
  if (count == 0) return;

  // temporary, to speed up updates and avoid the O(log(n)) time for [] operator
  std::unordered_map<LABEL, size_t> temp;
  count_clusters(labels, count, temp);

  // now copy temp to final
  clust_counts.assign(temp.begin(), temp.end());
  // sort
  std::sort(clust_counts.begin(), clust_counts.end(), 
  [](std::pair<LABEL, size_t> const & left, std::pair<LABEL, size_t> const & right){
      return left.first < right.first;
  });
}

template <typename LABEL>
void get_unique_clusters(std::vector<LABEL> const & lab, 
  std::vector<LABEL> & unique_sorted) {
    std::unordered_set<LABEL> unique_labels(lab.cbegin(), lab.cend());
    unique_sorted.clear();
    unique_sorted.insert(unique_sorted.end(), unique_labels.cbegin(), unique_labels.cend());
    
    std::sort(unique_sorted.begin(), unique_sorted.end());
}
