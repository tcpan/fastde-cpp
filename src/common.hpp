#pragma once


#include <map>
#include <unordered_map>
#include <vector>


#ifndef NROW
#define NROW(x) INTEGER(GET_DIM((x)))[0]
#endif

#ifndef NCOL
#define NCOL(x) INTEGER(GET_DIM((x)))[1]
#endif



// count the number of clusters.  does not assume the labels to be in consecutive range
// produce map UNORDERED by key, so further sort is needed to get ordered key
template <typename LABEL>
void count_clusters(LABEL const * labels, size_t const & count,
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


// count the number of clusters.  does not assume the labels to be in consecutive range
// produce map ORDERED by key, so no further sort is needed to get ordered key
template <typename LABEL>
void count_clusters(LABEL const * labels, size_t const & count,
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
template <typename LABEL>
void count_clusters(LABEL const * labels, size_t const & count,
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


