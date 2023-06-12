#pragma once

#include "fastde/cluster_utils.hpp"

// ==== function definition
#include <algorithm>
#include <unordered_set>
#include <omp.h>

#include "fastde/benchmark_utils.hpp"


// count the number of clusters.  does not assume the labels to be in consecutive range
// produce map UNORDERED by key, so further sort is needed to get ordered key
template <typename LABEL_ITER, typename LABEL>
void count_clusters(LABEL_ITER labels, size_t const & count,
  std::unordered_map<LABEL, size_t> & clust_counts, int const & threads) {

  clust_counts.clear();

  if (threads == 1) {
    LABEL key;

    // iterate over the input and labels.
    for (size_t i = 0; i < count; ++i) {
      key = *(labels + i);

      if (clust_counts.count(key) == 0) clust_counts[key] = 1;
      else ++clust_counts[key];
    }
  } else {
    std::vector<std::unordered_map<LABEL, size_t>> temp(threads);
    
#pragma omp parallel num_threads(threads)
{
    int tid = omp_get_thread_num();
    size_t block = count / threads;
    size_t rem = count - threads * block;
    size_t offset = tid * block + ((tid > static_cast<int>(rem)) ? rem : tid);
    int nid = tid + 1;
    size_t end = nid * block + ((nid > static_cast<int>(rem)) ? rem : nid);  

    LABEL key;
    // count per thread.
    for (; offset < end; ++offset) {
      key = *(labels + offset);

      if (temp[tid].count(key) == 0) temp[tid][key] = 1;
      else ++temp[tid][key];
    }
}
    // now summarize
    LABEL key;
    for (int t = 0; t < threads; ++t) {
      for (auto item : temp[t]) {
        key = item.first;
        if (clust_counts.count(key) == 0) clust_counts[key] = item.second;
        else clust_counts[key] += item.second;
      }
    }

  } // end parallel version.
}

template <typename LABEL_VEC, typename LABEL>
void count_clusters_vec(LABEL_VEC const & labels, size_t const & count,
  std::unordered_map<LABEL, size_t> & clust_counts, int const & threads){


  clust_counts.clear();

  if (threads == 1) {
    LABEL key;

    // iterate over the input and labels.
    for (size_t i = 0; i < count; ++i) {
      key = labels[i];

      if (clust_counts.count(key) == 0) clust_counts[key] = 1;
      else ++clust_counts[key];
    }
  } else {
    std::vector<std::unordered_map<LABEL, size_t>> temp(threads);
    
#pragma omp parallel num_threads(threads)
{
    int tid = omp_get_thread_num();
    size_t block = count / threads;
    size_t rem = count - threads * block;
    size_t offset = tid * block + ((tid > static_cast<int>(rem)) ? rem : tid);
    int nid = tid + 1;
    size_t end = nid * block + ((nid > static_cast<int>(rem)) ? rem : nid);  

    LABEL key;
    // count per thread.
    for (; offset < end; ++offset) {
      key = labels[offset];

      if (temp[tid].count(key) == 0) temp[tid][key] = 1;
      else ++temp[tid][key];
    }
}
    // now summarize
    LABEL key;
    for (int t = 0; t < threads; ++t) {
      for (auto item : temp[t]) {
        key = item.first;
        if (clust_counts.count(key) == 0) clust_counts[key] = item.second;
        else clust_counts[key] += item.second;
      }
    }

  } // end parallel version.

}

template <typename LABEL>
void count_clusters(std::vector<LABEL> const & labels,
  std::unordered_map<LABEL, size_t> & clust_counts, int const & threads) {
  count_clusters(labels.begin(), labels.size(), clust_counts, threads);
}


// count the number of clusters.  does not assume the labels to be in consecutive range
// produce map ORDERED by key, so no further sort is needed to get ordered key
template <typename LABEL_ITER, typename LABEL>
void count_clusters(LABEL_ITER labels, size_t const & count,
  std::map<LABEL, size_t> & clust_counts, int const & threads) {
  std::unordered_map<LABEL, size_t> temp;
  count_clusters(labels, count, temp, threads);
  // now copy temp to final
  for (auto item : temp) {
      clust_counts[item.first] = item.second;
  }
}

// count the number of clusters.  does not assume the labels to be in consecutive range
// produce map ORDERED by key, so no further sort is needed to get ordered key
template <typename LABEL_ITER,
  typename LABEL >
void count_clusters(LABEL_ITER labels, size_t const & count,
  std::vector<std::pair<LABEL, size_t> > & clust_counts, int const & threads) {
  std::unordered_map<LABEL, size_t> temp;
  count_clusters(labels, count, temp, threads);
  // now copy temp to final
  clust_counts.assign(temp.begin(), temp.end());
  // sort
  auto sorter = [](std::pair<LABEL, size_t> const & left, std::pair<LABEL, size_t> const & right){
      return left.first < right.first;
  };
  std::sort(clust_counts.begin(), clust_counts.end(), sorter);
}

template <typename LABEL_VEC, typename LABEL>
void count_clusters_vec(LABEL_VEC const & labels, size_t const & count,
  std::vector<std::pair<LABEL, size_t> > & clust_counts, int const & threads) {

  std::unordered_map<LABEL, size_t> temp;
  count_clusters_vec(labels, count, temp, threads);
  // now copy temp to final
  clust_counts.assign(temp.begin(), temp.end());
  // sort
  auto sorter = [](std::pair<LABEL, size_t> const & left, std::pair<LABEL, size_t> const & right){
      return left.first < right.first;
  };
  std::sort(clust_counts.begin(), clust_counts.end(), sorter);

}


template <typename LABEL >
void count_clusters(std::vector<LABEL> const & labels,
  std::vector<std::pair<LABEL, size_t> > & clust_counts, int const & threads) {
    count_clusters(labels.begin(), labels.size(), clust_counts, threads);
}

template <typename LABEL>
void get_unique_clusters(std::vector<LABEL> const & lab, 
  std::vector<LABEL> & unique_sorted) {
    std::unordered_set<LABEL> unique_labels(lab.cbegin(), lab.cend());
    unique_sorted.assign(unique_labels.cbegin(), unique_labels.cend());    
    std::sort(unique_sorted.begin(), unique_sorted.end());
}
