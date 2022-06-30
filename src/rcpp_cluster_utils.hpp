#pragma once

// ==== function declaration

#include <map>
#include <unordered_map>
#include <vector>


#ifndef NROW
#define NROW(x) INTEGER(GET_DIM((x)))[0]
#endif

#ifndef NCOL
#define NCOL(x) INTEGER(GET_DIM((x)))[1]
#endif

#define PVAL_TWO_SIDED 2
#define PVAL_GREATER 1
#define PVAL_LESS 0
#define TEST_VAL 3
#define TEST_PARAM 4


// count the number of clusters.  does not assume the labels to be in consecutive range
// produce map UNORDERED by key, so further sort is needed to get ordered key
template <typename LABEL_ITER,
  typename LABEL = typename std::iterator_traits<LABEL_ITER>::value_type >
void count_clusters(LABEL_ITER labels, size_t const & count,
  std::unordered_map<LABEL, size_t> & clust_counts);

template <typename LABEL >
void count_clusters(std::vector<LABEL> const & labels,
  std::unordered_map<LABEL, size_t> & clust_counts);

// count the number of clusters.  does not assume the labels to be in consecutive range
// produce map ORDERED by key, so no further sort is needed to get ordered key
template <typename LABEL_ITER,
  typename LABEL = typename std::iterator_traits<LABEL_ITER>::value_type >
void count_clusters(LABEL_ITER labels, size_t const & count,
  std::map<LABEL, size_t> & clust_counts);

// count the number of clusters.  does not assume the labels to be in consecutive range
// produce map ORDERED by key, so no further sort is needed to get ordered key
template <typename LABEL_ITER,
  typename LABEL = typename std::iterator_traits<LABEL_ITER>::value_type >
void count_clusters(LABEL_ITER labels, size_t const & count,
  std::vector<std::pair<LABEL, size_t> > & clust_counts);
template <typename LABEL >
void count_clusters(std::vector<LABEL> const & labels, 
  std::vector<std::pair<LABEL, size_t> > & clust_counts);

template <typename LABEL>
void get_unique_clusters(std::vector<LABEL> const & lab, 
  std::vector<LABEL> & unique_sorted);
