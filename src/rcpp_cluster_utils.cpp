// ==== explicit template declaration

#include "rcpp_cluster_utils.tpp"

template void count_clusters(int * labels, size_t const & count, std::unordered_map<int, size_t> & clust_counts);

template void count_clusters(std::vector<int> const & labels, std::unordered_map<int, size_t> & clust_counts);

template void count_clusters(int * labels, size_t const & count, std::map<int, size_t> & clust_counts);

template void count_clusters(int * labels, size_t const & count, std::vector<std::pair<int, size_t>> & clust_counts);

template void get_unique_clusters(std::vector<int> const & labels, std::vector<int> & unique_sorted);