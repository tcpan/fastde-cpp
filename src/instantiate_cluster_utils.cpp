// ==== explicit template declaration

#include "fastde/cluster_utils.tpp"
#include "cpp11/r_vector.hpp"
#include "cpp11/integers.hpp"


template void count_clusters(int * labels, size_t const & count, std::unordered_map<int, size_t> & clust_counts, int const & threads);

template void count_clusters(std::vector<int> const & labels, std::unordered_map<int, size_t> & clust_counts, int const & threads);

template void count_clusters(int * labels, size_t const & count, std::map<int, size_t> & clust_counts, int const & threads);

template void count_clusters(int * labels, size_t const & count, std::vector<std::pair<int, size_t> > & clust_counts, int const & threads);

template void count_clusters(std::vector<int> const & labels, std::vector<std::pair<int, size_t> > & clust_counts, int const & threads);

template void count_clusters(cpp11::r_vector<int>::const_iterator labels, size_t const & count, 
    std::vector<std::pair<int, size_t> > & clust_counts, int const & threads);


template void get_unique_clusters(std::vector<int> const & labels, std::vector<int> & unique_sorted);