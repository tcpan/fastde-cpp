#include "fastde/ttest.tpp"

// sum, sum of square, and 
template struct gaussian_stats<double>;


// ids points to start of the column's positive element row ids
// x  points to the start fo the columns positive element values
// count is hte number of positive elements in the column.
template void sparse_ttest_summary(double * in, int * ids, size_t const & nz_count, 
  int * labels, size_t const & count, double const & zero_val,
  std::vector<std::pair<int, size_t> > const & clust_counts,
  std::unordered_map<int, gaussian_stats<double> > & gaussian_sums);

// compute count, mean and variance for each subclass.
template void dense_ttest_summary(
  double * in, int * labels, size_t const & count, double const & zero_val,
  std::vector<std::pair<int, size_t> > const & clust_counts,
  std::unordered_map<int, gaussian_stats<double> > & gaussian_sums);

template class incomplete_beta<double, double, true>;
template class incomplete_beta<double, double, false>;


template class t_distribution<double, double, true>;
template class t_distribution<double, double, false>;


// types:  
template void two_sample_ttest(
  std::unordered_map<int, gaussian_stats<double> > const & gaussian_sums,
  std::vector<std::pair<int, size_t> > const & clust_counts,
  double * out, 
  int const & test_type, 
  bool const equal_variance) ;

// omp t-Test for dense matrix
template void omp_dense_ttest(
    double * mat, size_t const & nsamples, size_t const & nfeatures,
    int * lab, 
    int alternative, 
    bool var_equal, 
    std::vector<double> &pv,
    std::vector<std::pair<int, size_t> > &sorted_cluster_counts,
    int threads) ;

// =================================
template void omp_sparse_ttest(
    double * x, int * i, int * p, size_t nsamples, size_t nfeatures,
    int * lab,
    int alternative, 
    bool var_equal, 
    std::vector<double> &pv,
    std::vector<std::pair<int, size_t> > &sorted_cluster_counts,
    int threads);
template void omp_sparse_ttest(
    double * x, int * i, long * p, size_t nsamples, size_t nfeatures,
    int * lab,
    int alternative, 
    bool var_equal, 
    std::vector<double> &pv,
    std::vector<std::pair<int, size_t> > &sorted_cluster_counts,
    int threads);
