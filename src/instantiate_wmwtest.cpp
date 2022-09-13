#include "fastde/wmwtest.tpp"

template void sparse_sum_rank(std::pair<double, int> * temp, size_t const & nzc, 
  size_t const & count, double const & zero_val,
  std::unordered_map<int, size_t> const & z_cl_counts,
  std::unordered_map<int, size_t> & rank_sums,
  double & tie_sum
);

// ids points to start of the column's positive element row ids
// x  points to the start fo the columns positive element values
// count is hte number of positive elements in the column.
template void sparse_wmw_summary(double * in, int * ids,
  size_t const & nz_count, 
  int * labels, size_t const & count, double const & zero_val,
  std::vector<std::pair<int, size_t> > const & cl_counts,
  std::unordered_map<int, size_t> & rank_sums, double & tie_sum) ;

// rank_sums output:  map cluster to rank_sum.
// counts zeros and skip sorting of those.  this is faster than sorting all.
template void pseudosparse_wmw_summary(
  double * in, int * labels, size_t const & count, double const & zero_val,
  std::vector<std::pair<int, size_t> > const & cl_counts,
  std::unordered_map<int, size_t> & rank_sums, double & tie_sum);



// types:  0 = less, 1 = greater, 2 = twosided (default), 3 = U2
template void wmw(
  std::vector<std::pair<int, size_t> > const & cl_counts, 
  std::unordered_map<int, size_t> const & rank_sums, 
  double const & tie_sum,
  size_t const & count,
  double * out,
  int const & test_type, bool const & continuity) ;


template void omp_dense_wmw(
    double * mat, size_t const & nsamples, size_t const & nfeatures,
    int * lab, 
    int rtype, 
    bool continuity_correction, 
    std::vector<double> &pv,
    std::vector<std::pair<int, size_t> > &sorted_cluster_counts,
    int threads);


// =================================
template void omp_sparse_wmw(
    double * x, int * i, int * p, size_t nsamples, size_t nfeatures,
    int * lab,
    int rtype, 
    bool continuity_correction, 
    std::vector<double> &pv,
    std::vector<std::pair<int, size_t> > &sorted_cluster_counts,
    int threads);

template void omp_sparse_wmw(
    double * x, int * i, long * p, size_t nsamples, size_t nfeatures,
    int * lab,
    int rtype, 
    bool continuity_correction, 
    std::vector<double> &pv,
    std::vector<std::pair<int, size_t> > &sorted_cluster_counts,
    int threads);
