#include "fastde/wmwtest.hpp"

#include <cpp11/sexp.hpp>
#include <cpp11/matrix.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/strings.hpp>
#include <cpp11/doubles.hpp>

#include "fastde/benchmark_utils.hpp"
#include "utils_data.hpp"
#include "utils_sparsemat.hpp"


[[cpp11::register]]
extern cpp11::sexp cpp11_dense_wmw(
    cpp11::doubles_matrix<cpp11::by_column> const & input, 
    cpp11::strings const & features,
    cpp11::integers const & labels, 
    int rtype, 
    bool continuity_correction, 
    bool as_dataframe,
    int threads) {

  std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double>> start;
  start = std::chrono::steady_clock::now();
  // Rprintf("[TIME] WMW DN start Elapsed(ms)= %f\n", since(start).count());
  start = std::chrono::steady_clock::now();

  // ----------- copy to local
  // ---- input matrix
  int nsamples = input.nrow();
  int nfeatures = input.ncol();
  size_t nelem = nsamples * nfeatures;

  double * mat = reinterpret_cast<double *>(malloc(nelem * sizeof(double)));
  // Rprintf("malloc for mat : %dx%d %ld :  %x\n", nsamples, nfeatures, nelem, mat);

  copy_rmatrix_to_cppvector(input, mat);

  // ---- label vector
  int * lab = reinterpret_cast<int *>(malloc(nsamples * sizeof(int)));
  // Rprintf("malloc for lab : %d :  %x\n", nsamples, lab);
  copy_rvector_to_cppvector(labels, lab, nsamples);

  // ---- output pval matrix
  std::vector<double> pv;

  Rprintf("[TIME] WMW DN in copy Elapsed(ms)= %f\n", since(start).count());


  start = std::chrono::steady_clock::now();

  std::vector<std::pair<int, size_t> > sorted_cluster_counts;
  omp_dense_wmw(mat, nsamples, nfeatures, lab, 
    rtype, continuity_correction, pv, sorted_cluster_counts, threads);

  Rprintf("[TIME] WMW DN Elapsed(ms)= %f\n", since(start).count());

  free(mat);
  free(lab);

  // ------------------------ generate output
  if (as_dataframe) {
    return(cpp11::as_sexp(export_vec_to_r_dataframe(pv, "p_val", sorted_cluster_counts, features)));
  } else {
    // use clust for column names.
    return (cpp11::as_sexp(export_vec_to_r_matrix<cpp11::writable::doubles_matrix<cpp11::by_column>>(pv,
      sorted_cluster_counts.size(), pv.size() / sorted_cluster_counts.size())));
  }
}


template <typename PT>
extern cpp11::sexp _compute_wmwtest_sparse(
    cpp11::doubles const & _x, 
    cpp11::integers const & _i,
    cpp11::r_vector<PT> const & _p,
    cpp11::strings const & features,
    int const & rows, int const & cols,
    cpp11::integers const & labels,
    bool features_as_rows,
    int rtype, 
    bool continuity_correction, 
    bool as_dataframe,
    int threads) {
  // Rprintf("here 1\n");

    using PT2 = typename std::conditional<std::is_same<PT, double>::value, long, int>::type;

  std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double>> start;

  start = std::chrono::steady_clock::now();
  // ----------- copy to local
  // ---- input matrix
  PT2 nelem = _x.size();
  int nsamples = rows;
  int nfeatures = cols;

  double * x = reinterpret_cast<double *>(malloc(nelem * sizeof(double)));
  int * i = reinterpret_cast<int *>(malloc(nelem * sizeof(int)));
  PT2 * p;

  if (features_as_rows) {
    p = reinterpret_cast<PT2 *>(malloc((rows+1) * sizeof(PT2)));;
    // transpose
    if (threads == 1) {
      _sp_transpose(_x, _i, _p, rows, cols, x, i, p, threads);
    } else {
      _sp_transpose_par(_x, _i, _p, rows, cols, x, i, p, threads);
    }
    nsamples = cols;
    nfeatures = rows;
  } else {
    p = reinterpret_cast<PT2 *>(malloc((cols+1) * sizeof(PT2)));
    copy_rvector_to_cppvector(_x, x);
    copy_rvector_to_cppvector(_i, i);
    copy_rvector_to_cppvector(_p, p);
  }
  // Rprintf("Sparse DIM: samples %lu x features %lu, non-zeros %lu\n", nsamples, nfeatures, nelem); 

  // ---- label vector
  int * lab = reinterpret_cast<int *>(malloc(nsamples * sizeof(int)));
  copy_rvector_to_cppvector(labels, lab, nsamples);

  // ---- output pval matrix
  std::vector<double> pv;

  Rprintf("[TIME] WMW in copy Elapsed(ms)= %f\n", since(start).count());


  start = std::chrono::steady_clock::now();
  std::vector<std::pair<int, size_t> > sorted_cluster_counts;
  omp_sparse_wmw(x, i, p, nsamples, nfeatures, lab, 
    rtype, continuity_correction, 
    pv, sorted_cluster_counts, threads);

  Rprintf("[TIME] WMW Elapsed(ms)= %f\n", since(start).count());

  free(p);
  free(i);
  free(x);
  free(lab);
  // ------------------------ generate output
  if (as_dataframe) {
    return(cpp11::as_sexp(export_vec_to_r_dataframe(pv, "p_val", sorted_cluster_counts, features)));
  } else {
    // use clust for column names.
    return (cpp11::as_sexp(export_vec_to_r_matrix<cpp11::writable::doubles_matrix<cpp11::by_column>>(pv,
      sorted_cluster_counts.size(), pv.size() / sorted_cluster_counts.size())));
  }
}




[[cpp11::register]]
extern cpp11::sexp cpp11_sparse_wmw(
    cpp11::doubles const & x, 
    cpp11::integers const & i,
    cpp11::integers const & p,
    cpp11::strings const & features,
    int const & rows, int const & cols,
    cpp11::integers const & labels,
    bool features_as_rows,
    int rtype, 
    bool continuity_correction, 
    bool as_dataframe,
    int threads) {

    return _compute_wmwtest_sparse(x, i, p, features, rows, cols,
      labels, features_as_rows, rtype, continuity_correction, as_dataframe, threads);

}


[[cpp11::register]]
extern cpp11::sexp cpp11_sparse64_wmw(
    cpp11::doubles const & x, 
    cpp11::integers const & i,
    cpp11::doubles const & p,
    cpp11::strings const & features,
    int const & rows, int const & cols,
    cpp11::integers const & labels,
    bool features_as_rows,
    int rtype, 
    bool continuity_correction, 
    bool as_dataframe,
    int threads) {

    return _compute_wmwtest_sparse(x, i, p, features, rows, cols,
      labels, features_as_rows,  rtype, continuity_correction,  as_dataframe, threads);

}
