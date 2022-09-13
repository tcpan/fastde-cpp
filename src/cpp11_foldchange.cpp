#include "fastde/foldchange.hpp"

#include <chrono>

#include <cpp11/sexp.hpp>
#include <cpp11/matrix.hpp>
#include <cpp11/strings.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/doubles.hpp>

#include "utils_data.hpp"
#include "fastde/benchmark_utils.hpp"
#include "utils_sparsemat.hpp"


[[cpp11::register]]
extern cpp11::sexp cpp11_ComputeFoldChange(
  cpp11::doubles_matrix<cpp11::by_column> const & matrix, 
  cpp11::strings const & features,
  cpp11::integers const & labels,
  bool calc_percents, std::string fc_name, 
  bool use_expm1, double min_threshold, 
  bool use_log, double log_base, bool use_pseudocount, 
  bool as_dataframe,
  int threads) {


  // ----------- copy to local
  // ---- input matrix
  size_t nelem = matrix.size();
  int nsamples = matrix.nrow();
  int nfeatures = matrix.ncol();

  double * mat = reinterpret_cast<double *>(malloc(nelem * sizeof(double)));

    copy_rmatrix_to_cppvector(matrix, mat);

  // ---- label vector
  int * lab = reinterpret_cast<int *>(malloc(nsamples * sizeof(int)));
  copy_rvector_to_cppvector(labels, lab, nsamples);

  // ---- output pval matrix
  std::vector<double> fc;
  std::vector<double> p1;
  std::vector<double> p2;
  
  // Rprintf("[DEBUG] FC DN inpt size: %d %d \n", nsamples, nfeatures);


  std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double>> start = std::chrono::steady_clock::now();
  std::vector<std::pair<int, size_t> > sorted_cluster_counts;
  omp_dense_foldchange(mat, nsamples, nfeatures, lab, 
    calc_percents, fc_name, use_expm1, min_threshold,
    use_log, log_base, use_pseudocount,
    fc, p1, p2, sorted_cluster_counts, threads);

  Rprintf("[TIME] FC DN Elapsed(ms)= %f\n", since(start).count());
  // Rprintf("[DEBUG] FC DN output sizes: %d %d %d, $d\n", fc.size(), p1.size(), p2.size(), sorted_cluster_counts.size());

  free(mat);
  free(lab);
  // ------------------------ generate output
  // GET features.

  if (as_dataframe) {
    if (calc_percents)
      return(cpp11::as_sexp(export_fc_to_r_dataframe(
        fc, fc_name, 
        p1, "pct.1", 
        p2, "pct.2",
        sorted_cluster_counts, features)));
    else
      return(cpp11::as_sexp(export_vec_to_r_dataframe(fc, fc_name,
      sorted_cluster_counts, features)));
  } else {
    if (calc_percents)
      return (cpp11::as_sexp(export_fc_to_r_matrix(
        fc, fc_name, 
        p1, "pct.1", 
        p2, "pct.2",
        sorted_cluster_counts)));
    else
    // use clust for column names.
      return (cpp11::as_sexp(export_vec_to_r_matrix<cpp11::writable::doubles_matrix<cpp11::by_column>>(fc,
        sorted_cluster_counts.size(), 
        fc.size() / sorted_cluster_counts.size())));
  }
  
}


template<typename PT>
extern cpp11::sexp _compute_foldchange_sparse(
  cpp11::doubles const & _x, 
  cpp11::integers const & _i,
  cpp11::r_vector<PT> const & _p,
  cpp11::strings const & features,
  int const & rows, int const & cols,
  cpp11::integers const & labels,
  bool features_as_rows,
  bool calc_percents, std::string fc_name, 
  bool use_expm1, double min_threshold, 
  bool use_log, double log_base, bool use_pseudocount, 
  bool as_dataframe,
  int threads) {

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
  std::vector<double> fc;
  std::vector<double> p1;
  std::vector<double> p2;

  Rprintf("[TIME] FC 64 in copy Elapsed(ms)= %f\n", since(start).count());

  start = std::chrono::steady_clock::now();
  std::vector<std::pair<int, size_t> > sorted_cluster_counts;
  omp_sparse_foldchange(x, i, p, nsamples, nfeatures, lab, 
    calc_percents, fc_name, use_expm1, min_threshold, 
    use_log, log_base, use_pseudocount, 
    fc, p1, p2, sorted_cluster_counts, threads);

  Rprintf("[TIME] FC 64 Elapsed(ms)= %f\n", since(start).count());
  start = std::chrono::steady_clock::now();

  free(p);
  free(i);
  free(x);
  free(lab);

  // ------------------------ generate output
  // GET features.
  cpp11::sexp out;
  if (as_dataframe) {
    if (calc_percents)
      out = cpp11::as_sexp(export_fc_to_r_dataframe(
        fc, fc_name, 
        p1, "pct.1", 
        p2, "pct.2",
        sorted_cluster_counts, features));
    else
      out = cpp11::as_sexp(export_vec_to_r_dataframe(fc, fc_name,
      sorted_cluster_counts, features));
  } else {
    if (calc_percents)
      out = cpp11::as_sexp(export_fc_to_r_matrix(
        fc, fc_name, 
        p1, "pct.1", 
        p2, "pct.2",
        sorted_cluster_counts));
    else
    // use clust for column names.
      out = cpp11::as_sexp(export_vec_to_r_matrix<cpp11::writable::doubles_matrix<cpp11::by_column>>(fc, 
        sorted_cluster_counts.size(), 
        fc.size() / sorted_cluster_counts.size()));
  }
  Rprintf("[TIME] FC 64 out wrap Elapsed(ms)= %f\n", since(start).count());
  return out;

}



[[cpp11::register]]
extern cpp11::sexp cpp11_ComputeFoldChangeSparse(
  cpp11::doubles const & x, 
  cpp11::integers const & i,
  cpp11::integers const & p,
  cpp11::strings const & features,
  int const & rows, int const & cols,
  cpp11::integers const & labels,
  bool features_as_rows,
  bool calc_percents, std::string fc_name, 
  bool use_expm1, double min_threshold, 
  bool use_log, double log_base, bool use_pseudocount, 
  bool as_dataframe,
  int threads) {
    return _compute_foldchange_sparse(x, i, p, features, rows, cols, 
      labels, features_as_rows, calc_percents, fc_name, use_expm1, min_threshold, 
      use_log, log_base, use_pseudocount, as_dataframe, threads);
}



[[cpp11::register]]
extern cpp11::sexp cpp11_ComputeFoldChangeSparse64(
  cpp11::doubles const & x, 
  cpp11::integers const & i,
  cpp11::doubles const & p,
  cpp11::strings const & features,
  int const & rows, int const & cols,
  cpp11::integers const & labels,
  bool features_as_rows,
  bool calc_percents, std::string fc_name, 
  bool use_expm1, double min_threshold, 
  bool use_log, double log_base, bool use_pseudocount, 
  bool as_dataframe,
  int threads) {
    return _compute_foldchange_sparse(x, i, p, features, rows, cols, 
      labels, features_as_rows, calc_percents, fc_name, use_expm1, min_threshold, 
      use_log, log_base, use_pseudocount, as_dataframe, threads);
  }



[[cpp11::register]]
extern cpp11::writable::logicals cpp11_FilterFoldChange(
  cpp11::doubles const & fc,
  cpp11::doubles const & pct1, 
  cpp11::doubles const & pct2,
  cpp11::logicals const & init_mask,
  double min_pct, double min_diff_pct, double logfc_threshold, 
  bool only_pos, bool not_count,
  int threads) {
  
    // ----------- copy to local
    // ---- input matrix

    size_t nelem = fc.size();

  double * _fc = reinterpret_cast<double *>(malloc(nelem * sizeof(double)));
  double * _pct1 = reinterpret_cast<double *>(malloc(nelem * sizeof(double)));
  double * _pct2 = reinterpret_cast<double *>(malloc(nelem * sizeof(double)));

    copy_rvector_to_cppvector(fc, _fc);
    copy_rvector_to_cppvector(pct1, _pct1);
    copy_rvector_to_cppvector(pct2, _pct2);

    // ---- output matrix
  bool * mask = reinterpret_cast<bool *>(malloc(nelem * sizeof(bool)));
    if (init_mask.size() >= nelem) {
      copy_rvector_to_cppvector(init_mask, mask, nelem);
    } 
    else std::fill(mask, mask + nelem, true);

    // ------------------------ parameter
    // int threads;
    // double min_pct, min_diff_pct, logfc_threshold;
    // bool not_count, only_pos;
    // import_filterfc_common_params(min_pct, min_diff_pct, logfc_threshold,
    //   only_pos, min_pct, min_diff_pct, logfc_threshold, only_pos);
    // import_r_common_params(not_count, threads,
    //   not_count, threads);

      // Rprintf("C++ pct %f, diff pct %f, logfc_thresh %f, pos %d\n", min_pct, min_diff_pct, logfc_threshold, (only_pos ? 1 : 0));
      // Rprintf("C++ not count %d, threads %d\n", not_count, threads);

    // ----- compute

  std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double>> start = std::chrono::steady_clock::now();

  omp_filter_foldchange(_fc, _pct1, _pct2, mask, nelem,
    min_pct, min_diff_pct, logfc_threshold, only_pos, not_count, threads);

    Rprintf("[TIME] FC Filter Elapsed(ms)= %f\n", since(start).count());
    // Rprintf("FC filter counts:  init  %ld, fc %ld, non_scale %ld, pos %ld, final %ld\n", init_count, fc_count, scale_count, pos_count, final_count);
    
    cpp11::writable::logicals out = export_vec_to_rvec<cpp11::writable::logicals>(mask, nelem);

    free(_fc);
    free(_pct1);
    free(_pct2);
    free(mask);

    return out;
}



[[cpp11::register]]
extern cpp11::writable::logicals_matrix<cpp11::by_column> cpp11_FilterFoldChangeMat(
  cpp11::doubles_matrix<cpp11::by_column> const & fc,
  cpp11::doubles_matrix<cpp11::by_column> const & pct1, 
  cpp11::doubles_matrix<cpp11::by_column> const & pct2,
  cpp11::logicals_matrix<cpp11::by_column> const & init_mask,
  double min_pct, double min_diff_pct, double logfc_threshold, 
  bool only_pos, bool not_count,
  int threads) {
  
    // ----------- copy to local
    // ---- input matrix

  size_t nelem = fc.size();
  int nclusters = fc.nrow();
  int nfeatures = fc.ncol();

  double * _fc = reinterpret_cast<double *>(malloc(nelem * sizeof(double)));
  double * _pct1 = reinterpret_cast<double *>(malloc(nelem * sizeof(double)));
  double * _pct2 = reinterpret_cast<double *>(malloc(nelem * sizeof(double)));

    copy_rmatrix_to_cppvector(fc, _fc);
    copy_rmatrix_to_cppvector(pct1, _pct1);
    copy_rmatrix_to_cppvector(pct2, _pct2);

    // Rprintf("FilterFoldChange matrix: nclust %ld  nfeatures %ld  nelem %ld\n", nclusters, nfeatures, nelem);

    // ---- output matrix
  bool * mask = reinterpret_cast<bool *>(malloc(nelem * sizeof(bool)));  
    if ((init_mask.nrow() >= nclusters) || (init_mask.ncol() >= nfeatures)) {
      copy_rmatrix_to_cppvector(init_mask, mask);
    } 
    else std::fill(mask, mask + nelem, true);


    // ------------------------ parameter
    // int threads;
    // double min_pct, min_diff_pct, logfc_threshold;
    // bool not_count, only_pos;
    // import_filterfc_common_params(min_pct, min_diff_pct, logfc_threshold,
    //   only_pos, min_pct, min_diff_pct, logfc_threshold, only_pos);
    // import_r_common_params(not_count, threads,
    //   not_count, threads);

      // Rprintf("C++ pct %f, diff pct %f, logfc_thresh %f, pos %d\n", min_pct, min_diff_pct, logfc_threshold, (only_pos ? 1 : 0));
      // Rprintf("C++ not count %d, threads %d\n", not_count, threads);

    // ----- compute

  std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double>> start = std::chrono::steady_clock::now();

  omp_filter_foldchange(_fc, _pct1, _pct2, mask, nelem,
    min_pct, min_diff_pct, logfc_threshold, only_pos, not_count, threads);

    Rprintf("[TIME] FC Filter Elapsed(ms)= %f\n", since(start).count());
    // Rprintf("FC filter counts:  init  %ld, fc %ld, non_scale %ld, pos %ld, final %ld\n", init_count, fc_count, scale_count, pos_count, final_count);
    
  auto out = export_vec_to_r_matrix<cpp11::writable::logicals_matrix<cpp11::by_column>>(mask, nclusters, nfeatures);
    free(_fc);
    free(_pct1);
    free(_pct2);
    free(mask);

  return out;

}


