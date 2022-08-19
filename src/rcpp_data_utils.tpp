#pragma once

// ------- function definition

#include "rcpp_data_utils.hpp"

#include <algorithm>

// ------- class def
namespace Rcpp {

    // specialization of Rcpp::as
    template <> dgCMatrix as(SEXP mat) { return dgCMatrix(mat); };

    // specialization of Rcpp::wrap
    template <> SEXP wrap(const dgCMatrix& sm) {
        return sm.get_S4();
    };

    // specialization of Rcpp::as 
    template <> dgCMatrix64 as(SEXP mat) { return dgCMatrix64(mat); };

    template <> SEXP wrap(const dgCMatrix64& sm) {
        return sm.get_S4();
    };

    // specialization of Rcpp::as
    template <> dgMatrix as(SEXP mat) { return dgMatrix(mat); };

    // specialization of Rcpp::wrap
    template <> SEXP wrap(const dgMatrix& sm) {
        return sm.get_SEXP();
    };


}

Rcpp::dgCMatrix rttest_dgCMatrix(Rcpp::dgCMatrix const & mat){
    return mat;
}

Rcpp::dgCMatrix64 rttest_dgCMatrix64(Rcpp::dgCMatrix64 const & mat){
    return mat;
}


// ------- function def


Rcpp::StringVector copy_rmatrix_to_cppvector(Rcpp::NumericMatrix const & _matrix, std::vector<double> & mat,
    size_t & nrow, size_t & ncol, size_t & nelem) {
    nrow=_matrix.nrow(); // n
    ncol=_matrix.ncol(); // m
    nelem = nrow * ncol;

    mat.clear();
    mat.insert(mat.end(), _matrix.cbegin(), _matrix.cend());

    return Rcpp::colnames(_matrix);
}
Rcpp::StringVector copy_rmatrix_to_cppvector(Rcpp::LogicalMatrix const & _matrix, std::vector<bool> & mat,
    size_t & nrow, size_t & ncol, size_t & nelem) {
    nrow=_matrix.nrow(); // n
    ncol=_matrix.ncol(); // m
    nelem = nrow * ncol;

    mat.clear();
    mat.insert(mat.end(), _matrix.cbegin(), _matrix.cend());

    return Rcpp::colnames(_matrix);
}


size_t copy_rvector_to_cppvector(Rcpp::LogicalVector const & _vector, std::vector<bool> & vec, size_t const & length, size_t const & offset) {
    size_t veclen = static_cast<size_t>(_vector.size());
    if (offset >= veclen) return 0;  // offset is greater than vector length
    size_t len = std::min(length, veclen - offset);   // length to return
    vec.clear();
    vec.insert(vec.end(), _vector.cbegin() + offset, _vector.cbegin() + offset + len);

    return len;
}

size_t copy_rvector_to_cppvector(Rcpp::NumericVector const & _vector, std::vector<double> & vec, size_t const & length, size_t const & offset) {
    size_t veclen = static_cast<size_t>(_vector.size());
    if (offset >= veclen) return 0;  // offset is greater than vector length
    size_t len = std::min(length, veclen - offset);   // length to return
    vec.clear();
    vec.insert(vec.end(), _vector.cbegin() + offset, _vector.cbegin() + offset + len);

    return len;
}

size_t copy_rvector_to_cppvector(Rcpp::IntegerVector const & _vector, std::vector<int> & vec, size_t const & length, size_t const & offset) {
    size_t veclen = static_cast<size_t>(_vector.size());
    if (offset >= veclen) return 0;  // offset is greater than vector length
    size_t len = std::min(length, veclen - offset);   // length to return
    vec.clear();
    vec.insert(vec.end(), _vector.cbegin() + offset, _vector.cbegin() + offset + len);

    return len;
}
size_t copy_rvector_to_cppvector(Rcpp::NumericVector const & _vector, std::vector<long> & vec, size_t const & length, size_t const & offset) {
    size_t veclen = static_cast<size_t>(_vector.size());
    if (offset >= veclen) return 0;  // offset is greater than vector length
    size_t len = std::min(length, veclen - offset);   // length to return
    vec.clear();
    vec.insert(vec.end(), _vector.cbegin() + offset, _vector.cbegin() + offset + len);

    return len;
}

size_t copy_rvector_to_cppvector(SEXP _vector, std::vector<double> & vec, size_t const & length, size_t const & offset) {
    if (TYPEOF(_vector) != REALSXP) return 0;

    size_t veclen = Rf_xlength(_vector);
    if (offset >= veclen) return 0;  // offset is greater than vector length
    size_t len = std::min(length, veclen - offset);   // length to return
    
    double * val = REAL(_vector);
    val += offset;

    vec.clear();
    vec.insert(vec.end(), val, val + len);

    return len;
}

size_t copy_rvector_to_cppvector(SEXP _vector, std::vector<bool> & vec, size_t const & length, size_t const & offset) {
    if (TYPEOF(_vector) != LGLSXP) return 0;

    size_t veclen = Rf_xlength(_vector);
    if (offset >= veclen) return 0;  // offset is greater than vector length
    size_t len = std::min(length, veclen - offset);   // length to return
    
    int * val = LOGICAL(_vector);
    val += offset;

    vec.clear();
    vec.insert(vec.end(), val, val + len);

    return len;
}

size_t copy_rvector_to_cppvector(SEXP _vector, std::vector<int> & vec, size_t const & length, size_t const & offset) {
    if (TYPEOF(_vector) != INTSXP) return 0;

    size_t veclen = Rf_xlength(_vector);
    if (offset >= veclen) return 0;  // offset is greater than vector length
    size_t len = std::min(length, veclen - offset);   // length to return
    
    int * val = INTEGER(_vector);
    val += offset;

    vec.clear();
    vec.insert(vec.end(), val, val + len);

    return len;
}

size_t copy_rvector_to_cppvector(SEXP _vector, std::vector<long> & vec, size_t const & length, size_t const & offset) {
    if (TYPEOF(_vector) != REALSXP) return 0;

    size_t veclen = Rf_xlength(_vector);
    if (offset >= veclen) return 0;  // offset is greater than vector length
    size_t len = std::min(length, veclen - offset);   // length to return
    
    double * val = REAL(_vector);
    val += offset;


    vec.clear();
    vec.insert(vec.end(), val, val + len);

    return len;
}



Rcpp::StringVector copy_rsparsematrix_to_cppvectors(Rcpp::dgCMatrix const & obj, 
    std::vector<double> & x,
    std::vector<int> & i,
    std::vector<int> & p,
    size_t & nrow, size_t & ncol, size_t & nelem) {

    copy_rvector_to_cppvector(obj.get_colindices_SEXP(), i);
    copy_rvector_to_cppvector(obj.get_rowpointers_SEXP(), p);
    copy_rvector_to_cppvector(obj.get_entries_SEXP(), x);
    
    nrow = obj.get_nrow();   // sample count,  nrow
    ncol = obj.get_ncol();   // feature/gene count,  ncol
    nelem = p[ncol];   // since p is offsets, the ncol+1 entry has the total count

    // get the column names == features.
    // Rcpp::List dimnms = obj.Dimnames;
    // SEXP rownms = VECTOR_ELT(dimnms, 0);    // samples
    // Rcpp::StringVector features = dimnms[1];   // features_names = columnames
    // GET features.
    return Rcpp::as<Rcpp::StringVector>(obj.get_colnames());  // 
}

void copy_rsparsematrix_to_cppvectors(
    Rcpp::NumericVector const & _x, Rcpp::NumericVector const & _i, Rcpp::NumericVector const & _p, 
    std::vector<double> & x,
    std::vector<long> & i,
    std::vector<long> & p) {

    copy_rvector_to_cppvector(_i, i);
    copy_rvector_to_cppvector(_p, p);
    copy_rvector_to_cppvector(_x, x);
}

void copy_rsparsematrix_to_cppvectors(
    Rcpp::NumericVector const & _x, Rcpp::IntegerVector const & _i, Rcpp::IntegerVector const & _p, 
    std::vector<double> & x,
    std::vector<int> & i,
    std::vector<int> & p) {

    copy_rvector_to_cppvector(_i, i);
    copy_rvector_to_cppvector(_p, p);
    copy_rvector_to_cppvector(_x, x);
}




Rcpp::StringVector copy_rsparsematrix_to_cppvectors(Rcpp::dgCMatrix64 const & obj, 
    std::vector<double> & x,
    std::vector<long> & i,
    std::vector<long> & p,
    size_t & nrow, size_t & ncol, size_t & nelem) {

    copy_rsparsematrix_to_cppvectors(
        obj.get_entries_SEXP(), obj.get_colindices_SEXP(), obj.get_rowpointers_SEXP(), 
        x, i, p
    );

    nrow = obj.get_nrow();   // sample count,  nrow
    ncol = obj.get_ncol();   // feature/gene count,  ncol
    nelem = p[ncol];   // since p is offsets, the ncol+1 entry has the total count

    // Rcpp::List dimnms = obj.Dimnames;
    // SEXP rownms = VECTOR_ELT(dimnms, 0);    // samples
    // Rcpp::StringVector features = dimnms[1];   // features_names = columnames
    // GET features.
    return Rcpp::as<Rcpp::StringVector>(obj.get_colnames());  // 
}



SEXP copy_rsparsematrix_to_cppvectors(SEXP obj, 
    std::vector<double> & x,
    std::vector<int> & i,
    std::vector<int> & p,
    size_t & nrow, size_t & ncol, size_t & nelem) {
    
    if (TYPEOF(obj) != S4SXP) return NILSXP;

    copy_rvector_to_cppvector(R_do_slot(obj, Rf_install("i")), i);
    copy_rvector_to_cppvector(R_do_slot(obj, Rf_install("p")), p);
    copy_rvector_to_cppvector(R_do_slot(obj, Rf_install("x")), x);

    int * dim = INTEGER(R_do_slot(obj, Rf_install("Dim")));
    
    nrow = dim[0];   // sample count,  nrow
    ncol = dim[1];   // feature/gene count,  ncol
    nelem = p[ncol];   // since p is offsets, the ncol+1 entry has the total count

    // get the column names == features.
    // Rcpp::List dimnms = obj.Dimnames;
    // SEXP rownms = VECTOR_ELT(dimnms, 0);    // samples
    // Rcpp::StringVector features = dimnms[1];   // features_names = columnames
    // GET features.
    SEXP dimns = R_do_slot(obj, Rf_install("Dimnames"));  // 
    return VECTOR_ELT(dimns, 1);
}

void copy_rsparsematrix_to_cppvectors(
    SEXP _x, SEXP _i, SEXP _p, 
    std::vector<double> & x,
    std::vector<long> & i,
    std::vector<long> & p) {

    copy_rvector_to_cppvector(_i, i);
    copy_rvector_to_cppvector(_p, p);
    copy_rvector_to_cppvector(_x, x);
}

void copy_rsparsematrix_to_cppvectors(
    SEXP _x, SEXP _i, SEXP _p, 
    std::vector<double> & x,
    std::vector<int> & i,
    std::vector<int> & p) {

    copy_rvector_to_cppvector(_i, i);
    copy_rvector_to_cppvector(_p, p);
    copy_rvector_to_cppvector(_x, x);
}




SEXP copy_rsparsematrix_to_cppvectors(SEXP obj, 
    std::vector<double> & x,
    std::vector<long> & i,
    std::vector<long> & p,
    size_t & nrow, size_t & ncol, size_t & nelem) {

    if (TYPEOF(obj) != S4SXP) return NILSXP;

    copy_rvector_to_cppvector(R_do_slot(obj, Rf_install("i")), i);
    copy_rvector_to_cppvector(R_do_slot(obj, Rf_install("p")), p);
    copy_rvector_to_cppvector(R_do_slot(obj, Rf_install("x")), x);

    int * dim = INTEGER(R_do_slot(obj, Rf_install("Dim")));
    
    nrow = dim[0];   // sample count,  nrow
    ncol = dim[1];   // feature/gene count,  ncol
    nelem = p[ncol];   // since p is offsets, the ncol+1 entry has the total count

    // get the column names == features.
    // Rcpp::List dimnms = obj.Dimnames;
    // SEXP rownms = VECTOR_ELT(dimnms, 0);    // samples
    // Rcpp::StringVector features = dimnms[1];   // features_names = columnames
    // GET features.
    SEXP dimns = R_do_slot(obj, Rf_install("Dimnames"));  // 
    return VECTOR_ELT(dimns, 1);
}




// void import_r_common_params(SEXP as_dataframe, SEXP threads,
//     bool & _as_dataframe, int & nthreads) {
//     _as_dataframe = Rcpp::as<bool>(as_dataframe);
//     nthreads = Rcpp::as<int>(threads);
//     if (nthreads < 1) nthreads = 1;
// }
// void import_de_common_params(SEXP rtype,
//     SEXP bool_param, int & type, bool & bool_val) {

//     type = Rcpp::as<int>(rtype);
//     if (type > 4) Rprintf("ERROR: unsupported type: %d. Supports only less, greater, two.sided, stat, params\n", type);
  
//     bool_val = Rcpp::as<bool>(bool_param);
// }
// void import_fc_common_params(SEXP calc_percents, 
//     SEXP min_threshold, SEXP use_expm1,
//     SEXP use_log, SEXP log_base, SEXP use_pseudocount,
//     bool & perc, double & min_thresh, bool & _use_expm1, 
//     bool & _use_log, double & _log_base, bool & _use_pseudocount) {

//     perc = Rcpp::as<bool>(calc_percents);
//     min_thresh = Rcpp::as<double>(min_threshold);
//     _use_expm1 = Rcpp::as<bool>(use_expm1);
//     _use_log = Rcpp::as<bool>(use_log);
//     _log_base = Rcpp::as<double>(log_base);
//     _use_pseudocount = Rcpp::as<bool>(use_pseudocount);

// }

// void import_filterfc_common_params(SEXP min_pct, 
//     SEXP min_diff_pct, SEXP logfc_threshold,
//     SEXP only_pos, 
//     double & _min_pct, double & _min_diff_pct, 
//     double & _logfc_thresh, bool & _only_pos) {

//     _min_pct = Rcpp::as<double>(min_pct);
//     _min_diff_pct = Rcpp::as<double>(min_diff_pct);
//     _logfc_thresh = Rcpp::as<double>(logfc_threshold);
//     _only_pos = Rcpp::as<bool>(only_pos);
//     }



Rcpp::StringVector populate_feature_names(Rcpp::StringVector const & features,
    size_t const & nfeatures) {
    // GET features.
    // check if features is null.  if so, make a new one.
    // https://stackoverflow.com/questions/25172419/how-can-i-get-the-sexptype-of-an-sexp-value
    Rcpp::StringVector new_features;
    if (features.size() == 0) {
        new_features = Rcpp::StringVector(nfeatures);
        for (size_t j = 0; j < nfeatures; ++j) {
            // create string and set in clust.
            new_features[j] = std::to_string(j+1);
        }
    } else {
        new_features = Rcpp::clone(features);
    }
    return new_features;
}

Rcpp::NumericMatrix export_de_to_r_matrix(
    std::vector<double> const & pv,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    Rcpp::StringVector const & features
) {
    size_t label_count = sorted_labels.size();
    size_t nfeatures = features.size();

    // use clust for column names.
    Rcpp::StringVector clustr(label_count);
    size_t j = 0;
    for (auto item : sorted_labels) {
      clustr[j] = std::to_string(item.first);
      ++j;
    }

    Rcpp::NumericMatrix pvm(label_count, nfeatures, pv.begin());
    Rcpp::rownames(pvm) = clustr;
    Rcpp::colnames(pvm) = features;

    return pvm;
}

Rcpp::DataFrame export_de_to_r_dataframe(
    std::vector<double> const & pv, std::string const & name,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    Rcpp::StringVector const & features
) {
    // https://stackoverflow.com/questions/23547625/returning-a-data-frame-from-c-to-r

    size_t label_count = sorted_labels.size();
    size_t nfeatures = features.size();
    size_t el_count = label_count * nfeatures;

    // data frame has 3 columns:  cluster id, gene names, pvalue.
    
    // --- cluster id, gene names, row names
    Rcpp::IntegerVector clust = Rcpp::IntegerVector(el_count);
    Rcpp::StringVector genenames = Rcpp::StringVector(el_count);
    // Rcpp::StringVector rownames = Rcpp::StringVector(el_count);
    size_t j = 0, f;
    for (f = 0; f < nfeatures; ++f) {
      auto feature = features[f];
      for (auto item : sorted_labels) {
        // rotate through cluster labels for this feature.        
        clust[j] = item.first;
        // same feature name for the set of cluster
        genenames[j] = std::string(feature);
  
        ++j;
      }
    }
    
    // ------ create the dataframe
    Rcpp::NumericVector pvv(pv.size());
    pvv.assign(pv.cbegin(), pv.cend());
    Rcpp::DataFrame res = Rcpp::DataFrame::create(
      Rcpp::Named("cluster") = clust,
      Rcpp::Named("gene") = genenames,
      Rcpp::Named(name) = pvv
    );
    res.attr("row.names") = Rcpp::seq(1, el_count);  // set from 1 to n

    return res;
}

Rcpp::DataFrame export_fc_to_r_dataframe(
    std::vector<double> const & fc, std::string const & fcname,
    std::vector<double> const & p1, std::string const & p1name,
    std::vector<double> const & p2, std::string const & p2name,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    Rcpp::StringVector const & features
) {
    // https://stackoverflow.com/questions/23547625/returning-a-data-frame-from-c-to-r

    size_t label_count = sorted_labels.size();
    size_t nfeatures = features.size();
    size_t el_count = label_count * nfeatures;

    // data frame has 3 columns:  cluster id, gene names, pvalue.
    
    // --- cluster id, gene names, row names
    Rcpp::IntegerVector clust = Rcpp::IntegerVector(el_count);
    Rcpp::StringVector genenames = Rcpp::StringVector(el_count);
    // Rcpp::StringVector rownames = Rcpp::StringVector(el_count);
    size_t j = 0, f;
    for (f = 0; f < nfeatures; ++f) {
      auto feature = features[f];
      for (auto item : sorted_labels) {
        // rotate through cluster labels for this feature.        
        clust[j] = item.first;
        // same feature name for the set of cluster
        genenames[j] = std::string(feature);
  
        ++j;
      }
    }
    
    // ------ create the dataframe
    Rcpp::NumericVector fcv(fc.size());
    fcv.assign(fc.cbegin(), fc.cend());
    Rcpp::NumericVector p1v(p1.size());
    p1v.assign(p1.cbegin(), p1.cend());
    Rcpp::NumericVector p2v(p2.size());
    p2v.assign(p2.cbegin(), p2.cend());
    Rcpp::DataFrame res = Rcpp::DataFrame::create(
      Rcpp::Named("cluster") = clust,
      Rcpp::Named("gene") = genenames,
      Rcpp::Named(fcname) = fcv,
      Rcpp::Named(p1name) = p1v,
      Rcpp::Named(p2name) = p2v
    );
    res.attr("row.names") = Rcpp::seq(1, el_count);  // set from 1 to n

    return res;
}

Rcpp::List export_fc_to_r_matrix(
    std::vector<double> const & fc, std::string const & fcname,
    std::vector<double> const & p1, std::string const & p1name,
    std::vector<double> const & p2, std::string const & p2name,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    Rcpp::StringVector const & features
) {
    size_t label_count = sorted_labels.size();
    size_t nfeatures = features.size();

    // use clust for column names.
    Rcpp::StringVector clustr(label_count);
    size_t j = 0;
    for (auto item : sorted_labels) {
      clustr[j] = std::to_string(item.first);
      ++j;
    }

    Rcpp::NumericMatrix fcm(label_count, nfeatures, fc.begin());
    Rcpp::rownames(fcm) = clustr;
    Rcpp::colnames(fcm) = features;

    Rcpp::NumericMatrix p1m(label_count, nfeatures, p1.begin());
    Rcpp::rownames(p1m) = clustr;
    Rcpp::colnames(p1m) = features;
    
    Rcpp::NumericMatrix p2m(label_count, nfeatures, p2.begin());
    Rcpp::rownames(p2m) = clustr;
    Rcpp::colnames(p2m) = features;

    Rcpp::List res = Rcpp::List::create(
      Rcpp::Named(fcname) = fcm,
      Rcpp::Named(p1name) = p1m,
      Rcpp::Named(p2name) = p2m
    );

    return res;
}

Rcpp::List export_fc_to_r_matrix(
    std::vector<double> const & fc, std::string const & fcname,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    Rcpp::StringVector const & features
) {
    size_t label_count = sorted_labels.size();
    size_t nfeatures = features.size();

    // use clust for column names.
    Rcpp::StringVector clustr(label_count);
    size_t j = 0;
    for (auto item : sorted_labels) {
      clustr[j] = std::to_string(item.first);
      ++j;
    }

    Rcpp::NumericMatrix fcm(label_count, nfeatures, fc.begin());
    Rcpp::rownames(fcm) = clustr;
    Rcpp::colnames(fcm) = features;

    Rcpp::List res = Rcpp::List::create(
      Rcpp::Named(fcname) = fcm
    );

    return res;
}



SEXP export_de_to_r_dataframe(
    std::vector<double> const & pv, std::string const & name,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    SEXP features
) {
    // https://stackoverflow.com/questions/23547625/returning-a-data-frame-from-c-to-r
    // https://coolbutuseless.github.io/2020/09/16/creating-a-data.frame-in-c/


    int intlen = snprintf(NULL, 0, "%lu", std::numeric_limits<size_t>::max());
    char * str = reinterpret_cast<char *>(malloc(intlen + 1));
    int proc_count = 0;
    
    size_t label_count = sorted_labels.size();
    size_t el_count = pv.size();
    size_t nfeatures = el_count / label_count;

    if (TYPEOF(features) == NILSXP) {
        PROTECT(features = Rf_allocVector(STRSXP, nfeatures));
        ++proc_count;
        
        for (size_t j = 0; j < nfeatures; ++j) {
        // create string and set in clust.
        sprintf(str, "%lu", j);
        SET_STRING_ELT(features, j, Rf_mkChar(str));
        memset(str, 0, intlen + 1);
        }
    }

    SEXP fc, clust, genenames, res, names, rownames, cls, dimnames;
    int ncols = 3;
    int col_id = 0;

    PROTECT(res = Rf_allocVector(VECSXP, ncols));   // cluster id and gene names are repeated.
    PROTECT(names = Rf_allocVector(STRSXP, ncols));  // col names
    proc_count += 2;

    PROTECT(cls = Rf_allocVector(STRSXP, 1)); // class attribute
    SET_STRING_ELT(cls, 0, Rf_mkChar("data.frame"));
    Rf_classgets(res, cls);
    ++proc_count;
    // Rf_classgets(res, Rf_install("data.frame"));

    PROTECT(clust = Rf_allocVector(INTSXP, el_count));
    PROTECT(genenames = Rf_allocVector(STRSXP, el_count));
    PROTECT(fc = Rf_allocVector(REALSXP, el_count));
    proc_count += 3;

    double * fc_ptr = REAL(fc);
    std::copy(pv.begin(), pv.end(), fc_ptr);

    PROTECT(rownames = Rf_allocVector(STRSXP, label_count * nfeatures));  // dataframe column names.
    ++proc_count;

    // make the clusters vector.
    int * clust_ptr = INTEGER(clust);
    SEXP * features_ptr = STRING_PTR(features);
    size_t j = 0;
    // outer group = features, inner order = cluster
    for (size_t i = 0; i < nfeatures; ++i) {
      for (auto item : sorted_labels) {
        // rotate through cluster labels for this feature.        
        *clust_ptr = item.first;
        ++clust_ptr;

        // same feature name for the set of cluster
        SET_STRING_ELT(genenames, j, Rf_duplicate(*features_ptr));

        sprintf(str, "%lu", j);
        SET_STRING_ELT(rownames, j, Rf_mkChar(str));
        memset(str, 0, intlen + 1);
  
        ++j;
      }
      ++features_ptr;
    }

    SET_VECTOR_ELT(res, col_id, clust);  
    SET_STRING_ELT(names, col_id, Rf_mkChar("cluster"));
    ++col_id;
    SET_VECTOR_ELT(res, col_id, genenames);
    SET_STRING_ELT(names, col_id, Rf_mkChar("gene"));
    ++col_id;

    SET_VECTOR_ELT(res, col_id, fc);
    // convert from STRSXP to CHARSXP
    SET_STRING_ELT(names, col_id, Rf_mkChar(name.c_str())); 
    ++col_id;

    Rf_setAttrib(res, R_NamesSymbol, names);

    // set row names - NEEDED to print the dataframe!!!
    Rf_setAttrib(res, R_RowNamesSymbol, rownames);

    // //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // // Set the row.names on the list.
    // // Use the shortcut as used in .set_row_names() in R
    // // i.e. set rownames to c(NA_integer, -len) and it will
    // // take care of the rest. This is equivalent to rownames(x) <- NULL
    // //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // PROTECT(rownames = Rf_allocVector(INTSXP, 2));
    // SET_INTEGER_ELT(rownames, 0, NA_INTEGER);
    // SET_INTEGER_ELT(rownames, 1, -(static_cast<long>(el_count));
    // Rf_setAttrib(res, R_RowNamesSymbol, rownames);
    // ++proc_count;

    UNPROTECT(proc_count);
    free(str);
    // end = Clock::now();
    // Rprintf("cleaned up %ld ns\n", 
    // std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    return(res);

}

SEXP export_fc_to_r_dataframe(
    std::vector<double> const & fc, std::string const & fcname,
    std::vector<double> const & p1, std::string const & p1name,
    std::vector<double> const & p2, std::string const & p2name,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    SEXP features
) {
    // https://stackoverflow.com/questions/23547625/returning-a-data-frame-from-c-to-r
    // https://coolbutuseless.github.io/2020/09/16/creating-a-data.frame-in-c/


    int intlen = snprintf(NULL, 0, "%lu", std::numeric_limits<size_t>::max());
    char * str = reinterpret_cast<char *>(malloc(intlen + 1));
    int proc_count = 0;
    
    size_t label_count = sorted_labels.size();
    size_t el_count = fc.size();
    size_t nfeatures = el_count / label_count;

    if (TYPEOF(features) == NILSXP) {
        PROTECT(features = Rf_allocVector(STRSXP, nfeatures));
        ++proc_count;
        
        for (size_t j = 0; j < nfeatures; ++j) {
        // create string and set in clust.
        sprintf(str, "%lu", j);
        SET_STRING_ELT(features, j, Rf_mkChar(str));
        memset(str, 0, intlen + 1);
        }
    }

    SEXP _fc, _p1, _p2, clust, genenames, res, names, rownames, cls, dimnames;
    int ncols = 5;
    int col_id = 0;

    PROTECT(res = Rf_allocVector(VECSXP, ncols));   // cluster id and gene names are repeated.
    PROTECT(names = Rf_allocVector(STRSXP, ncols));  // col names
    proc_count += 2;

    PROTECT(cls = Rf_allocVector(STRSXP, 1)); // class attribute
    SET_STRING_ELT(cls, 0, Rf_mkChar("data.frame"));
    Rf_classgets(res, cls);
    ++proc_count;
    // Rf_classgets(res, Rf_install("data.frame"));

    PROTECT(clust = Rf_allocVector(INTSXP, el_count));
    PROTECT(genenames = Rf_allocVector(STRSXP, el_count));
    PROTECT(_fc = Rf_allocVector(REALSXP, el_count));
    proc_count += 3;

    PROTECT(_p1 = Rf_allocVector(REALSXP, el_count));
    PROTECT(_p2 = Rf_allocVector(REALSXP, el_count)); 
    proc_count += 2;


    double * fc_ptr = REAL(_fc);
    std::copy(fc.begin(), fc.end(), fc_ptr);

    double * p1_ptr = REAL(_p1);
    std::copy(p1.begin(), p1.end(), p1_ptr);
    double * p2_ptr = REAL(_p2);
    std::copy(p2.begin(), p2.end(), p2_ptr);

    PROTECT(rownames = Rf_allocVector(STRSXP, label_count * nfeatures));  // dataframe column names.
    ++proc_count;

    // make the clusters vector.
    int * clust_ptr = INTEGER(clust);
    SEXP * features_ptr = STRING_PTR(features);
    size_t j = 0;
    // outer group = features, inner order = cluster
    for (size_t i = 0; i < nfeatures; ++i) {
      for (auto item : sorted_labels) {
        // rotate through cluster labels for this feature.        
        *clust_ptr = item.first;
        ++clust_ptr;

        // same feature name for the set of cluster
        SET_STRING_ELT(genenames, j, Rf_duplicate(*features_ptr));

        sprintf(str, "%lu", j);
        SET_STRING_ELT(rownames, j, Rf_mkChar(str));
        memset(str, 0, intlen + 1);
  
        ++j;
      }
      ++features_ptr;
    }

    SET_VECTOR_ELT(res, col_id, clust);  
    SET_STRING_ELT(names, col_id, Rf_mkChar("cluster"));
    ++col_id;
    SET_VECTOR_ELT(res, col_id, genenames);
    SET_STRING_ELT(names, col_id, Rf_mkChar("gene"));
    ++col_id;

    SET_VECTOR_ELT(res, col_id, _fc);
    // convert from STRSXP to CHARSXP
    SET_STRING_ELT(names, col_id, Rf_mkChar(fcname.c_str())); 
    ++col_id;

    SET_VECTOR_ELT(res, col_id, _p1);
    SET_STRING_ELT(names, col_id, Rf_mkChar(p1name.c_str())); 
    ++col_id;
    SET_VECTOR_ELT(res, col_id, _p2);
    SET_STRING_ELT(names, col_id, Rf_mkChar(p2name.c_str())); 
    ++col_id;


    Rf_setAttrib(res, R_NamesSymbol, names);

    // set row names - NEEDED to print the dataframe!!!
    Rf_setAttrib(res, R_RowNamesSymbol, rownames);

    // //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // // Set the row.names on the list.
    // // Use the shortcut as used in .set_row_names() in R
    // // i.e. set rownames to c(NA_integer, -len) and it will
    // // take care of the rest. This is equivalent to rownames(x) <- NULL
    // //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // PROTECT(rownames = Rf_allocVector(INTSXP, 2));
    // SET_INTEGER_ELT(rownames, 0, NA_INTEGER);
    // SET_INTEGER_ELT(rownames, 1, -n);
    // Rf_setAttrib(res, R_RowNamesSymbol, rownames);
    // ++proc_count;

    UNPROTECT(proc_count);
    free(str);
    // end = Clock::now();
    // Rprintf("cleaned up %ld ns\n", 
    // std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    return(res);

}


SEXP export_fc_to_r_matrix(
    std::vector<double> const & fc, std::string const & fcname,
    std::vector<double> const & p1, std::string const & p1name,
    std::vector<double> const & p2, std::string const & p2name,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    SEXP features
) {
    
    int intlen = snprintf(NULL, 0, "%lu", std::numeric_limits<size_t>::max());
    char * str = reinterpret_cast<char *>(malloc(intlen + 1));
    int proc_count = 0;
    
    size_t label_count = sorted_labels.size();
    size_t nfeatures = Rf_xlength(features);
    size_t nelem = fc.size();

    if (TYPEOF(features) == NILSXP) {
        PROTECT(features = Rf_allocVector(STRSXP, nfeatures));
        ++proc_count;
        
        for (size_t j = 0; j < nfeatures; ++j) {
        // create string and set in clust.
        sprintf(str, "%lu", j);
        SET_STRING_ELT(features, j, Rf_mkChar(str));
        memset(str, 0, intlen + 1);
        }
    }

    SEXP _fc, _p1, _p2, clust, genenames, res, names, cls, dimnames;
    int ncols = 3;
    int col_id = 0;

    PROTECT(res = Rf_allocVector(VECSXP, ncols));   // cluster id and gene names are repeated.
    PROTECT(names = Rf_allocVector(STRSXP, ncols));  // col names
    proc_count += 2;

    // use clust for column names.
    PROTECT(clust = Rf_allocVector(STRSXP, label_count));
    ++proc_count;
    size_t j = 0;
    for (auto item : sorted_labels) {
      sprintf(str, "%d", item.first);
      SET_STRING_ELT(clust, j, Rf_mkChar(str));
      memset(str, 0, intlen + 1);
      ++j;
    }
    
    PROTECT(_fc = Rf_allocMatrix(REALSXP, label_count, nfeatures));
    ++proc_count;
    PROTECT(_p1 = Rf_allocMatrix(REALSXP, label_count, nfeatures));
    PROTECT(_p2 = Rf_allocMatrix(REALSXP, label_count, nfeatures)); 
    proc_count += 2;

    double * fc_ptr = REAL(_fc);
    std::copy(fc.begin(), fc.end(), fc_ptr);
    double * p1_ptr = REAL(_p1);
    std::copy(p1.begin(), p1.end(), p1_ptr);
    double * p2_ptr = REAL(_p2);
    std::copy(p2.begin(), p2.end(), p2_ptr);

    SET_VECTOR_ELT(res, col_id, _fc);
    // convert from STRSXP to CHARSXP
    SET_STRING_ELT(names, col_id, Rf_mkChar(fcname.c_str())); 
    ++col_id;

    SET_VECTOR_ELT(res, col_id, _p1);
    SET_STRING_ELT(names, col_id, Rf_mkChar(p1name.c_str()));
    ++col_id;
    SET_VECTOR_ELT(res, col_id, _p2);
    SET_STRING_ELT(names, col_id, Rf_mkChar(p2name.c_str()));
    ++col_id;

    // set list element names.
    Rf_setAttrib(res, R_NamesSymbol, names);

    PROTECT(dimnames = Rf_allocVector(VECSXP, 2));
    ++proc_count;
    SET_VECTOR_ELT(dimnames, 0, clust);  // rows = clusters
    SET_VECTOR_ELT(dimnames, 1, features);  // columns  = features (genes)

    // https://stackoverflow.com/questions/5709940/r-extension-in-c-setting-matrix-row-column-names
    Rf_setAttrib(_fc, R_DimNamesSymbol, dimnames);
    Rf_setAttrib(_p1, R_DimNamesSymbol, Rf_duplicate(dimnames));
    Rf_setAttrib(_p2, R_DimNamesSymbol, Rf_duplicate(dimnames));

    UNPROTECT(proc_count);
    free(str);
    // end = Clock::now();
    // Rprintf("cleaned up %ld ns\n", 
    // std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    return(res);


}

SEXP export_fc_to_r_matrix(
    std::vector<double> const & fc, std::string const & fcname,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    SEXP features
) {


    int intlen = snprintf(NULL, 0, "%lu", std::numeric_limits<size_t>::max());
    char * str = reinterpret_cast<char *>(malloc(intlen + 1));
    int proc_count = 0;
    
    size_t label_count = sorted_labels.size();
    size_t nfeatures = Rf_xlength(features);
    size_t nelem = fc.size();

    if (TYPEOF(features) == NILSXP) {
        PROTECT(features = Rf_allocVector(STRSXP, nfeatures));
        ++proc_count;
        
        for (size_t j = 0; j < nfeatures; ++j) {
        // create string and set in clust.
        sprintf(str, "%lu", j);
        SET_STRING_ELT(features, j, Rf_mkChar(str));
        memset(str, 0, intlen + 1);
        }
    }

    SEXP _fc, _p1, _p2, clust, genenames, res, names, cls, dimnames;
    int ncols = 1;
    int col_id = 0;

    PROTECT(res = Rf_allocVector(VECSXP, ncols));   // cluster id and gene names are repeated.
    PROTECT(names = Rf_allocVector(STRSXP, ncols));  // col names
    proc_count += 2;

    // use clust for column names.
    PROTECT(clust = Rf_allocVector(STRSXP, label_count));
    ++proc_count;
    size_t j = 0;
    for (auto item : sorted_labels) {
      sprintf(str, "%d", item.first);
      SET_STRING_ELT(clust, j, Rf_mkChar(str));
      memset(str, 0, intlen + 1);
      ++j;
    }
    
    PROTECT(_fc = Rf_allocMatrix(REALSXP, label_count, nfeatures));
    ++proc_count;

    double * fc_ptr = REAL(_fc);
    std::copy(fc.begin(), fc.end(), fc_ptr);

    SET_VECTOR_ELT(res, col_id, _fc);
    // convert from STRSXP to CHARSXP
    SET_STRING_ELT(names, col_id, Rf_mkChar(fcname.c_str())); 
    ++col_id;

    // set list element names.
    Rf_setAttrib(res, R_NamesSymbol, names);

    PROTECT(dimnames = Rf_allocVector(VECSXP, 2));
    ++proc_count;
    SET_VECTOR_ELT(dimnames, 0, clust);  // rows = clusters
    SET_VECTOR_ELT(dimnames, 1, features);  // columns  = features (genes)

    // https://stackoverflow.com/questions/5709940/r-extension-in-c-setting-matrix-row-column-names
    Rf_setAttrib(_fc, R_DimNamesSymbol, dimnames);

    UNPROTECT(proc_count);
    free(str);
    // end = Clock::now();
    // Rprintf("cleaned up %ld ns\n", 
    // std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    return(res);


}
