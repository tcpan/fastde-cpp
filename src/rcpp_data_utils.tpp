#pragma once

// ------- function definition

#include "rcpp_data_utils.hpp"

// ------- class def
namespace Rcpp {
    dgCMatrix::dgCMatrix(Rcpp::S4 const & mat) :
        i(mat.slot("i")),   // row id
        p(mat.slot("p")),   // offsets for starts of columns
        x(mat.slot("x")),
        Dim(mat.slot("Dim")),
        Dimnames(mat.slot("Dimnames")) {};
    dgCMatrix::dgCMatrix(int const & nrow, int const & ncol, int const & nelem) :
        i(nelem),  // row id
        p(ncol + 1),  // offsets for starts of columns
        x(nelem),
        Dim(2) {
        Dim[0] = nrow;
        Dim[1] = ncol;
    };


    // specialization of Rcpp::as
    template <> dgCMatrix as(SEXP mat) { return dgCMatrix(mat); };

    // specialization of Rcpp::wrap
    template <> SEXP wrap(const dgCMatrix& sm) {
        S4 s(std::string("dgCMatrix"));
        s.slot("i") = sm.i;
        s.slot("p") = sm.p;
        s.slot("x") = sm.x;
        s.slot("Dim") = sm.Dim;
        s.slot("Dimnames") = sm.Dimnames;
        return s;
    };

    dgCMatrix64::dgCMatrix64(Rcpp::S4 const & mat) :
        i(mat.slot("i")),
        p(mat.slot("p")),
        x(mat.slot("x")),
        dimension(mat.slot("dimension")),
        Dim(mat.slot("Dim")),
        Dimnames(mat.slot("Dimnames"))
    {};
    dgCMatrix64::dgCMatrix64(long const & nrow, long const & ncol, long const & nelem) :
        i(nelem),  // row id
        p(ncol + 1),  // offsets for starts of columns
        x(nelem),
        dimension(2),
        Dim(2) {
        dimension[0] = nrow;
        dimension[1] = ncol;
        Dim[0] = (nrow > 2147483647) ? -1 : nrow;
        Dim[1] = (ncol > 2147483647) ? -1 : ncol;
    };

    // specialization of Rcpp::as 
    template <> dgCMatrix64 as(SEXP mat) { return dgCMatrix64(mat); };

    template <> SEXP wrap(const dgCMatrix64& sm) {
        S4 s(std::string("dgCMatrix64"));
        s.slot("i") = sm.i;
        s.slot("p") = sm.p;
        s.slot("x") = sm.x;
        s.slot("dimension") = sm.dimension;
        s.slot("Dim") = sm.Dim;
        s.slot("Dimnames") = sm.Dimnames;
        return s;
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


Rcpp::StringVector copy_rsparsematrix_to_cppvectors(Rcpp::dgCMatrix const & obj, 
    std::vector<double> & x,
    std::vector<int> & i,
    std::vector<int> & p,
    size_t & nrow, size_t & ncol, size_t & nelem) {

    copy_rvector_to_cppvector(obj.i, i);
    copy_rvector_to_cppvector(obj.p, p);
    copy_rvector_to_cppvector(obj.x, x);
    
    Rcpp::IntegerVector dim = obj.Dim;
    nrow = dim[0];   // sample count,  nrow
    ncol = dim[1];   // feature/gene count,  ncol
    nelem = p[ncol];   // since p is offsets, the ncol+1 entry has the total count

    // get the column names == features.
    Rcpp::List dimnms = obj.Dimnames;
    // SEXP rownms = VECTOR_ELT(dimnms, 0);    // samples
    // Rcpp::StringVector features = dimnms[1];   // features_names = columnames
    // GET features.
    return dimnms[1];  // 
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
        obj.x, obj.i, obj.p, 
        x, i, p
    );

    Rcpp::NumericVector dim = obj.dimension;
    nrow = dim[0];   // sample count,  nrow
    ncol = dim[1];   // feature/gene count,  ncol
    nelem = p[ncol];   // since p is offsets, the ncol+1 entry has the total count

    Rcpp::List dimnms = obj.Dimnames;
    // SEXP rownms = VECTOR_ELT(dimnms, 0);    // samples
    // Rcpp::StringVector features = dimnms[1];   // features_names = columnames
    // GET features.
    return dimnms[1];  // 
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
