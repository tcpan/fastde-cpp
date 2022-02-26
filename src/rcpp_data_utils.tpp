#pragma once

// ------- function definition

#include "rcpp_data_utils.hpp"


Rcpp::StringVector rmatrix_to_vector(SEXP matrix, std::vector<double> & mat,
    size_t & nrow, size_t & ncol, size_t & nelem) {
    Rcpp::NumericMatrix _matrix(matrix);
    nrow=_matrix.nrow(); // n
    ncol=_matrix.ncol(); // m
    nelem = nrow * ncol;

    mat.clear();
    mat.insert(mat.end(), _matrix.cbegin(), _matrix.cend());

    return Rcpp::colnames(_matrix);
}

void rvector_to_vector(SEXP vect, std::vector<double> & vec) {
    Rcpp::NumericVector _vector(vect);
    
    vec.clear();
    vec.insert(vec.end(), _vector.cbegin(), _vector.cend());
}

void rvector_to_vector(SEXP vect, std::vector<int> & vec) {
    Rcpp::IntegerVector _vector(vect);
    
    vec.clear();
    vec.insert(vec.end(), _vector.cbegin(), _vector.cend());
}

size_t rvector_to_vector(SEXP vect, std::vector<double> & vec, size_t const & length, size_t const & offset) {
    Rcpp::NumericVector _vector(vect);
    size_t end_offset = std::min(offset + length, static_cast<size_t>(_vector.size()));
    vec.clear();
    vec.insert(vec.end(), _vector.cbegin() + offset, _vector.cbegin() + end_offset);

    return end_offset - offset;
}

size_t rvector_to_vector(SEXP vect, std::vector<int> & vec, size_t const & length, size_t const & offset) {
    Rcpp::IntegerVector _vector(vect);
    size_t end_offset = std::min(offset + length, static_cast<size_t>(_vector.size()));
    vec.clear();
    vec.insert(vec.end(), _vector.cbegin() + offset, _vector.cbegin() + end_offset);

    return end_offset - offset;
}


Rcpp::StringVector rsparsematrix_to_vectors(SEXP matrix, 
    std::vector<double> & x,
    std::vector<int> & i,
    std::vector<int> & p,
    size_t & nrow, size_t & ncol, size_t & nelem) {

    Rcpp::S4 obj(matrix);
    
    rvector_to_vector(obj.slot("i"), i);
    rvector_to_vector(obj.slot("p"), p);
    rvector_to_vector(obj.slot("x"), x);
    
    Rcpp::IntegerVector dim = obj.slot("Dim");
    nrow = dim[0];   // sample count,  nrow
    ncol = dim[1];   // feature/gene count,  ncol
    nelem = p[ncol];   // since p is offsets, the ncol+1 entry has the total count

    // get the column names == features.
    Rcpp::List dimnms = obj.slot("Dimnames");
    // SEXP rownms = VECTOR_ELT(dimnms, 0);    // samples
    // Rcpp::StringVector features = dimnms[1];   // features_names = columnames
    // GET features.
    return dimnms[1];  // 
}

void import_r_common_params(SEXP alternative,
    SEXP var_equal, SEXP as_dataframe, SEXP threads,
    int & type, bool & var_eq, bool & _as_dataframe, int & nthreads) {

    type = Rcpp::as<int>(alternative);
    if (type > 3) Rprintf("ERROR: unsupported type: %d. Supports only less, greater, two.sided, tstat\n", type);
  
    var_eq = Rcpp::as<bool>(var_equal);

    _as_dataframe = Rcpp::as<bool>(as_dataframe);

    nthreads = Rcpp::as<int>(threads);
    if (nthreads < 1) nthreads = 1;
}


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
    std::vector<std::pair<int, size_t>> const & sorted_labels,
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
    std::vector<double> const & pv,
    std::vector<std::pair<int, size_t>> const & sorted_labels,
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
      Rcpp::Named("p_val") = pvv
    );
    res.attr("row.names") = Rcpp::seq(1, el_count);  // set from 1 to n

    return res;
}