#pragma once

// ------- function definition

#include "utils_data.hpp"



// ------- function def


void copy_rmatrix_to_cppvector(
    cpp11::doubles_matrix<cpp11::by_column> const & _matrix, 
    std::vector<double> & mat) {
    size_t nrow=_matrix.nrow(); // n
    size_t ncol=_matrix.ncol(); // m

    mat.clear();
    mat.reserve(nrow * ncol);
    auto last = _matrix.end();
    for (auto slice = _matrix.begin(); slice != last; ++slice) {
        mat.insert(mat.end(), (*slice).begin(), (*slice).end());
    }
}
void copy_rmatrix_to_cppvector(
    cpp11::logicals_matrix<cpp11::by_column> const & _matrix, 
    std::vector<bool> & mat) {
    size_t nrow=_matrix.nrow(); // n
    size_t ncol=_matrix.ncol(); // m

    mat.clear();
    mat.reserve(nrow * ncol);
    auto last = _matrix.end();
    for (auto slice = _matrix.begin(); slice != last; ++slice) {
        mat.insert(mat.end(), (*slice).begin(), (*slice).end());
    }
}


void copy_rmatrix_to_cppvector(
    cpp11::doubles_matrix<cpp11::by_column> const & _matrix, 
    double * mat) {
    auto last = _matrix.end();
    for (auto slice = _matrix.begin(); slice != last; ++slice) {
        std::copy((*slice).begin(), (*slice).end(), mat);
        mat += (*slice).size();
    }
}
void copy_rmatrix_to_cppvector(
    cpp11::logicals_matrix<cpp11::by_column> const & _matrix, 
    bool * mat) {
    auto last = _matrix.end();
    for (auto slice = _matrix.begin(); slice != last; ++slice) {
        std::copy((*slice).begin(), (*slice).end(), mat);
        mat += (*slice).size();
    }
}

size_t copy_rvector_to_cppvector(cpp11::logicals const & _vector, 
    std::vector<bool> & vec, size_t const & length, size_t const & offset) {
    size_t veclen = static_cast<size_t>(_vector.size());
    if (offset >= veclen) return 0;  // offset is greater than vector length
    size_t len = std::min(length, veclen - offset);   // length to return
    vec.clear();
    vec.insert(vec.end(), _vector.begin() + offset, _vector.begin() + offset + len);

    return len;
}

size_t copy_rvector_to_cppvector(cpp11::doubles const & _vector, 
    std::vector<double> & vec, size_t const & length, size_t const & offset) {
    size_t veclen = static_cast<size_t>(_vector.size());
    if (offset >= veclen) return 0;  // offset is greater than vector length
    size_t len = std::min(length, veclen - offset);   // length to return
    vec.clear();
    vec.insert(vec.end(), _vector.begin() + offset, _vector.begin() + offset + len);

    return len;
}

size_t copy_rvector_to_cppvector(cpp11::integers const & _vector, 
    std::vector<int> & vec, size_t const & length, size_t const & offset) {
    size_t veclen = static_cast<size_t>(_vector.size());
    if (offset >= veclen) return 0;  // offset is greater than vector length
    size_t len = std::min(length, veclen - offset);   // length to return
    vec.clear();
    vec.insert(vec.end(), _vector.begin() + offset, _vector.begin() + offset + len);

    return len;
}
size_t copy_rvector_to_cppvector(cpp11::doubles const & _vector, 
    std::vector<long> & vec, size_t const & length, size_t const & offset) {
    size_t veclen = static_cast<size_t>(_vector.size());
    if (offset >= veclen) return 0;  // offset is greater than vector length
    size_t len = std::min(length, veclen - offset);   // length to return
    vec.clear();
    vec.insert(vec.end(), _vector.begin() + offset, _vector.begin() + offset + len);

    return len;
}


size_t copy_rvector_to_cppvector(cpp11::logicals const & _vector, 
    bool * vec, size_t const & length, size_t const & offset) {
    size_t veclen = static_cast<size_t>(_vector.size());
    if (offset >= veclen) return 0;  // offset is greater than vector length
    size_t len = std::min(length, veclen - offset);   // length to return
    std::copy(_vector.begin() + offset, _vector.begin() + offset + len, vec);

    return len;
}

size_t copy_rvector_to_cppvector(cpp11::doubles const & _vector, 
    double * vec, size_t const & length, size_t const & offset) {
    size_t veclen = static_cast<size_t>(_vector.size());
    if (offset >= veclen) return 0;  // offset is greater than vector length
    size_t len = std::min(length, veclen - offset);   // length to return
    std::copy(_vector.begin() + offset, _vector.begin() + offset + len, vec);

    return len;
}

size_t copy_rvector_to_cppvector(cpp11::integers const & _vector, 
    int * vec, size_t const & length, size_t const & offset) {
    size_t veclen = static_cast<size_t>(_vector.size());
    if (offset >= veclen) return 0;  // offset is greater than vector length
    size_t len = std::min(length, veclen - offset);   // length to return
    std::copy(_vector.begin() + offset, _vector.begin() + offset + len, vec);

    return len;
}
size_t copy_rvector_to_cppvector(cpp11::doubles const & _vector, 
    long * vec, size_t const & length, size_t const & offset) {
    size_t veclen = static_cast<size_t>(_vector.size());
    if (offset >= veclen) return 0;  // offset is greater than vector length
    size_t len = std::min(length, veclen - offset);   // length to return
    std::copy(_vector.begin() + offset, _vector.begin() + offset + len, vec);

    return len;
}

template <typename R_VEC, typename ITER>
R_VEC export_vec_to_rvec(ITER pv, size_t const & count) {
    return R_VEC(pv, pv + count);
}

// template <typename R_VEC, typename T>
// R_VEC export_vec_to_rvec(std::vector<T> const & pv) {
//     return export_vec_to_rvec<R_VEC>(pv.data(), pv.size());
// }




// no col or row names.
// pv is organized as a linearized matrix, with labels consecutive, and and strides of label_count.  
// output should be column major (i.e. rows are consecutive.)
// so the matrix would have labels in rows, and features in columns. 
template <typename R_MAT, typename ITER>
R_MAT export_vec_to_r_matrix(
    ITER pv, size_t const & nrow, size_t const & ncol
) {
    R_MAT pvm(nrow, ncol);
    auto slice_end = pvm.end();
    auto vec = pv;
    auto vec_end = pv + nrow;
    for (auto slice = pvm.begin(); slice != slice_end; ++slice) {
        std::copy(vec, vec_end, (*slice).begin());
        vec = vec_end;
        vec_end += nrow;
    }

    return pvm;
}


// no col or row names.
// pv is organized as a linearized matrix, with labels consecutive, and and strides of label_count.  
// output should be column major (i.e. rows are consecutive.)
// so the matrix would have labels in rows, and features in columns. 
template <typename R_MAT, typename T>
R_MAT export_vec_to_r_matrix(
    std::vector<T> const & pv, size_t const & nrow, size_t const & ncol
) {
    return export_vec_to_r_matrix<R_MAT>(pv.begin(), nrow, ncol);
}


template <typename ITER>
cpp11::writable::list export_fc_to_r_matrix(
    ITER fc, std::string const & fcname,
    ITER p1, std::string const & p1name,
    ITER p2, std::string const & p2name,
    size_t const & count,
    std::vector<std::pair<int, size_t> > const & sorted_labels
) {
    cpp11::writable::doubles_matrix<cpp11::by_column> fcm = export_vec_to_r_matrix<cpp11::writable::doubles_matrix<cpp11::by_column>>(fc, sorted_labels.size(), count/sorted_labels.size());
    cpp11::writable::doubles_matrix<cpp11::by_column> p1m = export_vec_to_r_matrix<cpp11::writable::doubles_matrix<cpp11::by_column>>(p1, sorted_labels.size(), count/sorted_labels.size());
    cpp11::writable::doubles_matrix<cpp11::by_column> p2m = export_vec_to_r_matrix<cpp11::writable::doubles_matrix<cpp11::by_column>>(p2, sorted_labels.size(), count/sorted_labels.size());

    cpp11::named_arg _fc(fcname.c_str()); _fc = fcm;
    cpp11::named_arg _p1(p1name.c_str()); _p1 = p1m;
    cpp11::named_arg _p2(p2name.c_str()); _p2 = p2m;
    return cpp11::writable::list( { _fc, _p1, _p2 } );
}

template <typename T>
cpp11::writable::list export_fc_to_r_matrix(
    std::vector<T> const & fc, std::string const & fcname,
    std::vector<T> const & p1, std::string const & p1name,
    std::vector<T> const & p2, std::string const & p2name,
    std::vector<std::pair<int, size_t> > const & sorted_labels
) {
    return export_fc_to_r_matrix(
        fc.data(), fcname,
        p1.data(), p1name,
        p2.data(), p2name,
        fc.size(), sorted_labels
    );

}



// pv is organized as a linearized matrix, with labels consecutive, and strides of label_count.   
// so data frame iteration needs to be nfeatures in outer loop, and label_count in inner.
template <typename ITER>
cpp11::writable::data_frame export_vec_to_r_dataframe(
    ITER pv, std::string const & name, size_t const & el_count,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    cpp11::strings const & features
) {
    // https://stackoverflow.com/questions/23547625/returning-a-data-frame-from-c-to-r

    size_t label_count = sorted_labels.size();
    size_t nfeatures = el_count / label_count;
    
    // data frame has 3 columns:  cluster id, gene names, pvalue.
    
    // --- cluster id, gene names, row names
    cpp11::writable::integers clust(el_count);
    cpp11::writable::strings genenames(el_count);

    auto clust_i = clust.begin();
    auto genenames_i = genenames.begin();
    
    auto features_end = features.end();
    for (auto features_i = features.begin(); features_i != features_end; ++features_i) {
      for (auto item : sorted_labels) {
        // rotate through cluster labels for this feature.        
        *clust_i = item.first;
        // same feature name for the set of cluster
        *genenames_i = *features_i;
  
        ++clust_i;
        ++genenames_i;
      }
    }
    
    cpp11::named_arg _fc(name.c_str()); _fc = export_vec_to_rvec<cpp11::writable::doubles>(pv, el_count);
    cpp11::named_arg _cl("cluster"); _cl = clust;
    cpp11::named_arg _gn("gene"); _gn = genenames;
    return cpp11::writable::data_frame( {_cl, _gn, _fc} );

    // ------ create the dataframe
}


// pv is organized as a linearized matrix, with labels consecutive, and strides of label_count.   
// so data frame iteration needs to be nfeatures in outer loop, and label_count in inner.
template <typename T>
cpp11::writable::data_frame export_vec_to_r_dataframe(
    std::vector<T> const & pv, std::string const & name,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    cpp11::strings const & features
) {
    return export_vec_to_r_dataframe(pv.data(), name, pv.size(),
        sorted_labels, features);
}

template <typename ITER>
cpp11::writable::data_frame export_fc_to_r_dataframe(
    ITER fc, std::string const & fcname,
    ITER p1, std::string const & p1name,
    ITER p2, std::string const & p2name,
    size_t const & el_count,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    cpp11::strings const & features
) {
    // https://stackoverflow.com/questions/23547625/returning-a-data-frame-from-c-to-r
    size_t label_count = sorted_labels.size();
    size_t nfeatures = el_count / label_count;
    
    // data frame has 3 columns:  cluster id, gene names, pvalue.
    
    // --- cluster id, gene names, row names
    cpp11::writable::integers clust(el_count);
    cpp11::writable::strings genenames(el_count);

    auto clust_i = clust.begin();
    auto genenames_i = genenames.begin();
    
    auto features_end = features.end();
    for (auto features_i = features.begin(); features_i != features_end; ++features_i) {
      for (auto item : sorted_labels) {
        // rotate through cluster labels for this feature.        
        *clust_i = item.first;
        // same feature name for the set of cluster
        *genenames_i = *features_i;
  
        ++clust_i;
        ++genenames_i;
      }
    }
    
    cpp11::named_arg _fc(fcname.c_str()); _fc = export_vec_to_rvec<cpp11::writable::doubles>(fc, el_count);
    cpp11::named_arg _p1(p1name.c_str()); _p1 = export_vec_to_rvec<cpp11::writable::doubles>(p1, el_count);
    cpp11::named_arg _p2(p2name.c_str()); _p2 = export_vec_to_rvec<cpp11::writable::doubles>(p2, el_count);
    cpp11::named_arg _cl("cluster"); _cl = clust;
    cpp11::named_arg _gn("gene"); _gn = genenames;
    return cpp11::writable::data_frame( {_cl, _gn, _fc, _p1, _p2} );
}

template <typename T>
cpp11::writable::data_frame export_fc_to_r_dataframe(
    std::vector<T> const & fc, std::string const & fcname,
    std::vector<T> const & p1, std::string const & p1name,
    std::vector<T> const & p2, std::string const & p2name,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    cpp11::strings const & features
) {
    return export_fc_to_r_dataframe(
        fc.data(), fcname,
        p1.data(), p1name,
        p2.data(), p2name,
        fc.size(), sorted_labels, features);
}


