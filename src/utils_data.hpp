#pragma once

// ------- function declaration
#include "cpp11/R.hpp"

#include "cpp11/r_vector.hpp"
#include "cpp11/matrix.hpp"
#include "cpp11/list.hpp"
#include "cpp11/integers.hpp"
#include "cpp11/doubles.hpp"
#include "cpp11/logicals.hpp"
#include "cpp11/data_frame.hpp"
#include "cpp11/strings.hpp"

#include <vector>



void copy_rmatrix_to_cppvector(
    cpp11::doubles_matrix<cpp11::by_column> const & _matrix, 
    std::vector<double> & mat);

void copy_rmatrix_to_cppvector(
    cpp11::logicals_matrix<cpp11::by_column> const & _matrix, 
    std::vector<bool> & mat);

void copy_rmatrix_to_cppvector(
    cpp11::doubles_matrix<cpp11::by_column> const & _matrix, 
    double * mat);

void copy_rmatrix_to_cppvector(
    cpp11::logicals_matrix<cpp11::by_column> const & _matrix, 
    bool * mat);


size_t copy_rvector_to_cppvector(cpp11::logicals const & _vector, 
    std::vector<bool> & vec, 
    size_t const & length = std::numeric_limits<size_t>::max(), 
    size_t const & offset = 0);

size_t copy_rvector_to_cppvector(cpp11::doubles const & _vector, 
    std::vector<double> & vec, 
    size_t const & length = std::numeric_limits<size_t>::max(), 
    size_t const & offset = 0);

size_t copy_rvector_to_cppvector(cpp11::integers const & _vector, 
    std::vector<int> & vec, 
    size_t const & length = std::numeric_limits<size_t>::max(), 
    size_t const & offset = 0);

size_t copy_rvector_to_cppvector(cpp11::doubles const & _vector, 
    std::vector<long> & vec, 
    size_t const & length = std::numeric_limits<size_t>::max(), 
    size_t const & offset = 0);

size_t copy_rvector_to_cppvector(cpp11::logicals const & _vector, 
    bool * vec, 
    size_t const & length = std::numeric_limits<size_t>::max(), 
    size_t const & offset = 0);

size_t copy_rvector_to_cppvector(cpp11::doubles const & _vector, 
    double * vec, 
    size_t const & length = std::numeric_limits<size_t>::max(), 
    size_t const & offset = 0);

size_t copy_rvector_to_cppvector(cpp11::integers const & _vector, 
    int * vec, 
    size_t const & length = std::numeric_limits<size_t>::max(), 
    size_t const & offset = 0);

size_t copy_rvector_to_cppvector(cpp11::doubles const & _vector, 
    long * vec, 
    size_t const & length = std::numeric_limits<size_t>::max(), 
    size_t const & offset = 0);


template <typename R_VEC, typename ITER>
R_VEC export_vec_to_rvec(ITER pv, size_t const & count);

// template <typename R_VEC, typename T>
// R_VEC export_vec_to_rvec(std::vector<T> const & pv);

template <typename R_MAT, typename ITER>
R_MAT export_vec_to_r_matrix(
    ITER pv, size_t const & nrow, size_t const & ncol
);

template <typename R_MAT, typename T>
R_MAT export_vec_to_r_matrix(
    std::vector<T> const & pv, size_t const & nrow, size_t const & ncol
);

template <typename T>
cpp11::writable::list export_fc_to_r_matrix(
    std::vector<T> const & fc, std::string const & fcname,
    std::vector<T> const & p1, std::string const & p1name,
    std::vector<T> const & p2, std::string const & p2name,
    std::vector<std::pair<int, size_t> > const & sorted_labels
);

template <typename ITER>
cpp11::writable::list export_fc_to_r_matrix(
    ITER fc, std::string const & fcname,
    ITER p1, std::string const & p1name,
    ITER p2, std::string const & p2name,
    size_t const & count,
    std::vector<std::pair<int, size_t> > const & sorted_labels
);

template <typename T>
cpp11::writable::data_frame export_vec_to_r_dataframe(
    std::vector<T> const & pv, std::string const & name,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    cpp11::strings const & features
);
template <typename ITER>
cpp11::writable::data_frame export_vec_to_r_dataframe(
    ITER pv, std::string const & name, size_t const & count,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    cpp11::strings const & features
);

template <typename T>
cpp11::writable::data_frame export_fc_to_r_dataframe(
    std::vector<T> const & fc, std::string const & fcname,
    std::vector<T> const & p1, std::string const & p1name,
    std::vector<T> const & p2, std::string const & p2name,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    cpp11::strings const & features
);

template <typename ITER>
cpp11::writable::data_frame export_fc_to_r_dataframe(
    ITER fc, std::string const & fcname,
    ITER p1, std::string const & p1name,
    ITER p2, std::string const & p2name,
    size_t const & count,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    cpp11::strings const & features
);
