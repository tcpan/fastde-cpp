#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <algorithm>
#include <iterator>
#include <type_traits>
#include <random>

#ifdef USE_MPI
#include <mpi.h>
#endif

// #include "utils/HDF5MatrixWriter.hpp"
#include "utils/HDF5Writer.hpp"

namespace utils {

template <typename T, typename std::enable_if<std::is_integral<T>::value, int>::type dummy = 1 >
std::vector<T> make_random_vector(
	long const & seed, T const & rmin, T const & rmax,
	size_t const & cols
) {
	// allocate input.

    // First create an instance of an engine.
    std::random_device rnd_device;
    // Specify the engine and distribution.
    std::mt19937 mersenne_engine { rnd_device() };  // Generates random integers
    std::uniform_int_distribution<T> dist(rmin, rmax);
    
    auto gen = [&dist, &mersenne_engine](){
                   return dist(mersenne_engine);
               };

	std::vector<T> vec(cols);
    std::generate(begin(vec), end(vec), gen);

	return vec;
};

template <typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type dummy = 1 >
std::vector<T> make_random_vector(
	long const & seed, T const & rmin, T const & rmax,
	size_t const & cols
) {
	// allocate input.

    // First create an instance of an engine.
    std::random_device rnd_device;
    // Specify the engine and distribution.
    std::mt19937 mersenne_engine { rnd_device() };  // Generates random integers
    std::uniform_real_distribution<T> dist(rmin, rmax);
    
    auto gen = [&dist, &mersenne_engine](){
                   return dist(mersenne_engine);
               };

	std::vector<T> vec(cols);
    std::generate(begin(vec), end(vec), gen);

	return vec;
};

// produce csc file.
template <typename T, typename S, typename std::enable_if<std::is_integral<S>::value, int>::type dummy = 1 >
std::tuple<std::vector<T>, std::vector<int>, std::vector<S>> make_random_sparse_matrix(
	long const & seed, T const & rmin, T const & rmax,
	size_t const & rows, size_t const & cols, double const & sparsity
) {

	// first figure out how many per row.
	std::vector<S> p(cols+1, 0);

	size_t col_target = static_cast<size_t>(static_cast<double>(rows) * sparsity);
	
	std::random_device rd;
	std::mt19937 gen(rd());
	std::poisson_distribution<S> d(col_target);

	for (size_t c = 1; c <= cols; ++c) {
		p[c] = p[c-1] + d(gen);
	}

	// generate x
	std::vector<T> x = make_random_vector(seed, rmin, rmax, rows * cols);
	// generate row id
	std::vector<int> i = make_random_vector(seed, static_cast<int>(0), static_cast<int>(rows-1), rows*cols);
	// sort the row id within each columns

	for (size_t c = 0; c < cols; ++c) {
		std::sort(i.begin() + p[c], i.begin() + p[c+1]);
	}

	return std::make_tuple(x, i, p);
};


template <typename T, typename S>
void write_hdf5_sparse_matrix(std::string const & filename, 
	std::string const & datasetname,
	std::vector<T> const & x, std::vector<int> const & i, std::vector<S> const & p, size_t const & rows, size_t const & cols,
	bool const & is_transposed = true) {
	if (filename.length() > 0) {
		// write to file.  MPI enabled.  Not thread enabled.
		// FMT_ROOT_PRINT("name sizes: {}, {}\n", row_names.size(), col_names.size());
		// FMT_ROOT_PRINT("outputing matrix size: {}, {}\n", data.rows(), data.columns());
		utils::HDF5Writer writer(filename);
		// sequential writes only.
		writer.storeSparseMatrixData(datasetname, rows, cols, x, i, p, is_transposed);
	}
};

template <typename T>
void write_hdf5_vector(std::string const & filename, 
	std::string const & datasetname,
	std::vector<T> const & x, size_t const & count) {
	if (filename.length() > 0) {
		// write to file.  MPI enabled.  Not thread enabled.
		// FMT_ROOT_PRINT("name sizes: {}, {}\n", row_names.size(), col_names.size());
		// FMT_ROOT_PRINT("outputing matrix size: {}, {}\n", data.rows(), data.columns());
		utils::HDF5Writer writer(filename);
		// sequential writes only.
		writer.storeVectorData(datasetname, count, x);
	}
};


template <typename T>
void write_hdf5_matrix(std::string const & filename, 
	std::string const & datasetname,
	std::vector<T> const & data, size_t const & rows, size_t const & cols, bool const & is_transposed = true) {
	if (filename.length() > 0) {
		// write to file.  MPI enabled.  Not thread enabled.
		// FMT_ROOT_PRINT("name sizes: {}, {}\n", row_names.size(), col_names.size());
		// FMT_ROOT_PRINT("outputing matrix size: {}, {}\n", data.rows(), data.columns());
		utils::HDF5Writer writer(filename);

		writer.storeMatrixData(datasetname, rows, cols, data.data(), cols * sizeof(T), false);

	}
};

}