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
#include "utils/HDF5Reader.hpp"
#include "utils/HDF5Writer.hpp"

namespace utils {

template <typename T, typename std::enable_if<std::is_integral<T>::value, int>::type dummy = 1 >
std::vector<T> make_random_vector(
	long const & seed, T const & rmin, T const & rmax,
	size_t const & count
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

	std::vector<T> vec(count);
    std::generate(begin(vec), end(vec), gen);

	return vec;
};

template <typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type dummy = 1 >
std::vector<T> make_random_vector(
	long const & seed, T const & rmin, T const & rmax,
	size_t const & count
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

	std::vector<T> vec(count);
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

	// average number of elements per column.
	size_t col_target = static_cast<size_t>(static_cast<double>(rows) * sparsity);
	
	std::random_device rd;
	std::mt19937 gen(rd());
	std::poisson_distribution<S> d(col_target);  // poisson, mean is col_target.

	// per column count, partial summed.
	for (size_t c = 1; c <= cols; ++c) {
		p[c] = p[c-1] + d(gen);
	}

	// total count is p[cols]
	// generate x
	std::vector<T> x = make_random_vector(seed, rmin, rmax, p[cols] );
	// generate row id
	std::vector<int> i = make_random_vector(seed, static_cast<int>(0), static_cast<int>(rows-1), p[cols]);
	// sort the row id within each columns

	for (size_t c = 0; c < cols; ++c) {
		std::sort(i.begin() + p[c], i.begin() + p[c+1]);
	}

	return std::make_tuple(x, i, p);
};



template <typename T, typename S>
std::vector<T> make_random_matrix(
	long const & seed, T const & rmin, T const & rmax,
	S const & rows, S const & cols
) {
	return make_random_vector(seed, rmin, rmax, rows * cols);
};


template <typename T, typename S>
void read_hdf5_sparse_matrix(std::string const & filename, 
	std::string const & datasetname,
	std::vector<T> & x, std::vector<int> & i, std::vector<S> & p, size_t & rows, size_t & cols,
	bool & is_transposed) {
	if (filename.length() > 0) {
		// write to file.  MPI enabled.  Not thread enabled.
		// FMT_ROOT_PRINT("name sizes: {}, {}\n", row_names.size(), col_names.size());
		// FMT_ROOT_PRINT("outputing matrix size: {}, {}\n", data.rows(), data.columns());
		utils::HDF5Reader reader(filename);
		// sequential writes only.
		ssize_t nr = 0, nc = 0, nz = 0;
		reader.getSparseMatrixSize(datasetname, nr, nc, nz);

		rows = nr;
		cols = nc;

		x.resize(nz);
		i.resize(nz);
		p.resize(cols + 1, 0);
		
		reader.loadSparseMatrixData(datasetname, nr, nc, nz, x.data(), i.data(), p.data());

	}
};

template <typename T>
void read_hdf5_vector(std::string const & filename, 
	std::string const & datasetname,
	std::vector<T> & out, size_t & count) {
	if (filename.length() > 0) {
		// write to file.  MPI enabled.  Not thread enabled.
		// FMT_ROOT_PRINT("name sizes: {}, {}\n", row_names.size(), col_names.size());
		// FMT_ROOT_PRINT("outputing matrix size: {}, {}\n", data.rows(), data.columns());
		utils::HDF5Reader reader(filename);
		// sequential writes only.

		ssize_t cnt;
		reader.getVectorSize(datasetname, cnt);

		out.resize(cnt);

		count = cnt;
		reader.loadVectorData(datasetname, count, out.data());

	}
};


template <typename T>
void read_hdf5_matrix(std::string const & filename, 
	std::string const & datasetname,
	std::vector<T> & data, size_t & rows, size_t & cols, bool & is_transposed) {
	if (filename.length() > 0) {
		// write to file.  MPI enabled.  Not thread enabled.
		// FMT_ROOT_PRINT("name sizes: {}, {}\n", row_names.size(), col_names.size());
		// FMT_ROOT_PRINT("outputing matrix size: {}, {}\n", data.rows(), data.columns());
		utils::HDF5Reader reader(filename);

		ssize_t nr, nc;
		reader.getMatrixSize(datasetname, nr, nc);

		data.resize(nr * nc);

		rows = nr;
		cols = nc;
		reader.loadMatrixData(datasetname, nr, nc, data.data(), nc * sizeof(T));

	}
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
	std::vector<T> const & data, 
	size_t const & rows, size_t const & cols, bool const & is_transposed = true) {
	if (filename.length() > 0) {
		// write to file.  MPI enabled.  Not thread enabled.
		// FMT_ROOT_PRINT("name sizes: {}, {}\n", row_names.size(), col_names.size());
		// FMT_ROOT_PRINT("outputing matrix size: {}, {}\n", data.rows(), data.columns());
		utils::HDF5Writer writer(filename);

		writer.storeMatrixData(datasetname, rows, cols, data.data(), cols * sizeof(T), false);

	}
};

}