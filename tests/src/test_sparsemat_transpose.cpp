/*
 * test_sparsemat_transpose.cpp
 *
 *  Created on: June 29, 2020
 *      Author: Tony Pan
 *		  School of Computational Science & Engineering
 *		  Georgia Institute of Technology, USA
 */

#include <sstream>
#include <string>
#include <iostream>
#include <fstream>

#include "utils/CLIParserCommon.hpp"
#include "utils/benchmark.hpp"
#include "utils/matrix_io.hpp"
#include "fastde/benchmark_utils.hpp"
#include "fastde/sparsemat.tpp"


#ifdef USE_OPENMP
#include <omp.h>
#endif


#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

TEST_CASE( "transpose_iter", "[sparsemat]" ) {
	// source
	std::vector<double> x;
	std::vector<int> i;
	std::vector<int> p;
	size_t rows, cols;
	std::string format;
	utils::read_hdf5_sparse_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc.h5", "sparse_matrix", x, i, p, rows, cols, format);

	REQUIRE(rows > 0);
	REQUIRE(cols > 0);
	REQUIRE(format == "csc");
	REQUIRE(x.size() > 0);
	REQUIRE(i.size() > 0);
	REQUIRE(p.size() > 0);

	// gold
	std::vector<double> goldx;
	std::vector<int> goldi;
	std::vector<int> goldp;
	size_t goldrows, goldcols;
	std::string goldformat;
	utils::read_hdf5_sparse_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc_transposed_csc.h5", "sparse_matrix", goldx, goldi, goldp, goldrows, goldcols, goldformat);

	REQUIRE(goldrows > 0);
	REQUIRE(goldcols > 0);
	REQUIRE(goldformat == "csc");
	REQUIRE(goldx.size() > 0);
	REQUIRE(goldi.size() > 0);
	REQUIRE(goldp.size() > 0);

	// compute
	std::vector<double> ox(x.size());
	std::vector<int> oi(i.size());
	std::vector<int> op(rows + 1);
	csc_transpose_csc(x.cbegin(), i.cbegin(), p.cbegin(), rows, cols, ox.begin(), oi.begin(), op.begin(), 1);
		
	// compare
	REQUIRE(goldrows == cols);
	REQUIRE(goldcols == rows);
	REQUIRE(ox == goldx);
	REQUIRE(oi == goldi);
	REQUIRE(op == goldp);

	// compare
	ox.clear(); ox.resize(x.size());
	oi.clear(); oi.resize(i.size());
	op.clear(); op.resize(rows + 1);
	csc_transpose_csc(x.cbegin(), i.cbegin(), p.cbegin(), rows, cols, ox.begin(), oi.begin(), op.begin(), 4);

	REQUIRE(goldrows == cols);
	REQUIRE(goldcols == rows);
	REQUIRE(ox == goldx);
	REQUIRE(oi == goldi);
	REQUIRE(op == goldp);

}


TEST_CASE( "transpose_vec", "[sparsemat]" ) {
	// source
	std::vector<double> x;
	std::vector<int> i;
	std::vector<int> p;
	size_t rows, cols;
	std::string format;
	utils::read_hdf5_sparse_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc.h5", "sparse_matrix", x, i, p, rows, cols, format);

	REQUIRE(rows > 0);
	REQUIRE(cols > 0);
	REQUIRE(format == "csc");
	REQUIRE(x.size() > 0);
	REQUIRE(i.size() > 0);
	REQUIRE(p.size() > 0);

	// gold
	std::vector<double> goldx;
	std::vector<int> goldi;
	std::vector<int> goldp;
	size_t goldrows, goldcols;
	utils::read_hdf5_sparse_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc_transposed_csc.h5", "sparse_matrix", goldx, goldi, goldp, goldrows, goldcols, format);

	// compute
	std::vector<double> ox(x.size());
	std::vector<int> oi(i.size());
	std::vector<int> op(rows + 1);
	csc_transpose_csc_vec(x, i, p, rows, cols, ox, oi, op, 1);
		
	// compare
	REQUIRE(goldrows == cols);
	REQUIRE(goldcols == rows);
	REQUIRE(ox == goldx);
	REQUIRE(oi == goldi);
	REQUIRE(op == goldp);

	// compare
	ox.clear(); ox.resize(x.size());
	oi.clear(); oi.resize(i.size());
	op.clear(); op.resize(rows + 1);
	csc_transpose_csc_vec(x, i, p, rows, cols, ox, oi, op, 4);

	REQUIRE(goldrows == cols);
	REQUIRE(goldcols == rows);
	REQUIRE(ox == goldx);
	REQUIRE(oi == goldi);
	REQUIRE(op == goldp);

}


TEST_CASE( "transpose_to_dense_iter", "[sparsemat]" ) {
	// source
	std::vector<double> x;
	std::vector<int> i;
	std::vector<int> p;
	size_t rows, cols;
	std::string format;
	utils::read_hdf5_sparse_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc.h5", "sparse_matrix", x, i, p, rows, cols, format);

	REQUIRE(rows > 0);
	REQUIRE(cols > 0);
	REQUIRE(format == "csc");
	REQUIRE(x.size() > 0);
	REQUIRE(i.size() > 0);
	REQUIRE(p.size() > 0);

	// gold
	std::vector<double> gold;
	size_t goldrows, goldcols;
	bool rowMajor;
	utils::read_hdf5_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc_transpose_to_dense_c.h5", "matrix", gold, goldrows, goldcols, rowMajor);

	// compute
	std::vector<double> out(rows * cols);
	csc_to_dense_transposed_c(x.cbegin(), i.cbegin(), p.cbegin(), rows, cols, out.begin(), 1);
		
	// compare
	REQUIRE(goldrows == cols);
	REQUIRE(goldcols == rows);
	REQUIRE(out.size() == gold.size());
	REQUIRE(out == gold);

	// compute
	out.clear(); out.resize(rows * cols);
	csc_to_dense_transposed_c(x.cbegin(), i.cbegin(), p.cbegin(), rows, cols, out.begin(), 4);
		
	// compare
	REQUIRE(goldrows == cols);
	REQUIRE(goldcols == rows);
	REQUIRE(out == gold);
}


TEST_CASE( "transpose_to_dense_vec", "[sparsemat]" ) {
	// source
	std::vector<double> x;
	std::vector<int> i;
	std::vector<int> p;
	size_t rows, cols;
	std::string format;
	utils::read_hdf5_sparse_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc.h5", "sparse_matrix", x, i, p, rows, cols, format);

	REQUIRE(rows > 0);
	REQUIRE(cols > 0);
	REQUIRE(format == "csc");
	REQUIRE(x.size() > 0);
	REQUIRE(i.size() > 0);
	REQUIRE(p.size() > 0);

	// gold
	std::vector<double> gold;
	size_t goldrows, goldcols;
	bool rowMajor;
	utils::read_hdf5_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc_transpose_to_dense_c.h5", "matrix", gold, goldrows, goldcols, rowMajor);

	// compute
	std::vector<double> out(rows * cols);
	csc_to_dense_transposed_c_vec(x, i, p, rows, cols, out, 1);
		
	// compare
	REQUIRE(goldrows == cols);
	REQUIRE(goldcols == rows);
	REQUIRE(out.size() == gold.size());
	REQUIRE(out == gold);

	// compute
	out.clear(); out.resize(rows * cols);
	csc_to_dense_transposed_c_vec(x, i, p, rows, cols, out, 4);
		
	// compare
	REQUIRE(goldrows == cols);
	REQUIRE(goldcols == rows);
	REQUIRE(out == gold);
}

