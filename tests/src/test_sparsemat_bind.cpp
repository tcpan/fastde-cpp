/*
 * test_sparse.cpp
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


TEST_CASE( "cbind_iter", "[sparsemat]" ) {

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
	utils::read_hdf5_sparse_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc_cbind_csc.h5", "sparse_matrix", goldx, goldi, goldp, goldrows, goldcols, format);


	int repeats = 3;

	std::vector<std::vector<double>::const_iterator> ix(repeats);
	std::vector<std::vector<int>::const_iterator> ii(repeats);
	std::vector<std::vector<int>::const_iterator> ip(repeats);
	std::vector<int> vrows(repeats, rows);
	std::vector<int> vcols(repeats, cols);

	for (int v = 0; v < repeats; ++v) {
		ix[v] = x.cbegin();
		ii[v] = i.cbegin();
		ip[v] = p.cbegin();
	}

	std::vector<double> ox(x.size() * repeats);
	std::vector<int> oi(i.size() * repeats);
	std::vector<int> op(cols * repeats + 1, 0);

	csc_cbind(ix, ii, ip, vrows, vcols, ox.begin(), oi.begin(), op.begin(), 1);

	// compare
	REQUIRE(goldcols == cols * repeats);
	REQUIRE(goldrows == rows);
	REQUIRE(ox == goldx);
	REQUIRE(oi == goldi);
	REQUIRE(op == goldp);


	ox.clear(); ox.resize(x.size() * repeats);
	oi.clear(); oi.resize(i.size() * repeats);
	op.clear(); op.resize(cols * repeats + 1, 0);

	csc_cbind(ix, ii, ip, vrows, vcols, ox.begin(), oi.begin(), op.begin(), 4);

	// compare
	REQUIRE(goldcols == cols * repeats);
	REQUIRE(goldrows == rows);
	REQUIRE(ox == goldx);
	REQUIRE(oi == goldi);
	REQUIRE(op == goldp);

}


TEST_CASE( "cbind_vec", "[sparsemat]" ) {

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
	utils::read_hdf5_sparse_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc_cbind_csc.h5", "sparse_matrix", goldx, goldi, goldp, goldrows, goldcols, format);


	int repeats = 3;

	std::vector<std::vector<double>> ix(repeats);
	std::vector<std::vector<int>> ii(repeats);
	std::vector<std::vector<int>> ip(repeats);
	std::vector<int> vrows(repeats, rows);
	std::vector<int> vcols(repeats, cols);

	for (int v = 0; v < repeats; ++v) {
		ix[v] = x;
		ii[v] = i;
		ip[v] = p;
	}

	std::vector<double> ox(x.size() * repeats);
	std::vector<int> oi(i.size() * repeats);
	std::vector<int> op(cols * repeats + 1, 0);

	csc_cbind_vec(ix, ii, ip, vrows, vcols, ox, oi, op, 1);

	// compare
	REQUIRE(goldcols == cols * repeats);
	REQUIRE(goldrows == rows);
	REQUIRE(ox == goldx);
	REQUIRE(oi == goldi);
	REQUIRE(op == goldp);


	ox.clear(); ox.resize(x.size() * repeats);
	oi.clear(); oi.resize(i.size() * repeats);
	op.clear(); op.resize(cols * repeats + 1, 0);

	csc_cbind_vec(ix, ii, ip, vrows, vcols, ox, oi, op, 4);

	// compare
	REQUIRE(goldcols == cols * repeats);
	REQUIRE(goldrows == rows);
	REQUIRE(ox == goldx);
	REQUIRE(oi == goldi);
	REQUIRE(op == goldp);

}


TEST_CASE( "rbind_iter", "[sparsemat]" ) {

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
	utils::read_hdf5_sparse_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc_rbind_csc.h5", "sparse_matrix", goldx, goldi, goldp, goldrows, goldcols, format);


	int repeats = 3;

	std::vector<std::vector<double>::const_iterator> ix(repeats);
	std::vector<std::vector<int>::const_iterator> ii(repeats);
	std::vector<std::vector<int>::const_iterator> ip(repeats);
	std::vector<int> vrows(repeats, rows);
	std::vector<int> vcols(repeats, cols);

	for (int v = 0; v < repeats; ++v) {
		ix[v] = x.cbegin();
		ii[v] = i.cbegin();
		ip[v] = p.cbegin();
	}

	std::vector<double> ox(x.size() * repeats);
	std::vector<int> oi(i.size() * repeats);
	std::vector<int> op(cols + 1, 0);

	csc_rbind(ix, ii, ip, vrows, vcols, ox.begin(), oi.begin(), op.begin(), 1);

	// compare
	REQUIRE(goldcols == cols);
	REQUIRE(goldrows == rows * repeats);
	REQUIRE(ox.size() == goldx.size());
	REQUIRE(oi.size() == goldi.size());
	REQUIRE(op.size() == goldp.size());
	REQUIRE(ox == goldx);
	REQUIRE(oi == goldi);
	REQUIRE(op == goldp);


	ox.clear(); ox.resize(x.size() * repeats);
	oi.clear(); oi.resize(i.size() * repeats);
	op.clear(); op.resize(cols + 1, 0);

	csc_rbind(ix, ii, ip, vrows, vcols, ox.begin(), oi.begin(), op.begin(), 4);

	// compare
	REQUIRE(goldcols == cols);
	REQUIRE(goldrows == rows * repeats);
	REQUIRE(ox == goldx);
	REQUIRE(oi == goldi);
	REQUIRE(op == goldp);

}


TEST_CASE( "rbind_vec", "[sparsemat]" ) {

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
	utils::read_hdf5_sparse_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc_rbind_csc.h5", "sparse_matrix", goldx, goldi, goldp, goldrows, goldcols, format);


	int repeats = 3;

	std::vector<std::vector<double>> ix(repeats);
	std::vector<std::vector<int>> ii(repeats);
	std::vector<std::vector<int>> ip(repeats);
	std::vector<int> vrows(repeats, rows);
	std::vector<int> vcols(repeats, cols);

	for (int v = 0; v < repeats; ++v) {
		ix[v] = x;
		ii[v] = i;
		ip[v] = p;
	}

	std::vector<double> ox(x.size() * repeats);
	std::vector<int> oi(i.size() * repeats);
	std::vector<int> op(cols + 1, 0);

	csc_rbind_vec(ix, ii, ip, vrows, vcols, ox, oi, op, 1);

	// compare
	REQUIRE(goldcols == cols);
	REQUIRE(goldrows == rows * repeats);
	REQUIRE(ox.size() == goldx.size());
	REQUIRE(oi.size() == goldi.size());
	REQUIRE(op.size() == goldp.size());
	REQUIRE(ox == goldx);
	REQUIRE(oi == goldi);
	REQUIRE(op == goldp);


	ox.clear(); ox.resize(x.size() * repeats);
	oi.clear(); oi.resize(i.size() * repeats);
	op.clear(); op.resize(cols + 1, 0);

	csc_rbind_vec(ix, ii, ip, vrows, vcols, ox, oi, op, 4);

	// compare
	REQUIRE(goldcols == cols);
	REQUIRE(goldrows == rows * repeats);
	REQUIRE(ox == goldx);
	REQUIRE(oi == goldi);
	REQUIRE(op == goldp);
}

