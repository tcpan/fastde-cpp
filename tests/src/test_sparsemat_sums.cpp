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

TEST_CASE( "colsums_iter", "[sparsemat]" ) {

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
	size_t counts = 0;
	utils::read_hdf5_vector(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc_csum.h5", "vector", gold, counts);

	std::vector<double> out(cols);
	csc_colsums_iter(x.cbegin(), p.cbegin(), cols, out.begin(), 1);

	// compare
	REQUIRE(counts == cols);
    REQUIRE_THAT(out, Catch::Matchers::Approx(gold));

	out.clear(); out.resize(cols);
	csc_colsums_iter(x.cbegin(), p.cbegin(), cols, out.begin(), 4);

	// compare
	REQUIRE(counts == cols);
    REQUIRE_THAT(out, Catch::Matchers::Approx(gold));

}


TEST_CASE( "colsums_vec", "[sparsemat]" ) {

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
	size_t counts;
	utils::read_hdf5_vector(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc_csum.h5", "vector", gold, counts);

	std::vector<double> out(cols);
	csc_colsums_vec(x, p, cols, out, 1);

	// compare
	REQUIRE(counts == cols);
    REQUIRE_THAT(out, Catch::Matchers::Approx(gold));

	out.clear(); out.resize(cols);
	csc_colsums_vec(x, p, cols, out, 4);

	// compare
	REQUIRE(counts == cols);
    REQUIRE_THAT(out, Catch::Matchers::Approx(gold));
}


TEST_CASE( "rowsums_iter", "[sparsemat]" ) {

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
	size_t counts;
	utils::read_hdf5_vector(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc_rsum.h5", "vector", gold, counts);

	std::vector<double> out(rows);
	csc_rowsums_iter(x.cbegin(), i.cbegin(), rows, x.size(), out.begin(), 1);

	// compare
	REQUIRE(counts == rows);
    REQUIRE_THAT(out, Catch::Matchers::Approx(gold));

	out.clear(); out.resize(rows);
	csc_rowsums_iter(x.cbegin(), i.cbegin(), rows, x.size(), out.begin(), 4);

	// compare
	REQUIRE(counts == rows);
    REQUIRE_THAT(out, Catch::Matchers::Approx(gold));

}


TEST_CASE( "rowsums_vec", "[sparsemat]" ) {

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
	size_t counts;
	utils::read_hdf5_vector(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc_rsum.h5", "vector", gold, counts);

	std::vector<double> out(rows);
	csc_rowsums_vec(x, i, rows, out, 1);

	// compare
	REQUIRE(counts == rows);
    REQUIRE_THAT(out, Catch::Matchers::Approx(gold));

	out.clear(); out.resize(rows);
	csc_rowsums_vec(x, i, rows, out, 4);

	// compare
	REQUIRE(counts == rows);
    REQUIRE_THAT(out, Catch::Matchers::Approx(gold));
}



