/*
 * test_normalize.cpp
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
#include "fastde/normalize.tpp"


#ifdef USE_OPENMP
#include <omp.h>
#endif


#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

TEST_CASE( "lognorm_iter", "[normalize]" ) {

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
	utils::read_hdf5_sparse_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc_lognorm.h5", "sparse_matrix", goldx, goldi, goldp, goldrows, goldcols, goldformat);

	REQUIRE(goldrows == rows);
	REQUIRE(goldcols == cols);
	REQUIRE(goldi == i);
	REQUIRE(goldp == p);
	REQUIRE(goldformat == format);

	std::vector<double> ox(x.size(), 0);
	csc_log_normalize_iter(x.cbegin(), p.cbegin(), cols, 0.5, ox.data(), 1);

	// compare
    REQUIRE_THAT(ox, Catch::Matchers::Approx(goldx));

	ox.clear(); ox.resize(x.size(), 0);
	csc_log_normalize_iter(x.cbegin(), p.cbegin(), cols, 0.5, ox.data(), 4);

	// compare
    REQUIRE_THAT(ox, Catch::Matchers::Approx(goldx));

}



TEST_CASE( "lognorm_vec", "[normalize]" ) {

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
	utils::read_hdf5_sparse_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc_lognorm.h5", "sparse_matrix", goldx, goldi, goldp, goldrows, goldcols, goldformat);

	REQUIRE(goldrows == rows);
	REQUIRE(goldcols == cols);
	REQUIRE(goldi == i);
	REQUIRE(goldp == p);
	REQUIRE(goldformat == format);

	std::vector<double> ox(x.size(), 0);
	csc_log_normalize_vec(x, p, cols, 0.5, ox, 1);

	// compare
    REQUIRE_THAT(ox, Catch::Matchers::Approx(goldx));

	ox.clear(); ox.resize(x.size(), 0);
	csc_log_normalize_vec(x, p, cols, 0.5, ox, 4);

	// compare
    REQUIRE_THAT(ox, Catch::Matchers::Approx(goldx));

}



TEST_CASE( "clr_col_iter", "[normalize]" ) {

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
	utils::read_hdf5_sparse_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc_clr.h5", "sparse_matrix", goldx, goldi, goldp, goldrows, goldcols, goldformat);

	REQUIRE(goldrows == rows);
	REQUIRE(goldcols == cols);
	REQUIRE(goldi == i);
	REQUIRE(goldp == p);
	REQUIRE(goldformat == format);

	std::vector<double> ox(x.size(), 0);
	csc_clr_cols_iter(x.cbegin(), i.cbegin(), p.cbegin(), rows, cols, ox.data(), 1);

	// compare
    REQUIRE_THAT(ox, Catch::Matchers::Approx(goldx));

	ox.clear(); ox.resize(x.size(), 0);
	csc_clr_cols_iter(x.cbegin(), i.cbegin(), p.cbegin(), rows, cols, ox.data(), 4);

	// compare
    REQUIRE_THAT(ox, Catch::Matchers::Approx(goldx));

}



TEST_CASE( "clr_col_vec", "[normalize]" ) {

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
	utils::read_hdf5_sparse_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc_clr.h5", "sparse_matrix", goldx, goldi, goldp, goldrows, goldcols, goldformat);

	REQUIRE(goldrows == rows);
	REQUIRE(goldcols == cols);
	REQUIRE(goldi == i);
	REQUIRE(goldp == p);
	REQUIRE(goldformat == format);

	std::vector<double> ox(x.size(), 0);
	csc_clr_cols_vec(x, i, p, rows, cols, ox, 1);

	// compare
    REQUIRE_THAT(ox, Catch::Matchers::Approx(goldx));

	ox.clear(); ox.resize(x.size(), 0);
	csc_clr_cols_vec(x, i, p, rows, cols, ox, 4);

	// compare
    REQUIRE_THAT(ox, Catch::Matchers::Approx(goldx));

}



TEST_CASE( "clr_row_iter", "[normalize]" ) {

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
	utils::read_hdf5_sparse_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc_clr_r.h5", "sparse_matrix", goldx, goldi, goldp, goldrows, goldcols, goldformat);

	REQUIRE(goldrows == rows);
	REQUIRE(goldcols == cols);
	REQUIRE(goldi == i);
	REQUIRE(goldp == p);
	REQUIRE(goldformat == format);

	std::vector<double> ox(x.size(), 0);
	csc_clr_rows_iter(x.cbegin(), i.cbegin(), p.cbegin(), rows, cols, ox.data(), 1);

	// compare
    REQUIRE_THAT(ox, Catch::Matchers::Approx(goldx));

	ox.clear(); ox.resize(x.size(), 0);
	csc_clr_rows_iter(x.cbegin(), i.cbegin(), p.cbegin(), rows, cols, ox.data(), 4);

	// compare
    REQUIRE_THAT(ox, Catch::Matchers::Approx(goldx));

}



TEST_CASE( "clr_row_vec", "[normalize]" ) {

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
	utils::read_hdf5_sparse_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc_clr_r.h5", "sparse_matrix", goldx, goldi, goldp, goldrows, goldcols, goldformat);

	REQUIRE(goldrows == rows);
	REQUIRE(goldcols == cols);
	REQUIRE(goldi == i);
	REQUIRE(goldp == p);
	REQUIRE(goldformat == format);

	std::vector<double> ox(x.size(), 0);
	csc_clr_rows_vec(x, i, p, rows, cols, ox, 1);

	// compare
    REQUIRE_THAT(ox, Catch::Matchers::Approx(goldx));

	ox.clear(); ox.resize(x.size(), 0);
	csc_clr_rows_vec(x, i, p, rows, cols, ox, 4);

	// compare
    REQUIRE_THAT(ox, Catch::Matchers::Approx(goldx));

}



TEST_CASE( "relcount_iter", "[normalize]" ) {

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
	utils::read_hdf5_sparse_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc_relative_count.h5", "sparse_matrix", goldx, goldi, goldp, goldrows, goldcols, goldformat);

	REQUIRE(goldrows == rows);
	REQUIRE(goldcols == cols);
	REQUIRE(goldi == i);
	REQUIRE(goldp == p);
	REQUIRE(goldformat == format);

	std::vector<double> ox(x.size(), 0);
	csc_relative_count_iter(x.cbegin(), p.cbegin(), cols, 0.5, ox.data(), 1);

	// compare
    REQUIRE_THAT(ox, Catch::Matchers::Approx(goldx));

	ox.clear(); ox.resize(x.size(), 0);
	csc_relative_count_iter(x.cbegin(), p.cbegin(), cols, 0.5, ox.data(), 4);

	// compare
    REQUIRE_THAT(ox, Catch::Matchers::Approx(goldx));
	
}



TEST_CASE( "relcount_vec", "[normalize]" ) {

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
	utils::read_hdf5_sparse_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc_relative_count.h5", "sparse_matrix", goldx, goldi, goldp, goldrows, goldcols, goldformat);

	REQUIRE(goldrows == rows);
	REQUIRE(goldcols == cols);
	REQUIRE(goldi == i);
	REQUIRE(goldp == p);
	REQUIRE(goldformat == format);

	std::vector<double> ox(x.size(), 0);
	csc_relative_count_vec(x, p, cols, 0.5, ox, 1);

	// compare
    REQUIRE_THAT(ox, Catch::Matchers::Approx(goldx));

	ox.clear(); ox.resize(x.size(), 0);
	csc_relative_count_vec(x, p, cols, 0.5, ox, 4);

	// compare
    REQUIRE_THAT(ox, Catch::Matchers::Approx(goldx));

}

