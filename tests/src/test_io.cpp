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


TEST_CASE( "row_matrix", "[io]" ) {
	// read
	std::vector<double> gold;
	bool rowMajor;
	size_t rows, cols;
	utils::read_hdf5_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_matrix_r.h5", "matrix", gold, rows, cols, rowMajor);

	// write			
	utils::write_hdf5_matrix("temp_matrix_r.h5", "matrix", gold, rows, cols, rowMajor);

	// read
	std::vector<double> test;
	size_t rows2, cols2;
	bool rowMajor2;
	utils::read_hdf5_matrix("temp_matrix_r.h5", "matrix", test, rows2, cols2, rowMajor2);

	// compare
	REQUIRE(rows2 > 0);
	REQUIRE(cols2 > 0);
	REQUIRE(rowMajor == rowMajor2);
	REQUIRE(rowMajor == true);
	REQUIRE(test.size() > 0);
	REQUIRE(rows2 == rows);
	REQUIRE(cols2 == cols);
	REQUIRE(test == gold);
}


TEST_CASE( "col_matrix", "[io]" ) {
	// read
	std::vector<double> gold;
	bool rowMajor;
	size_t rows, cols;
	utils::read_hdf5_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_matrix_c.h5", "matrix", gold, rows, cols, rowMajor);

	// write			
	utils::write_hdf5_matrix("temp_matrix_c.h5", "matrix", gold, rows, cols, rowMajor);

	// read
	std::vector<double> test;
	size_t rows2, cols2;
	bool rowMajor2;
	utils::read_hdf5_matrix("temp_matrix_c.h5", "matrix", test, rows2, cols2, rowMajor2);

	// compare
	REQUIRE(rows2 > 0);
	REQUIRE(cols2 > 0);
	REQUIRE(rowMajor == rowMajor2);
	REQUIRE(rowMajor == false);
	REQUIRE(test.size() > 0);
	REQUIRE(rows2 == rows);
	REQUIRE(cols2 == cols);
	REQUIRE(test == gold);
}



TEST_CASE( "vector", "[io]" ) {
	// read
	std::vector<int> gold;
	size_t count;
	utils::read_hdf5_vector(std::string(FDE_TEST_DATA_DIR) + "/test_int_vec.h5", "vector", gold, count);

	// write			
	utils::write_hdf5_vector("temp_int_vec.h5", "vector", gold, count);

	// read
	std::vector<int> test;
	size_t count2;
	utils::read_hdf5_vector("temp_int_vec.h5", "vector", test, count2);

	// compare
	REQUIRE(count2 > 0);
	REQUIRE(test.size() > 0);
	REQUIRE(count2 == count);
	REQUIRE(test == gold);

	std::vector<double> dgold;
	utils::read_hdf5_vector(std::string(FDE_TEST_DATA_DIR) + "/test_double_vec.h5", "vector", dgold, count);

	// write			
	utils::write_hdf5_vector("temp_double_vec.h5", "vector", dgold, count);

	// read
	std::vector<double> dtest;
	utils::read_hdf5_vector("temp_double_vec.h5", "vector", dtest, count2);

	// compare
	REQUIRE(count2 == count);
	REQUIRE(dtest == dgold);
}

TEST_CASE( "sparsemat_csc", "[io]" ) {
	// read
	std::vector<double> goldx;
	std::vector<int> goldi;
	std::vector<int> goldp;
	size_t rows, cols;
	std::string format;
	utils::read_hdf5_sparse_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc.h5", "sparse_matrix", goldx, goldi, goldp, rows, cols, format);

	// write			
	utils::write_hdf5_sparse_matrix("temp_spmat_csc.h5", "sparse_matrix", goldx, goldi, goldp, rows, cols, format);




	// read
	std::vector<double> testx;
	std::vector<int> testi;
	std::vector<int> testp;
	size_t rows2, cols2;
	std::string format2;
	utils::read_hdf5_sparse_matrix("temp_spmat_csc.h5", "sparse_matrix", testx, testi, testp, rows2, cols2, format2);

	// compare
	REQUIRE(rows2 > 0);
	REQUIRE(cols2 > 0);
	REQUIRE(format == format2);
	REQUIRE(format == "csc");
	REQUIRE(testx.size() > 0);
	REQUIRE(testi.size() > 0);
	REQUIRE(testp.size() > 0);
	REQUIRE(testx.size() == goldx.size());
	REQUIRE(testi.size() == goldi.size());
	REQUIRE(testp.size() == goldp.size());
	REQUIRE(rows2 == rows);
	REQUIRE(cols2 == cols);
	REQUIRE(testx == goldx);
	REQUIRE(testi == goldi);
	REQUIRE(testp == goldp);
}

TEST_CASE( "sparsemat_csr", "[io]" ) {
	// read
	std::vector<double> goldx;
	std::vector<int> goldi;
	std::vector<int> goldp;
	size_t rows, cols;
	std::string format;
	utils::read_hdf5_sparse_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csr.h5", "sparse_matrix", goldx, goldi, goldp, rows, cols, format);

	// write			
	utils::write_hdf5_sparse_matrix("temp_spmat_csr.h5", "sparse_matrix", goldx, goldi, goldp, rows, cols, format);

	// read
	std::vector<double> testx;
	std::vector<int> testi;
	std::vector<int> testp;
	size_t rows2, cols2;
	std::string format2;
	utils::read_hdf5_sparse_matrix("temp_spmat_csr.h5", "sparse_matrix", testx, testi, testp, rows2, cols2, format2);

	// compare
	REQUIRE(rows2 > 0);
	REQUIRE(cols2 > 0);
	REQUIRE(format == "csr");
	REQUIRE(format == format2);
	REQUIRE(testx.size() > 0);
	REQUIRE(testi.size() > 0);
	REQUIRE(testp.size() > 0);
	REQUIRE(testx.size() == goldx.size());
	REQUIRE(testi.size() == goldi.size());
	REQUIRE(testp.size() == goldp.size());
	REQUIRE(rows2 == rows);
	REQUIRE(cols2 == cols);
	REQUIRE(testx == goldx);
	REQUIRE(testi == goldi);
	REQUIRE(testp == goldp);
}

