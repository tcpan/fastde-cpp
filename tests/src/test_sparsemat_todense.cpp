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


TEST_CASE( "to_dense_iter", "[sparsemat]" ) {
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
	utils::read_hdf5_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc_to_dense_c.h5", "matrix", gold, goldrows, goldcols, rowMajor);

	// compute
	std::vector<double> out(rows * cols, 0);
	csc_to_dense_c(x.cbegin(), i.cbegin(), p.cbegin(), rows, cols, out.begin(), 1);
	// note:  the test_spmat_csc_to_dense_c.h5 (written by python) has weird data.  temp_spmat_csc_to_dense_c.h5 (written by c++) have correct data when using h5dump:  
	//  to dense: 0.8972061247 (computed) == 1.5560668779 (in h5 file),  transposed_to_dense:  0.6131512398 == 1.3604524348, and there are multiple entries in python version of array that are over 1 and under 1.
	// using just python to do round trip read with files (some val > 1) seems okay.   using c++ for roundtrip seems fine too.
	// data types are the same.  is the problem in both h5dump and my c++ matrix read code? probably python code.  problem is only with dense matrix, not sparse, or vec.

	// recap. python write is not correct, and c++ write seems correct.  but the exact issue and condition is not clear.
	// cause:   duplicate row ids in the same column in the randomly generated data.  When scipy.sparse.matrix.todense() is called, these value are added together, silently....
	utils::write_hdf5_matrix("temp_spmat_csc_to_dense_c.h5", "matrix", out, rows, cols, false);
	
	// compare
	REQUIRE(goldrows == rows);
	REQUIRE(goldcols == cols);
	REQUIRE(out.size() == gold.size());
	REQUIRE(out == gold);

	// compute
	out.clear(); out.resize(rows * cols, 0);
	csc_to_dense_c(x.cbegin(), i.cbegin(), p.cbegin(), rows, cols, out.begin(), 4);
		
	// compare
	REQUIRE(goldrows == rows);
	REQUIRE(goldcols == cols);
	REQUIRE(out == gold);
}


TEST_CASE( "to_dense_vec", "[sparsemat]" ) {
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

	utils::read_hdf5_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_csc_to_dense_c.h5", "matrix", gold, goldrows, goldcols, rowMajor);

	// compute
	std::vector<double> out(rows * cols);
	csc_to_dense_c_vec(x, i, p, rows, cols, out, 1);   // note that this output is column major.
		
	// compare
	REQUIRE(goldrows == rows);
	REQUIRE(goldcols == cols);
	REQUIRE(out.size() == gold.size());
	REQUIRE(out == gold);

	// compute
	out.clear(); out.resize(rows * cols);
	csc_to_dense_c_vec(x, i, p, rows, cols, out, 4);
		
	// compare
	REQUIRE(goldrows == rows);
	REQUIRE(goldcols == cols);
	REQUIRE(out == gold);
}


