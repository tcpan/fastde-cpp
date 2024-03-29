/*
 * test_ttest.cpp
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
#include "fastde/benchmark_utils.tpp"
#include "fastde/cluster_utils.tpp"
#include "fastde/sparsemat.tpp"
#include "fastde/ttest.tpp"


#ifdef USE_OPENMP
#include <omp.h>
#endif


#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>


TEST_CASE( "ttest_sparse_iter", "[ttest]" ) {

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

	// labels
	std::vector<int> labels;
	size_t nlabels;
	utils::read_hdf5_vector(std::string(FDE_TEST_DATA_DIR) + "/test_label_vec.h5", "vector", labels, nlabels);

	std::unordered_set<int> unique;
	for (auto l : labels) {
		unique.insert(l);
	}
	size_t nuniq = unique.size();

	// gold
	std::vector<double> goldpv;
	size_t goldrows, goldcols;
	bool rowMajor;
	utils::read_hdf5_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_ttest_pval_matrix_c.h5", "matrix", goldpv, goldrows, goldcols, rowMajor);


	REQUIRE(goldcols == cols);
	REQUIRE(goldrows == nuniq);
	REQUIRE(rowMajor == false);

	std::vector<double> opv(nuniq * cols, 0);
	std::vector<std::pair<int, size_t> > cl_counts;

	omp_sparse_ttest(x.cbegin(), i.cbegin(), p.cbegin(), rows, cols, labels.cbegin(), 
		2, true, opv, cl_counts, 1);

	// compare
    REQUIRE_THAT(opv, Catch::Matchers::Approx(goldpv));

	opv.clear(); opv.resize(nuniq * cols, 0);
	cl_counts.clear();
	omp_sparse_ttest(x.cbegin(), i.cbegin(), p.cbegin(), rows, cols, labels.cbegin(), 
		2, true, opv, cl_counts, 4);

	// compare
    REQUIRE_THAT(opv, Catch::Matchers::Approx(goldpv));

}




TEST_CASE( "ttest_pseudosparse_iter", "[ttest]" ) {

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


	// labels
	std::vector<int> labels;
	size_t nlabels;
	utils::read_hdf5_vector(std::string(FDE_TEST_DATA_DIR) + "/test_label_vec.h5", "vector", labels, nlabels);

	std::unordered_set<int> unique;
	for (auto l : labels) {
		unique.insert(l);
	}
	size_t nuniq = unique.size();

	// gold
	std::vector<double> goldpv;
	size_t goldrows, goldcols;
	bool rowMajor;
	utils::read_hdf5_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_ttest_pval_matrix_c.h5", "matrix", goldpv, goldrows, goldcols, rowMajor);


	REQUIRE(goldcols == cols);
	REQUIRE(goldrows == nuniq);
	REQUIRE(rowMajor == false);


	// convert to dense
	std::vector<double> dense(rows * cols, 0);
	csc_to_dense_c_vec(x, i, p, rows, cols, dense, 1);


	std::vector<double> opv(nuniq * cols, 0);
	std::vector<std::pair<int, size_t> > cl_counts;
	omp_dense_ttest(dense.cbegin(), rows, cols, labels.cbegin(), 2, true, opv, cl_counts, 1);

	// compare
    REQUIRE_THAT(opv, Catch::Matchers::Approx(goldpv));


	opv.clear(); opv.resize(nuniq * cols, 0);
	cl_counts.clear();
	omp_dense_ttest(dense.cbegin(), rows, cols, labels.cbegin(), 2, true, opv, cl_counts, 4);

	// compare
    REQUIRE_THAT(opv, Catch::Matchers::Approx(goldpv));

}
