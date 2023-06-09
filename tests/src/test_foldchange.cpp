/*
 * test_foldchange.cpp
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
#include "fastde/foldchange.tpp"


#ifdef USE_OPENMP
#include <omp.h>
#endif


#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>


TEST_CASE( "foldchange_sparse_iter", "[fc]" ) {

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
	std::vector<double> goldfc;
	std::vector<double> goldp1;
	std::vector<double> goldp2;
	size_t goldrows, goldcols;
	bool rowMajor;
	utils::read_hdf5_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_foldchange_matrix_c.h5", "matrix", goldfc, goldrows, goldcols, rowMajor);
	REQUIRE(goldcols == cols);
	REQUIRE(goldrows == nuniq);
	REQUIRE(rowMajor == false);

	utils::read_hdf5_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_foldchange_p1_matrix_c.h5", "matrix", goldp1, goldrows, goldcols, rowMajor);
	REQUIRE(goldcols == cols);
	REQUIRE(goldrows == nuniq);
	REQUIRE(rowMajor == false);

	utils::read_hdf5_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_foldchange_p2_matrix_c.h5", "matrix", goldp2, goldrows, goldcols, rowMajor);
	REQUIRE(goldcols == cols);
	REQUIRE(goldrows == nuniq);
	REQUIRE(rowMajor == false);

	std::vector<double> ofc(nuniq * cols, 0);
	std::vector<double> op1(nuniq * cols, 0);
	std::vector<double> op2(nuniq * cols, 0);
	std::vector<std::pair<int, size_t> > cl_counts;

	omp_sparse_foldchange(x.cbegin(), i.cbegin(), p.cbegin(), rows, cols, 
		labels.cbegin(), 
		true, "fc", false, 0.25, true, 2, true, ofc, op1, op2, cl_counts, 1);

	// compare
    REQUIRE_THAT(ofc, Catch::Matchers::Approx(goldfc));
    REQUIRE_THAT(op1, Catch::Matchers::Approx(goldp1));
    REQUIRE_THAT(op2, Catch::Matchers::Approx(goldp2));

	ofc.clear(); ofc.resize(nuniq * cols, 0);
	op1.clear(); op1.resize(nuniq * cols, 0);
	op2.clear(); op2.resize(nuniq * cols, 0);
	cl_counts.clear();
	omp_sparse_foldchange(x.cbegin(), i.cbegin(), p.cbegin(), rows, cols,
			labels.cbegin(), 
		true, "fc", false, 0.25, true, 2, true, ofc, op1, op2, cl_counts, 4);


	// compare
    REQUIRE_THAT(ofc, Catch::Matchers::Approx(goldfc));
    REQUIRE_THAT(op1, Catch::Matchers::Approx(goldp1));
    REQUIRE_THAT(op2, Catch::Matchers::Approx(goldp2));

}





TEST_CASE( "foldchange_pseudosparse_iter", "[fc]" ) {

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
	std::vector<double> goldfc;
	std::vector<double> goldp1;
	std::vector<double> goldp2;
	size_t goldrows, goldcols;
	bool rowMajor;
	utils::read_hdf5_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_foldchange_matrix_c.h5", "matrix", goldfc, goldrows, goldcols, rowMajor);
	REQUIRE(goldcols == cols);
	REQUIRE(goldrows == nuniq);
	REQUIRE(rowMajor == false);

	utils::read_hdf5_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_foldchange_p1_matrix_c.h5", "matrix", goldp1, goldrows, goldcols, rowMajor);
	REQUIRE(goldcols == cols);
	REQUIRE(goldrows == nuniq);
	REQUIRE(rowMajor == false);

	utils::read_hdf5_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_foldchange_p2_matrix_c.h5", "matrix", goldp2, goldrows, goldcols, rowMajor);
	REQUIRE(goldcols == cols);
	REQUIRE(goldrows == nuniq);
	REQUIRE(rowMajor == false);


	REQUIRE(goldcols == cols);
	REQUIRE(goldrows == nuniq);
	REQUIRE(rowMajor == false);


	// convert to dense
	std::vector<double> dense(rows * cols, 0);
	csc_to_dense_c_vec(x, i, p, rows, cols, dense, 1);


	std::vector<double> ofc(nuniq * cols, 0);
	std::vector<double> op1(nuniq * cols, 0);
	std::vector<double> op2(nuniq * cols, 0);
	std::vector<std::pair<int, size_t> > cl_counts;
	omp_dense_foldchange(dense.cbegin(), rows, cols, labels.cbegin(), 
		true, "fc", false, 0.25, true, 2, true, ofc, op1, op2, cl_counts, 1);

	// compare
    REQUIRE_THAT(ofc, Catch::Matchers::Approx(goldfc));
    REQUIRE_THAT(op1, Catch::Matchers::Approx(goldp1));
    REQUIRE_THAT(op2, Catch::Matchers::Approx(goldp2));


	ofc.clear(); ofc.resize(nuniq * cols, 0);
	op1.clear(); op1.resize(nuniq * cols, 0);
	op2.clear(); op2.resize(nuniq * cols, 0);
	cl_counts.clear();
	omp_dense_foldchange(dense.cbegin(), rows, cols, labels.cbegin(),
		true, "fc", false, 0.25, true, 2, true, ofc, op1, op2, cl_counts, 1);

	// compare
    REQUIRE_THAT(ofc, Catch::Matchers::Approx(goldfc));
    REQUIRE_THAT(op1, Catch::Matchers::Approx(goldp1));
    REQUIRE_THAT(op2, Catch::Matchers::Approx(goldp2));

}
