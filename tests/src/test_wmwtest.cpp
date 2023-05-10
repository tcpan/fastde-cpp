/*
 * test_wmwtest.cpp
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
#include "fastde/wmwtest.tpp"


#ifdef USE_OPENMP
#include <omp.h>
#endif


#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>




TEST_CASE( "wmwtest_sparse_vec", "[wilcox]" ) {

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

	// count number of unique labels
	std::vector<std::pair<int, size_t> > cl_counts;
	count_clusters_vec(labels, labels.size(), cl_counts, 1);
	size_t nuniq = cl_counts.size();


	// gold
	std::vector<double> goldpv;
	size_t goldrows, goldcols;
	bool rowMajor;
	utils::read_hdf5_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_wilcox_pval_matrix_c.h5", "matrix", goldpv, goldrows, goldcols, rowMajor);

	REQUIRE(goldcols == cols);
	REQUIRE(goldrows == nuniq);
	REQUIRE(rowMajor == false);  // column major data.

	std::vector<double> opv(nuniq * cols, 0);
	// output should have cols rows and nuniq cols, and in row major.
	csc_wmw_vecsc(x, i, p, rows, cols, labels, 2, true, opv, cl_counts, 1);

	// compare
    REQUIRE_THAT(opv, Catch::Matchers::Approx(goldpv));  

	opv.clear(); opv.resize(nuniq * cols, 0);
	csc_wmw_vecsc(x, i, p, rows, cols, labels, 2, true, opv, cl_counts, 4);

	// compare
    REQUIRE_THAT(opv, Catch::Matchers::Approx(goldpv));


	// convert gold:
	// size_t offset = 0;
	// std::vector<std::vector<double>> gold(cols);
	// std::vector<std::vector<double>> pvm(cols);
	// for (size_t c = 0; c < cols; ++c) {
	// 	gold[c].resize(nuniq, 0);
	// 	for (size_t r = 0; r < nuniq; ++r, ++offset) {
	// 		gold[c][r] = goldpv[offset];
	// 	}

	// 	pvm[c].resize(nuniq, 0);
	// }

	// // test with vector of vector.
	// csc_wmw_matc(x, i, p, rows, cols, labels, 2, true, pvm, cl_counts, 1);

	// for (size_t c = 0; c < cols; ++c) {
	//     REQUIRE_THAT(pvm[c], Catch::Matchers::Approx(gold[c]));
	// }



	// for (size_t c = 0; c < cols; ++c) {
	// 	pvm[c].clear();
	// 	pvm[c].resize(nuniq, 0);
	// }

	// // test with vector of vector.
	// csc_wmw_matc(x, i, p, rows, cols, labels, 2, true, pvm, cl_counts, 4);

	// for (size_t c = 0; c < cols; ++c) {
	//     REQUIRE_THAT(pvm[c], Catch::Matchers::Approx(gold[c]));
	// }

	

}

TEST_CASE( "wmwtest_sparse_iter", "[wilcox]" ) {

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
	utils::read_hdf5_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_wilcox_pval_matrix_c.h5", "matrix", goldpv, goldrows, goldcols, rowMajor);


	REQUIRE(goldcols == cols);
	REQUIRE(goldrows == nuniq);
	REQUIRE(rowMajor == false);

	std::vector<double> opv(nuniq * cols, 0);
	std::vector<std::pair<int, size_t> > cl_counts;

	omp_sparse_wmw(x.cbegin(), i.cbegin(), p.cbegin(), rows, cols, labels.cbegin(), 
		2, true, opv, cl_counts, 1);

	// compare
    REQUIRE_THAT(opv, Catch::Matchers::Approx(goldpv));

	opv.clear(); opv.resize(nuniq * cols, 0);
	cl_counts.clear();
	omp_sparse_wmw(x.cbegin(), i.cbegin(), p.cbegin(), rows, cols, labels.cbegin(), 
		2, true, opv, cl_counts, 4);

	// compare
    REQUIRE_THAT(opv, Catch::Matchers::Approx(goldpv));

}



TEST_CASE( "wmwtest_pseudosparse_vec", "[wilcox]" ) {

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

	// count number of unique labels
	std::vector<std::pair<int, size_t> > cl_counts;
	count_clusters_vec(labels, labels.size(), cl_counts, 1);
	size_t nuniq = cl_counts.size();

	// gold
	std::vector<double> goldpv;
	size_t goldrows, goldcols;
	bool rowMajor;
	utils::read_hdf5_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_wilcox_pval_matrix_c.h5", "matrix", goldpv, goldrows, goldcols, rowMajor);

	REQUIRE(goldcols == cols);
	REQUIRE(goldrows == nuniq);
	REQUIRE(rowMajor == false);


	// convert to dense
	std::vector<double> dense(rows * cols, 0);
	csc_to_dense_c_vec(x, i, p, rows, cols, dense, 1);


	std::vector<double> opv(nuniq * cols, 0);
	vecsc_wmw_vecsc(dense, rows, cols, labels, 2, true, opv, cl_counts, 1);

	// compare
    REQUIRE_THAT(opv, Catch::Matchers::Approx(goldpv));


	opv.clear(); opv.resize(nuniq * cols, 0);
	vecsc_wmw_vecsc(dense, rows, cols, labels, 2, true, opv, cl_counts, 4);

	// compare
    REQUIRE_THAT(opv, Catch::Matchers::Approx(goldpv));


	// convert gold:
	// test with vector of vector.
	


	// convert gold:
	// size_t offset = 0;
	// size_t inoffset = 0;
	// std::vector<std::vector<double>> gold(cols);
	// std::vector<std::vector<double>> inmat(cols);
	// std::vector<std::vector<double>> pvm(cols);
	// for (size_t c = 0; c < cols; ++c) {
	// 	inmat[c].resize(rows, 0);
	// 	for (size_t r = 0; r < rows; ++r, ++inoffset) {
	// 		inmat[c][r] = dense[inoffset];
	// 	}

	// 	gold[c].resize(nuniq, 0);
	// 	for (size_t r = 0; r < nuniq; ++r, ++offset) {
	// 		gold[c][r] = goldpv[offset];
	// 	}

	// 	pvm[c].resize(nuniq, 0);
	// }

	// // test with vector of vector.
	// matc_wmw_matc(inmat, rows, cols, labels, 2, true, pvm, cl_counts, 1);

	// for (size_t c = 0; c < cols; ++c) {
	//     REQUIRE_THAT(pvm[c], Catch::Matchers::Approx(gold[c]));
	// }


	// for (size_t c = 0; c < cols; ++c) {
	// 	pvm[c].clear();
	// 	pvm[c].resize(nuniq, 0);
	// }

	// // test with vector of vector.
	// matc_wmw_matc(inmat, rows, cols, labels, 2, true, pvm, cl_counts, 4);

	// for (size_t c = 0; c < cols; ++c) {
	//     REQUIRE_THAT(pvm[c], Catch::Matchers::Approx(gold[c]));
	// }

}



TEST_CASE( "wmwtest_pseudosparse_iter", "[wilcox]" ) {

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
	utils::read_hdf5_matrix(std::string(FDE_TEST_DATA_DIR) + "/test_spmat_wilcox_pval_matrix_c.h5", "matrix", goldpv, goldrows, goldcols, rowMajor);


	REQUIRE(goldcols == cols);
	REQUIRE(goldrows == nuniq);
	REQUIRE(rowMajor == false);


	// convert to dense
	std::vector<double> dense(rows * cols, 0);
	csc_to_dense_c_vec(x, i, p, rows, cols, dense, 1);


	std::vector<double> opv(nuniq * cols, 0);
	std::vector<std::pair<int, size_t> > cl_counts;
	omp_dense_wmw(dense.cbegin(), rows, cols, labels.cbegin(), 2, true, opv, cl_counts, 1);

	// compare
    REQUIRE_THAT(opv, Catch::Matchers::Approx(goldpv));


	opv.clear(); opv.resize(nuniq * cols, 0);
	cl_counts.clear();
	omp_dense_wmw(dense.cbegin(), rows, cols, labels.cbegin(), 2, true, opv, cl_counts, 4);

	// compare
    REQUIRE_THAT(opv, Catch::Matchers::Approx(goldpv));

}
