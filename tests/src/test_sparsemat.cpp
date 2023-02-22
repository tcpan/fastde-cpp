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



