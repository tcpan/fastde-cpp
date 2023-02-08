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

class app_parameters : public parameters_base {
	public:
		enum atof_type : int { std = 0, fast = 1 };
		enum test_type : int { transpose = 0, to_dense = 1, transpose_to_dense = 2, cbind = 3, rbind = 4, row_sum = 5, col_sum = 6, all = 99 };

		atof_type atof_method;
		test_type test_method;

		app_parameters() : atof_method(fast), test_method(test_type::all) {}
		virtual ~app_parameters() {}

		virtual void config(CLI::App& app) {
            app.add_option("-a,--atof", atof_method, "atof impl: std=0, fast=1");
            app.add_option("-m,--method", test_method, "reader impl: lightpcc=0, delim_str=1, delim_char=2");
		}
		virtual void print(const char * prefix) {
            FMT_ROOT_PRINT("{} atof method: {}\n", prefix, (atof_method == std ? "std::atof" : 
				"fast"));
            FMT_ROOT_PRINT("{} reader method: {}\n", prefix, test_method);
		}
};

int main(int argc, char* argv[]) {

	//==============  PARSE INPUT =====================
	CLI::App app{"EXP Reader Test"};

	// handle MPI (TODO: replace with MXX later)
	utils::mpi_parameters mpi_params(argc, argv);
	mpi_params.config(app);

	// set up CLI parsers.
	utils::common_parameters common_params;
	app_parameters app_params;

	common_params.config(app);
	app_params.config(app);

	// parse
	CLI11_PARSE(app, argc, argv);

	// print out, for fun.
	FMT_ROOT_PRINT_RT("command line: ");
	for (int i = 0; i < argc; ++i) {
		FMT_ROOT_PRINT("{} ", argv[i]);
	}
	FMT_ROOT_PRINT("\n");


#ifdef USE_OPENMP
	// omp_set_dynamic(0);
	omp_set_num_threads(common_params.num_threads);
	FMT_PRINT_RT("omp num threads %d.  user threads %lu\n", omp_get_max_threads(), common_params.num_threads);
#endif

	// =============== SETUP INPUT ===================
	// NOTE: input data is replicated on all MPI procs.
	std::vector<double> x;
	std::vector<int> i;
	std::vector<int> p;

	auto stime = getSysTime();
	auto etime = getSysTime();

	size_t rows = common_params.num_vectors;
	size_t cols = common_params.vector_size;


	stime = getSysTime();
	if (common_params.random) {
		std::tie(x, i, p) = utils::make_random_sparse_matrix<double, int>(common_params.rseed, 
			common_params.rmin, common_params.rmax,
			common_params.num_vectors, common_params.vector_size,
			common_params.sparsity);
	}
	etime = getSysTime();
	FMT_ROOT_PRINT("Load data in {} sec\n", get_duration_s(stime, etime));
	// input.print("INPUT: ");

	utils::write_hdf5_sparse_matrix(common_params.output + ".colsum_iter.t4.h5", "vector", x, i, p, rows, cols, true);


	if (mpi_params.rank == 0) {
		mpi_params.print("[PARAM] ");
		common_params.print("[PARAM] ");
		app_params.print("[PARAM] ");
	}

	// ===== DEBUG ====== WRITE OUT INPUT =========
		// write to file.  MPI enabled.  Not thread enabled.
	if ((app_params.test_method == app_parameters::test_type::transpose) || 
		(app_params.test_method == app_parameters::test_type::all)) {
		std::vector<double> ox(x.size(), 0);
		std::vector<int> oi(i.size(), 0);
		std::vector<int> op(rows+1, 0);

		stime = getSysTime();
		csc_transpose_csc(x.cbegin(), i.cbegin(), p.cbegin(), rows, cols, ox.begin(), oi.begin(), op.begin(), 1);
		etime = getSysTime();
		FMT_ROOT_PRINT("transposed in {} sec\n", get_duration_s(stime, etime));

		utils::write_hdf5_sparse_matrix(common_params.output + ".transpose.t1.h5", "sparse_matrix", ox, oi, op, rows, cols, true);

		stime = getSysTime();
		csc_transpose_csc(x.cbegin(), i.cbegin(), p.cbegin(),  rows, cols, ox.begin(), oi.begin(), op.begin(), 4);
		etime = getSysTime();
		FMT_ROOT_PRINT("transposed in {} sec, 4 threads\n", get_duration_s(stime, etime));

		utils::write_hdf5_sparse_matrix(common_params.output + ".transpose.t4.h5", "sparse_matrix", ox, oi, op, rows, cols, true);



		stime = getSysTime();
		csc_transpose_csc_vec(x, i, p, rows, cols, ox, oi, op, 1);
		etime = getSysTime();
		FMT_ROOT_PRINT("transposed in {} sec\n", get_duration_s(stime, etime));

		utils::write_hdf5_sparse_matrix(common_params.output + ".transpose_vec.t1.h5", "sparse_matrix", ox, oi, op, rows, cols, true);

		stime = getSysTime();
		csc_transpose_csc_vec(x, i, p, rows, cols, ox, oi, op, 4);
		etime = getSysTime();
		FMT_ROOT_PRINT("transposed in {} sec, 4 threads\n", get_duration_s(stime, etime));

		utils::write_hdf5_sparse_matrix(common_params.output + ".transpose_vec.t4.h5", "sparse_matrix", ox, oi, op, rows, cols, true);

	}

	if ((app_params.test_method == app_parameters::test_type::to_dense) || 
		(app_params.test_method == app_parameters::test_type::all)) {

		std::vector<double> out(rows * cols, 0);

		stime = getSysTime();
		csc_to_dense_c(x.cbegin(), i.cbegin(), p.cbegin(),  rows, cols, out.begin(), 1);
		etime = getSysTime();
		FMT_ROOT_PRINT("to dense in {} sec\n", get_duration_s(stime, etime));

		utils::write_hdf5_matrix(common_params.output + ".to_dense.t1.h5", "matrix", out, rows, cols, true);
		std::fill(out.begin(), out.end(), 0);

		stime = getSysTime();
		csc_to_dense_c(x.cbegin(), i.cbegin(), p.cbegin(),  rows, cols, out.begin(), 4);
		etime = getSysTime();
		FMT_ROOT_PRINT("to dense in {} sec, 4 threads\n", get_duration_s(stime, etime));

		utils::write_hdf5_matrix(common_params.output + ".to_dense.t4.h5", "matrix", out, rows, cols, true);


		
		stime = getSysTime();
		csc_to_dense_c_vec(x, i, p,  rows, cols, out, 1);
		etime = getSysTime();
		FMT_ROOT_PRINT("to dense in {} sec\n", get_duration_s(stime, etime));

		utils::write_hdf5_matrix(common_params.output + ".to_dense_vec.t1.h5", "matrix", out, rows, cols, true);
		std::fill(out.begin(), out.end(), 0);

		stime = getSysTime();
		csc_to_dense_c_vec(x, i, p,  rows, cols, out, 4);
		etime = getSysTime();
		FMT_ROOT_PRINT("to dense in {} sec, 4 threads\n", get_duration_s(stime, etime));

		utils::write_hdf5_matrix(common_params.output + ".to_dense_vec.t4.h5", "matrix", out, rows, cols, true);
	}
	if ((app_params.test_method == app_parameters::test_type::transpose_to_dense) || 
		(app_params.test_method == app_parameters::test_type::all)) {
		std::vector<double> out(rows * cols);

		stime = getSysTime();
		csc_to_dense_transposed_c(x.cbegin(), i.cbegin(), p.cbegin(), rows, cols, out.begin(), 1);
		etime = getSysTime();
		FMT_ROOT_PRINT("to dense transposed in {} sec\n", get_duration_s(stime, etime));

		utils::write_hdf5_matrix(common_params.output + ".transpose_to_dense.t1.h5", "matrix", out, cols, rows, true);
		std::fill(out.begin(), out.end(), 0);

		stime = getSysTime();
		csc_to_dense_transposed_c(x.cbegin(), i.cbegin(), p.cbegin(), rows, cols, out.begin(), 4);
		etime = getSysTime();
		FMT_ROOT_PRINT("to dense transposed in {} sec, 4 threads\n", get_duration_s(stime, etime));

		utils::write_hdf5_matrix(common_params.output + ".transpose_to_dense.t4.h5", "matrix", out, cols, rows, true);

		
		stime = getSysTime();
		csc_to_dense_transposed_c_vec(x, i, p,  rows, cols, out, 1);
		etime = getSysTime();
		FMT_ROOT_PRINT("to dense transposed in {} sec\n", get_duration_s(stime, etime));

		utils::write_hdf5_matrix(common_params.output + ".transpose_to_dense_vec.t1.h5", "matrix", out, cols, rows, true);
		std::fill(out.begin(), out.end(), 0);

		stime = getSysTime();
		csc_to_dense_transposed_c_vec(x, i, p,  rows, cols, out, 4);
		etime = getSysTime();
		FMT_ROOT_PRINT("to dense transposed in {} sec, 4 threads\n", get_duration_s(stime, etime));

		utils::write_hdf5_matrix(common_params.output + ".transpose_to_dense_vec.t4.h5", "matrix", out, cols, rows, true);

	}
	if ((app_params.test_method == app_parameters::test_type::cbind) || 
		(app_params.test_method == app_parameters::test_type::all)) {

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

		stime = getSysTime();
		csc_cbind_vec(ix, ii, ip, vrows, vcols, ox, oi, op, 1);
		etime = getSysTime();
		FMT_ROOT_PRINT("cbind in {} sec\n", get_duration_s(stime, etime));

		utils::write_hdf5_sparse_matrix(common_params.output + ".cbind_vec.t1.h5", "sparse_matrix", ox, oi, op, rows, cols * repeats, true);

		stime = getSysTime();
		csc_cbind_vec(ix, ii, ip, vrows, vcols, ox, oi, op, 4);
		etime = getSysTime();
		FMT_ROOT_PRINT("cbind in {} sec, 4 threads\n", get_duration_s(stime, etime));

		utils::write_hdf5_sparse_matrix(common_params.output + ".cbind_vec.t4.h5", "sparse_matrix", ox, oi, op, rows, cols * repeats, true);


		std::vector<std::vector<double>::const_iterator> ixp(repeats);
		std::vector<std::vector<int>::const_iterator> iip(repeats);
		std::vector<std::vector<int>::const_iterator> ipp(repeats);
		for (int v = 0; v < repeats; ++v) {
			ixp[v] = x.cbegin();
			iip[v] = i.cbegin();
			ipp[v] = p.cbegin();
		}



		stime = getSysTime();
		csc_cbind(ixp, iip, ipp, vrows, vcols, ox.begin(), oi.begin(), op.begin(), 1);
		etime = getSysTime();
		FMT_ROOT_PRINT("cbind iter in {} sec\n", get_duration_s(stime, etime));

		utils::write_hdf5_sparse_matrix(common_params.output + ".cbind_iter.t1.h5", "sparse_matrix", ox, oi, op, rows, cols * repeats, true);

		stime = getSysTime();
		csc_cbind(ixp, iip, ipp,  vrows, vcols, ox.begin(), oi.begin(), op.begin(), 4);
		etime = getSysTime();
		FMT_ROOT_PRINT("cbind iter in {} sec, 4 threads\n", get_duration_s(stime, etime));

		utils::write_hdf5_sparse_matrix(common_params.output + ".cbind_iter.t4.h5", "sparse_matrix", ox, oi, op, rows, cols * repeats, true);


	}
	if ((app_params.test_method == app_parameters::test_type::rbind) || 
		(app_params.test_method == app_parameters::test_type::all)) {

		int repeats = 3;

		std::vector<std::vector<double>> ix(repeats);
		std::vector<std::vector<int>> ii(repeats);
		std::vector<std::vector<int>> ip(repeats);
		
		for (int v = 0; v < repeats; ++v) {
			ix[v] = x;
			ii[v] = i;
			ip[v] = p;
		}

		std::vector<int> vrows(repeats, rows);
		std::vector<int> vcols(repeats, cols);

		std::vector<double> ox(x.size() * repeats);
		std::vector<int> oi(i.size() * repeats);
		std::vector<int> op(cols + 1, 0);

		stime = getSysTime();
		csc_rbind_vec(ix, ii, ip, vrows, vcols, ox, oi, op, 1);
		etime = getSysTime();
		FMT_ROOT_PRINT("rbind in {} sec\n", get_duration_s(stime, etime));

		utils::write_hdf5_sparse_matrix(common_params.output + ".rbind_vec.t1.h5", "sparse_matrix", ox, oi, op, rows * repeats, cols, true);

		stime = getSysTime();
		csc_rbind_vec(ix, ii, ip, vrows, vcols, ox, oi, op, 4);
		etime = getSysTime();
		FMT_ROOT_PRINT("rbind in {} sec, 4 threads\n", get_duration_s(stime, etime));

		utils::write_hdf5_sparse_matrix(common_params.output + ".rbind_vec.t4.h5", "sparse_matrix", ox, oi, op, rows * repeats, cols, true);


		std::vector<std::vector<double>::const_iterator> ixp(repeats);
		std::vector<std::vector<int>::const_iterator> iip(repeats);
		std::vector<std::vector<int>::const_iterator> ipp(repeats);
		for (int v = 0; v < repeats; ++v) {
			ixp[v] = x.cbegin();
			iip[v] = i.cbegin();
			ipp[v] = p.cbegin();
		}

		stime = getSysTime();
		csc_rbind(ixp, iip, ipp,  vrows, vcols, ox.begin(), oi.begin(), op.begin(), 1);
		etime = getSysTime();
		FMT_ROOT_PRINT("rbind iter in {} sec\n", get_duration_s(stime, etime));

		utils::write_hdf5_sparse_matrix(common_params.output + ".rbind_iter.t1.h5", "sparse_matrix", ox, oi, op, rows * repeats, cols, true);

		stime = getSysTime();
		csc_rbind(ixp, iip, ipp,  vrows, vcols, ox.begin(), oi.begin(), op.begin(), 4);
		etime = getSysTime();
		FMT_ROOT_PRINT("rbind iter in {} sec, 4 threads\n", get_duration_s(stime, etime));

		utils::write_hdf5_sparse_matrix(common_params.output + ".rbind_iter.t4.h5", "sparse_matrix", ox, oi, op, rows * repeats, cols, true);

	}
	if ((app_params.test_method == app_parameters::test_type::col_sum) || 
		(app_params.test_method == app_parameters::test_type::all)) {


		std::vector<double> out(cols, 0);

		stime = getSysTime();
		csc_colsums_vec(x, p, cols, out, 1);
		etime = getSysTime();
		FMT_ROOT_PRINT("colsums vec in {} sec\n", get_duration_s(stime, etime));

		utils::write_hdf5_vector(common_params.output + ".colsum_vec.t1.h5", "vector", out, cols);

		stime = getSysTime();
		csc_colsums_vec(x, p, cols, out, 4);
		etime = getSysTime();
		FMT_ROOT_PRINT("colsums vec in {} sec, 4 threads\n", get_duration_s(stime, etime));

		utils::write_hdf5_vector(common_params.output + ".colsum_vec.t4.h5", "vector", out, cols);

		stime = getSysTime();
		csc_colsums_iter(x.cbegin(), p.cbegin(), cols, out.begin(), 1);
		etime = getSysTime();
		FMT_ROOT_PRINT("colsums iter in {} sec\n", get_duration_s(stime, etime));

		utils::write_hdf5_vector(common_params.output + ".colsum_iter.t1.h5", "vector", out, cols);

		stime = getSysTime();
		csc_colsums_iter(x.cbegin(), p.cbegin(), cols, out.begin(), 4);
		etime = getSysTime();
		FMT_ROOT_PRINT("colsums iter in {} sec, 4 threads\n", get_duration_s(stime, etime));

		utils::write_hdf5_vector(common_params.output + ".colsum_iter.t4.h5", "vector", out, cols);

	}
	if ((app_params.test_method == app_parameters::test_type::row_sum) || 
		(app_params.test_method == app_parameters::test_type::all)) {

		std::vector<double> out(rows, 0);
		size_t nz = x.size();

		stime = getSysTime();
		csc_rowsums_vec(x, i, rows, out, 1);
		etime = getSysTime();
		FMT_ROOT_PRINT("rowsums vec in {} sec\n", get_duration_s(stime, etime));
		utils::write_hdf5_vector(common_params.output + ".rowsum_vec.t1.h5", "vector", out, rows);

		stime = getSysTime();
		csc_rowsums_vec(x, i, rows, out, 4);
		etime = getSysTime();
		FMT_ROOT_PRINT("rowsums vec in {} sec, 4 threads\n", get_duration_s(stime, etime));
		utils::write_hdf5_vector(common_params.output + ".rowsum_vec.t4.h5", "vector", out, rows);

		stime = getSysTime();
		csc_rowsums_iter(x.cbegin(), i.cbegin(), rows, nz, out.begin(), 1);
		etime = getSysTime();
		FMT_ROOT_PRINT("rowsums iter in {} sec\n", get_duration_s(stime, etime));
		utils::write_hdf5_vector(common_params.output + ".rowsum_iter.t1.h5", "vector", out, rows);

		stime = getSysTime();
		csc_rowsums_iter(x.cbegin(), i.cbegin(), rows, nz, out.begin(), 4);
		etime = getSysTime();
		FMT_ROOT_PRINT("rowsums iter in {} sec, 4 threads\n", get_duration_s(stime, etime));
		utils::write_hdf5_vector(common_params.output + ".rowsum_iter.t4.h5", "vector", out, rows);
	}

	return 0;
}
