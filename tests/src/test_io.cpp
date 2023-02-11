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
		enum test_type : int { vector = 0, matrix = 1, sparse_matrix = 2 };

		test_type test_method;

		app_parameters() : test_method(test_type::vector) {}
		virtual ~app_parameters() {}

		virtual void config(CLI::App& app) {
            app.add_option("-m,--method", test_method, "reader impl: vector=0, matrix=1, sparse_matrix=2");
		}
		virtual void print(const char * prefix) {
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


	auto stime = getSysTime();
	auto etime = getSysTime();

	size_t rows = common_params.num_vectors;
	size_t cols = common_params.vector_size;
	bool dummy = true; 
	if (app_params.test_method == app_parameters::test_type::matrix) {

		std::vector<double> x;

		stime = getSysTime();
		if (common_params.random) {
			x = utils::make_random_matrix<double, int>(common_params.rseed, 
				common_params.rmin, common_params.rmax,
				common_params.num_vectors, common_params.vector_size);
		}
		etime = getSysTime();
		FMT_ROOT_PRINT_RT("[TIME] Load data in {} sec\n", get_duration_s(stime, etime));
		// input.print("INPUT: ");

		utils::write_hdf5_matrix("temp_matrix.h5", "matrix", x, rows, cols, dummy);

		std::vector<double> x2;
		size_t rows2 = 0, cols2 = 0;
		utils::read_hdf5_matrix("temp_matrix.h5", "matrix", x2, rows2, cols2, dummy);

		// now compare.
		if ((rows != rows2) || (cols != cols2)) FMT_PRINT_ERR("mat io test dimension mismatch");
		if ((x == x2) == false) FMT_PRINT_ERR("mat io test x mismatch");

 	} if (app_params.test_method == app_parameters::test_type::vector) {

		std::vector<double> x;

		stime = getSysTime();
		if (common_params.random) {
			x = utils::make_random_vector<double>(common_params.rseed, 
				common_params.rmin, common_params.rmax,
				common_params.vector_size);
		}
		etime = getSysTime();
		FMT_ROOT_PRINT_RT("[TIME] Load data in {} sec\n", get_duration_s(stime, etime));
		// input.print("INPUT: ");
		size_t count = common_params.vector_size;

		utils::write_hdf5_vector("temp_double_vec.h5", "vector", x, count);

		std::vector<double> x2;
		size_t count2 = 0;
		utils::read_hdf5_vector("temp_double_vec.h5", "vector", x2, count2);

		// compare
		if (count != count2) FMT_PRINT_ERR("vec io test double dimension mismatch");
		if ((x == x2) == false) FMT_PRINT_ERR("vec io test double mismatch");


		std::vector<int> y;

		stime = getSysTime();
		if (common_params.random) {
			y = utils::make_random_vector<int>(common_params.rseed, 
				common_params.rmin, common_params.rmax,
				common_params.vector_size);
		}
		etime = getSysTime();
		FMT_ROOT_PRINT_RT("[TIME] Load data in {} sec\n", get_duration_s(stime, etime));
		// input.print("INPUT: ");

		utils::write_hdf5_vector("temp_int_vec.h5", "vector", y, count);

		std::vector<int> y2;
		count2 = 0;
		utils::read_hdf5_vector("temp_int_vec.h5", "vector", y2, count2);

		if (count != count2) FMT_PRINT_ERR("vec io test int dimension mismatch");
		if ((y == y2) == false) FMT_PRINT_ERR("vec io test int mismatch");

 	} if (app_params.test_method == app_parameters::test_type::sparse_matrix) {

		std::vector<double> x;
		std::vector<int> i;
		std::vector<int> p;

		stime = getSysTime();
		if (common_params.random) {
			std::tie(x, i, p) = utils::make_random_sparse_matrix<double, int>(common_params.rseed, 
				common_params.rmin, common_params.rmax,
				common_params.num_vectors, common_params.vector_size,
				common_params.sparsity);
		}
		etime = getSysTime();
		FMT_ROOT_PRINT_RT("[TIME] Load data with elements {} in {} sec\n", x.size(), get_duration_s(stime, etime));
		// input.print("INPUT: ");

		utils::write_hdf5_sparse_matrix("temp_spmat.h5", "sparse_matrix", x, i, p, rows, cols, dummy);

		std::vector<double> x2;
		std::vector<int> i2;
		std::vector<int> p2;
		size_t rows2 = 0, cols2 = 0;

		utils::read_hdf5_sparse_matrix("temp_spmat.h5", "sparse_matrix", x2, i2, p2, rows2, cols2, dummy);

		if ((rows != rows2) || (cols != cols2)) FMT_PRINT_ERR("spmat io test dimension mismatch");
		if ((x == x2) == false) FMT_PRINT_ERR("spmat io test x mismatch");
		if ((i == i2) == false) FMT_PRINT_ERR("spmat io test i mismatch");
		if ((p == p2) == false) FMT_PRINT_ERR("spmat io test p mismatch");

	}

	return 0;
}
