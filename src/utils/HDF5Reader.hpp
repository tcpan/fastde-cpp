/*
 *  HDFMatrixReader.hpp
 *
 *  Created on: Aug 21, 2020
 *  Author: Tony Pan
 *  Affiliation: Institute for Data Engineering and Science
 *  						Georgia Institute of Technology, Atlanta, GA 30332
 *  
 * serial and parallel hdf5 output.
 * following group structure is used:
 * 	root: no datasets
 *  top level children: contains source expression data as "dataset", and any other data (e.g. gene names, sample names, regulatory factors, etc) 
 * 		labeled with dataset names . e.g. /flower, /leaf, etc.
 *  derived data appears as children immediately below its source data.
 * 		e.g. /flower/mi, /leaf/pearson, /flower/mi/dpi
 *  aggregated results appear under deepest level of its aggregate sources.  naming should be enforced programatically.
 * 		e.g. /flower/mi/dpi/dpi+mi.  
 * 
 * Need tool to ingest/extract datasets
 */

#pragma once

#include <string>  // string
#include <algorithm> 
#include <vector>
#include <string>

#include "utils/benchmark.hpp"
#include "utils/report.hpp"



#ifdef USE_MPI
#include <mpi.h>
#include "utils/mpi_types.hpp"
#endif  // with mpi

#include <hdf5.h>
#include "utils/hdf5_types.hpp"


namespace utils { 


class HDF5Reader {

	protected:
		std::string filename;


		bool readVectorSize(hid_t obj_id, std::string const & path, ssize_t & count) {
			// https://stackoverflow.com/questions/15786626/get-the-dimensions-of-a-hdf5-dataset#:~:text=First%20you%20need%20to%20get,int%20ndims%20%3D%20H5Sget_simple_extent_ndims(dspace)%3B
			auto exists = H5Lexists(obj_id, path.c_str(), H5P_DEFAULT);
			if (exists <= 0) return false;

			hid_t dataset_id = H5Dopen(obj_id, path.c_str(), H5P_DEFAULT);

			hid_t dataspace_id = H5Dget_space(dataset_id);
			const int ndims = H5Sget_simple_extent_ndims(dataspace_id);
			hsize_t dims[ndims];
			H5Sget_simple_extent_dims(dataspace_id, dims, NULL);

			count = dims[0]; 

			H5Sclose(dataspace_id);
			H5Dclose(dataset_id);

			return true;
		}

		template <typename T>
		bool readSparseMatrixSize(hid_t obj_id, std::string const & path, T & rows, T & cols, T& elements) {
			// https://stackoverflow.com/questions/15786626/get-the-dimensions-of-a-hdf5-dataset#:~:text=First%20you%20need%20to%20get,int%20ndims%20%3D%20H5Sget_simple_extent_ndims(dspace)%3B
			auto exists = H5Lexists(obj_id, path.c_str(), H5P_DEFAULT);
			if (exists <= 0) return false;

			// get the group
			hid_t group_id = H5Gopen(obj_id, path.c_str(), H5P_DEFAULT);

			// then get the attribute of the group
			exists = H5Aexists(group_id, "shape");
        	if (exists <= 0)  return false;

           	// open
            hid_t attr_id = H5Aopen(group_id, "shape", H5P_DEFAULT);

			hid_t atype  = H5Aget_type(attr_id);
  			hid_t space_id = H5Aget_space(attr_id);
  			size_t npoints    = H5Sget_simple_extent_npoints(space_id);
			int64_t *dims = (int64_t *)malloc(sizeof(int64_t) * npoints);
			H5Aread(attr_id, atype, dims);

			rows = dims[0]; 
			cols = dims[1];

			free(dims);
			readVectorSize(group_id, "x", elements);

			H5Tclose(atype);
			H5Aclose(attr_id);
			H5Sclose(space_id);
			H5Gclose(group_id);

			return true;
		}

		bool readMatrixSize(hid_t obj_id, std::string const & path, ssize_t & rows, ssize_t & cols) {
			// https://stackoverflow.com/questions/15786626/get-the-dimensions-of-a-hdf5-dataset#:~:text=First%20you%20need%20to%20get,int%20ndims%20%3D%20H5Sget_simple_extent_ndims(dspace)%3B
			auto exists = H5Lexists(obj_id, path.c_str(), H5P_DEFAULT);
			if (exists <= 0) return false;

			hid_t dataset_id = H5Dopen(obj_id, path.c_str(), H5P_DEFAULT);

			hid_t dataspace_id = H5Dget_space(dataset_id);
			const int ndims = H5Sget_simple_extent_ndims(dataspace_id);
			hsize_t dims[ndims];
			H5Sget_simple_extent_dims(dataspace_id, dims, NULL);

			rows = dims[0]; // rows
			cols = dims[1]; // cols.

			H5Sclose(dataspace_id);
			H5Dclose(dataset_id);

			return true;
		}

		// max length is set.
		inline size_t strnlen(char * str, size_t n) {
			size_t i = 0;
			while ((i < n) && (str[i] != 0)) ++i;
			return i;
		}

		bool readStrings(hid_t obj_id, std::string const & path, std::vector<std::string> & out ) {
			// from https://stackoverflow.com/questions/581209/how-to-best-write-out-a-stdvector-stdstring-container-to-a-hdf5-dataset
			// MODIFIED to use C api.

			// open data set
			auto exists = H5Lexists(obj_id, path.c_str(), H5P_DEFAULT);
			if (exists <= 0) return false;

			hid_t dataset_id = H5Dopen(obj_id, path.c_str(), H5P_DEFAULT);

			hid_t filetype_id = H5Dget_type(dataset_id);
			size_t max_len = H5Tget_size(filetype_id);

			// open data space and get dimensions
			hid_t dataspace_id = H5Dget_space(dataset_id);
			const int ndims = H5Sget_simple_extent_ndims(dataspace_id);
			hsize_t dims[ndims];
			H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
			
			// get continuous space now...
			char * data = reinterpret_cast<char *>(calloc( (dims[0] + 1), max_len * sizeof(char)));
			data[dims[0]*max_len] = 0;
			//FMT_ROOT_PRINT("In read STRING dataset, got number of strings: [{}].  temp array at {:p}\n", dims[0], data );

			// prepare output

			// Variable length string type
			// read data.  use mem_type is variable length string.  memspace is same as file space. 
			// NOTE: data should be a continuous memory block.  
			// NOTE: Not sure why online docs show as char** with automatic allocation assumption. API behavior change?
			//auto status = 
			H5Dread(dataset_id, filetype_id, H5S_ALL, dataspace_id, H5P_DEFAULT, data);

			// convert to string objects
			out.clear();
			out.reserve(dims[0]);
			char * ptr = data;
			for(size_t x=0; x < dims[0]; ++x, ptr += max_len)
			{
				// auto l = strlen(ptr);
				// FMT_ROOT_PRINT("GOT STRING {} {:p} {} \"{}\"\n", x, ptr, l, std::string(ptr, l) );
				out.emplace_back(ptr, strnlen(ptr, max_len));
			}
			// H5Dvlen_reclaim (filetype_id, dataspace_id, H5P_DEFAULT, data);
			free(data);
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);
			H5Tclose(filetype_id);
			
			return true;
		}

		// specify number of rows and cols to read.
		template <typename T>
		bool readVector(hid_t obj_id, std::string const & path,
			size_t const & count, T * vector) {

			// open data set
			auto exists = H5Lexists(obj_id, path.c_str(), H5P_DEFAULT);
			if (exists <= 0) return false;

			hid_t dataset_id = H5Dopen(obj_id, path.c_str(), H5P_DEFAULT);

			// open data space and get dimensions
			hid_t filespace_id = H5Dget_space(dataset_id);
			const int ndims = H5Sget_simple_extent_ndims(filespace_id);
			hsize_t file_dims[ndims];
			H5Sget_simple_extent_dims(filespace_id, file_dims, NULL);

			// create target space
			hsize_t mem_dims[ndims] = {static_cast<hsize_t>(count)};
			hid_t memspace_id = H5Screate_simple(ndims, mem_dims, NULL);
			// select hyperslab of memory, for row by row traversal
			hsize_t mstart = 0;  // element offset for first block
			hsize_t mcount = count; // # of blocks
			hsize_t mstride = 1;  // element stride to get to next block
			hsize_t mblock = 1;  // block size  1xcols
			H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, &mstart, &mstride, &mcount, &mblock);

			// float type
			utils::hdf5::datatype<T> type;
			hid_t type_id = type.value;

			// read data.  use mem_type is variable length string.  memspace is same as file space. 
			//auto status = 
			H5Dread(dataset_id, type_id, memspace_id, filespace_id, H5P_DEFAULT, vector);

			H5Sclose(memspace_id);
			H5Sclose(filespace_id);
			H5Dclose(dataset_id);

			return true;
		}


		// specify number of rows and cols to read.
		template <typename T>
		bool readMatrix(hid_t obj_id, std::string const & path,
			size_t const & rows, size_t const & cols, 
			T * vectors, size_t const & stride_bytes) {

			if ((stride_bytes % sizeof(T)) > 0) {
				// unsupported.  This means having to write row by row, and some procs may have to write 0 bytes - complicating PHDF5 write.
            	FMT_PRINT_ERR("ERROR: column stride not a multiple of element data type.  This is not support and will be deprecated.\n");
            	return false;
			}
			
			// open data set
			auto exists = H5Lexists(obj_id, path.c_str(), H5P_DEFAULT);
			if (exists <= 0) return false;

			hid_t dataset_id = H5Dopen(obj_id, path.c_str(), H5P_DEFAULT);
			
			// open data space and get dimensions
			hid_t filespace_id = H5Dget_space(dataset_id);
			const int ndims = H5Sget_simple_extent_ndims(filespace_id);
			hsize_t file_dims[ndims];
			H5Sget_simple_extent_dims(filespace_id, file_dims, NULL);
			
			// create target space
			hsize_t mem_dims[ndims] = {rows, stride_bytes / sizeof(T)};
			hid_t memspace_id = H5Screate_simple(ndims, mem_dims, NULL);
			// select hyperslab of memory, for row by row traversal
			hsize_t mstart[2] = {0, 0};  // element offset for first block
			hsize_t mcount[2] = {rows, 1}; // # of blocks
			hsize_t mstride[2] = {1, stride_bytes / sizeof(T)};  // element stride to get to next block
			hsize_t mblock[2] = {1, cols};  // block size  1xcols
			H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, mstart, mstride, mcount, mblock);

			// float type
			utils::hdf5::datatype<T> type;
			hid_t type_id = type.value;

			// read data.  use mem_type is variable length string.  memspace is same as file space. 
			//auto status = 
			H5Dread(dataset_id, type_id, memspace_id, filespace_id, H5P_DEFAULT, vectors);

			H5Sclose(memspace_id);
			H5Sclose(filespace_id);
			H5Dclose(dataset_id);

			return true;
		}

#ifdef USE_MPI
		template <typename T>
		bool readMatrix(hid_t obj_id, std::string const & path,
			size_t const & rows, size_t const & cols, 
			T * vectors, size_t const & stride_bytes,
			MPI_Comm const & comm) {

			if ((stride_bytes % sizeof(T)) > 0) {
				// unsupported.  This means having to write row by row, and some procs may have to write 0 bytes - complicating PHDF5 write.
            	FMT_PRINT_ERR("ERROR: column stride not a multiple of element data type.  This is not support and will be deprecated.\n");
            	return false;
			}

			int procs, rank;
			MPI_Comm_size(comm, &procs);
			MPI_Comm_rank(comm, &rank);
			
			size_t row_offset = rows;
			MPI_Exscan(MPI_IN_PLACE, &row_offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, comm);
			if (rank == 0) row_offset = 0;

			// open data set
			auto exists = H5Lexists(obj_id, path.c_str(), H5P_DEFAULT);
			if (exists <= 0) return false;

			hid_t dataset_id = H5Dopen(obj_id, path.c_str(), H5P_DEFAULT);

			// open data space and get dimensions
			hid_t filespace_id = H5Dget_space(dataset_id);
			const int ndims = H5Sget_simple_extent_ndims(filespace_id);
			hsize_t file_dims[ndims];
			H5Sget_simple_extent_dims(filespace_id, file_dims, NULL);
			
			// each process defines dataset in memory and hyperslab in file.
			hsize_t start[2] = {row_offset, 0};  // starting offset, row, then col.
			hsize_t count[2] = {rows, cols};   // number of row and col blocks.
			H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, start, NULL, count, NULL);

			hsize_t mem_dims[ndims] = {rows, stride_bytes / sizeof(T) };
			hid_t memspace_id = H5Screate_simple(ndims, mem_dims, NULL);
			// select hyperslab of memory, for row by row traversal
			hsize_t mstart[2] = {0, 0};  // element offset for first block
			hsize_t mcount[2] = {rows, 1}; // # of blocks
			hsize_t mstride[2] = {1, stride_bytes / sizeof(T)};  // element stride to get to next block
			hsize_t mblock[2] = {1, cols};  // block size  1xcols
			H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, mstart, mstride, mcount, mblock);

			// float type
			utils::hdf5::datatype<T> type;
			hid_t type_id = type.value;

			// read data.  use mem_type is variable length string.  memspace is same as file space. 
			hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
			H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
			//auto status = 
			H5Dread(dataset_id, type_id, memspace_id, filespace_id, plist_id, vectors);

			// convert to string objects
			H5Sclose(memspace_id);
			H5Sclose(filespace_id);
			H5Dclose(dataset_id);

			H5Pclose(plist_id);

			return true;
		}
#endif


	public:

		HDF5Reader(std::string const & _filename) : filename(_filename) {}
		virtual ~HDF5Reader() {}

			// // open the file for reading only.
			// fid = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);


		/*get gene expression matrix size*/
		bool getVectorSize(std::string const & path, ssize_t& count) {
			auto stime = getSysTime();

			// open the file for reading only.
			hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
			if (file_id < 0) return false;
			
			hid_t group_id;
            auto status = H5Lexists(file_id, path.c_str(), H5P_DEFAULT);
            if (status > 0) {
                group_id = H5Gopen(file_id, path.c_str(), H5P_DEFAULT);
            } else {
                FMT_PRINT_ERR("WARN: unable to get group {} in file {}\n", path, filename);
				H5Fclose(file_id);
                return false;
            }

			bool res = readVectorSize(group_id, "block0_values", count);

			H5Gclose(group_id);
			H5Fclose(file_id);

			auto etime = getSysTime();
			FMT_ROOT_PRINT("get vector size {} in {} sec\n", count, get_duration_s(stime, etime));
			return res;
		}
		bool getSparseMatrixSize(std::string const & path, ssize_t& rows, ssize_t& cols, ssize_t& nelements) {
						auto stime = getSysTime();

			// open the file for reading only.
			hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
			if (file_id < 0) return false;
			
			bool res = readSparseMatrixSize(file_id, path.c_str(), rows, cols, nelements);


			H5Fclose(file_id);

			auto etime = getSysTime();
			FMT_ROOT_PRINT("get sparse matrix size {}x{} in {} sec\n", rows, cols, get_duration_s(stime, etime));
			return res;

		}
		bool getMatrixSize(std::string const & path, ssize_t& rows, ssize_t& cols) {
			auto stime = getSysTime();

			// open the file for reading only.
			hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
			if (file_id < 0) return false;
			
			hid_t group_id;
            auto status = H5Lexists(file_id, path.c_str(), H5P_DEFAULT);
            if (status > 0) {
                group_id = H5Gopen(file_id, path.c_str(), H5P_DEFAULT);
            } else {
                FMT_PRINT_ERR("WARN: unable to get group {} in file {}\n", path, filename);
				H5Fclose(file_id);
                return false;
            }

			bool res = readMatrixSize(group_id, "block0_values", rows, cols);


			H5Gclose(group_id);
			H5Fclose(file_id);

			auto etime = getSysTime();
			FMT_ROOT_PRINT("get matrix size {}x{} in {} sec\n", rows, cols, get_duration_s(stime, etime));
			return res;
		}

		template <typename T>
		bool loadVectorData(std::string const & path, 
			ssize_t const & count, 
			T* vectors) {

			// open the file for reading only.
			hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
			if (file_id < 0) return false;

			hid_t group_id;
            auto status = H5Lexists(file_id, path.c_str(), H5P_DEFAULT);
            if (status > 0) {
                group_id = H5Gopen(file_id, path.c_str(), H5P_DEFAULT);
            } else {
                FMT_PRINT_ERR("WARN: unable to get group {} in file {}\n", path, filename);
				H5Fclose(file_id);
                return false;
            }

			// read the data.
			readVector(group_id, "block0_values", count, vectors);

			H5Gclose(group_id);
			H5Fclose(file_id);
			return true;
		}

		template <typename T, typename I, typename P>
		bool loadSparseMatrixData(std::string const & path,
			const ssize_t rows, const ssize_t cols, const ssize_t nz,
			T* x, I* i, P* p) {

			// open the file for reading only.
			hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
			if (file_id < 0) return false;

			hid_t group_id;
            auto status = H5Lexists(file_id, path.c_str(), H5P_DEFAULT);
            if (status > 0) {
                group_id = H5Gopen(file_id, path.c_str(), H5P_DEFAULT);
            } else {
                FMT_PRINT_ERR("WARN: unable to get group {} in file {}\n", path, filename);
				H5Fclose(file_id);
                return false;
            }

			// read the data.

			readVector(group_id, "x", nz, x);
			readVector(group_id, "i", nz, i);
			readVector(group_id, "p", cols + 1, p);

			H5Gclose(group_id);
			H5Fclose(file_id);
			return true;
		}


		template <typename T>
		bool loadMatrixData(std::string const & path, ssize_t const & rows,
			ssize_t const & cols, 
			T* vectors, ssize_t const & stride_bytes) {

			// open the file for reading only.
			hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
			if (file_id < 0) return false;

			hid_t group_id;
            auto status = H5Lexists(file_id, path.c_str(), H5P_DEFAULT);
            if (status > 0) {
                group_id = H5Gopen(file_id, path.c_str(), H5P_DEFAULT);
            } else {
                FMT_PRINT_ERR("WARN: unable to get group {} in file {}\n", path, filename);
				H5Fclose(file_id);
                return false;
            }

			// // read the names.
			// readStrings(group_id, "axis1", genes);
			// readStrings(group_id, "axis0", samples);

			// read the data.
			readMatrix(group_id, "block0_values", rows, cols, vectors, stride_bytes);

			H5Gclose(group_id);
			H5Fclose(file_id);
			return true;
		}


};
}