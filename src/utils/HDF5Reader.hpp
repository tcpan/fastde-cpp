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



#include <hdf5.h>
#include "utils/hdf5_types.hpp"


// HDF5 stores data in ROW MAJOR order:
// HDF5 uses C storage conventions, assuming that the last listed dimension is the fastest-changing dimension and the first-listed dimension is the slowest changing.
// from the HDF5 User's Guide.
// ACTUALLY, the output is exactly matching to memory, either col major or row major.  so we need some annotation of "C" and "F"

namespace utils { 


class HDF5Reader {

	protected:
		std::string filename;

		template <typename T>
		bool readVectorSize(hid_t obj_id, std::string const & path, T & count) {
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

		template <typename T, typename S>
		bool readSparseMatrixSize(hid_t obj_id, std::string const & path, T & rows, T & cols, S& elements) {
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
  			size_t npoints  = H5Sget_simple_extent_npoints(space_id);
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

		// template <typename T>
		// bool readMatrixSize(hid_t obj_id, std::string const & path, T & rows, T & cols) {
		// 	// https://stackoverflow.com/questions/15786626/get-the-dimensions-of-a-hdf5-dataset#:~:text=First%20you%20need%20to%20get,int%20ndims%20%3D%20H5Sget_simple_extent_ndims(dspace)%3B
		// 	auto exists = H5Lexists(obj_id, path.c_str(), H5P_DEFAULT);
		// 	if (exists <= 0) return false;

		// 	hid_t dataset_id = H5Dopen(obj_id, path.c_str(), H5P_DEFAULT);

		// 	hid_t dataspace_id = H5Dget_space(dataset_id);
		// 	const int ndims = H5Sget_simple_extent_ndims(dataspace_id);
		// 	hsize_t dims[ndims];
		// 	H5Sget_simple_extent_dims(dataspace_id, dims, NULL);

		// 	rows = dims[0]; // rows
		// 	cols = dims[1]; // cols.

		// 	H5Sclose(dataspace_id);
		// 	H5Dclose(dataset_id);

		// 	return true;
		// }

		template <typename T>
		bool readMatrixSize(hid_t obj_id, std::string const & path, T & rows, T & cols) {
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
  			size_t npoints  = H5Sget_simple_extent_npoints(space_id);
			int64_t *dims = (int64_t *)malloc(sizeof(int64_t) * npoints);
			H5Aread(attr_id, atype, dims);

			rows = dims[0]; 
			cols = dims[1];

			free(dims);

			H5Tclose(atype);
			H5Aclose(attr_id);
			H5Sclose(space_id);
			H5Gclose(group_id);

			return true;
		}


		// max length is set.
		inline size_t strnlen(char * str, size_t n) {
			size_t i = 0;
			while ((i < n) && (str[i] != 0)) ++i;
			return i;
		}

		//https://support.hdfgroup.org/ftp/HDF5/examples/misc-examples/attrvstr.c
		bool readStringAttribute(hid_t obj_id, std::string const & name, std::string & value) {
			auto exists = H5Aexists(obj_id, name.c_str());
			if (exists <= 0) return false;

			// open
			hid_t attr_id = H5Aopen(obj_id, name.c_str(), H5P_DEFAULT);
			hid_t type_id = H5Aget_type(attr_id);
		

			// string could be fixed length or variable. 
			//  variable length would have H5Tget_size(type) == sizeof(char *)
			//  fixed length would have H5Tget_size(type) == actual size.
			//  H5Aget_storage_size would also return actual size.
			// not sure what to do when the sizes all equal....


			hsize_t ptrsize = H5Tget_size(type_id);

			hsize_t sz = H5Aget_storage_size(attr_id);
			// FMT_PRINT("size found {} {} \n", ptrsize, sz);

			if (ptrsize == sz) {
				// fixed length (unless sz is same as sizeof(char *)) in which case I am not sure....
				char * val = new char[sz + 1];  val[sz] = 0;
				H5Aread(attr_id, type_id, val);
				value = std::string(val);
				delete [] val;
			} else {
				// variable length.  ptrsize is the size of sizeof(char *)
				char * val =0;
				H5Aread(attr_id, type_id, &val);
				value = std::string(val);
				H5free_memory(val);
			}


			// H5T_class_t type_class = H5Tget_class (type_id);   

			// if (type_class == H5T_STRING) printf ("File datatype has class H5T_STRING\n");
			// size_t size = H5Tget_size(type_id);
		    // printf(" Size is of the file datatype returned by H5Tget_size %d \n This is a size of char pointer\n Use H5Tis_variable_str call instead \n", size);

			// htri_t size_var;
			// if((size_var = H5Tis_variable_str(type_id)) == 1)  printf(" to find if string has variable size \n");

		    // hid_t type = H5Tget_native_type(type_id, H5T_DIR_ASCEND);

			H5Tclose(type_id);
			H5Aclose(attr_id);
			// H5Tclose(type);

			return true;
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
		template <typename T, typename S>
		bool readVector(hid_t obj_id, std::string const & path,
			S const & count, T * vector) {

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
			hsize_t mstart[ndims] = {0};  // element offset for first block
			hsize_t mcount[ndims] = {1}; // # of blocks
			hsize_t mstride[ndims] = {static_cast<hsize_t>(count)};  // element stride to get to next block
			hsize_t mblock[ndims] = {static_cast<hsize_t>(count)};  // block size  1xcols
			H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, mstart, mstride, mcount, mblock);

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
		template <typename T, typename S>
		bool readMatrix(hid_t obj_id, std::string const & path,
			S const & rows, S const & cols, bool const & rowMajor, 
			T * vectors) {
			
			// FMT_PRINT("READING MATRIX in {} major\n", (rowMajor ? "row" : "column"));

			// open data set
			auto exists = H5Lexists(obj_id, path.c_str(), H5P_DEFAULT);
			if (exists <= 0) return false;

			hid_t dataset_id = H5Dopen(obj_id, path.c_str(), H5P_DEFAULT);
			
			// open data space and get FILE dimensions.  this is contiguous.
			hid_t filespace_id = H5Dget_space(dataset_id);
			const int ndims = H5Sget_simple_extent_ndims(filespace_id);
			hsize_t file_dims[ndims];
			H5Sget_simple_extent_dims(filespace_id, file_dims, NULL);
			
			// create target space.  determined by the rowMajor flag.
			// see http://davis.lbl.gov/Manuals/HDF5-1.8.7/UG/12_Dataspaces.html
			hid_t memspace_id;
			if (rowMajor) {
				hsize_t mem_dims[ndims] = {static_cast<hsize_t>(rows), static_cast<hsize_t>(cols)};
				memspace_id = H5Screate_simple(ndims, mem_dims, NULL);
				// select hyperslab of memory, for row by row traversal
				hsize_t mstart[2] = {0, 0};  // element offset for first block
				hsize_t mcount[2] = {static_cast<hsize_t>(rows), 1}; // # of blocks
				hsize_t mstride[2] = {1, static_cast<hsize_t>(cols)};  // element stride to get to next block
				hsize_t mblock[2] = {1, static_cast<hsize_t>(cols)};  // block size  1xcols
				H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, mstart, mstride, mcount, mblock);
			} else {
				hsize_t mem_dims[ndims] = {static_cast<hsize_t>(cols), static_cast<hsize_t>(rows)};
				memspace_id = H5Screate_simple(ndims, mem_dims, NULL);
				// select hyperslab of memory, for row by row traversal
				hsize_t mstart[2] = {0, 0};  // element offset for first block
				hsize_t mcount[2] = {static_cast<hsize_t>(cols), 1}; // # of blocks
				hsize_t mstride[2] = {1, static_cast<hsize_t>(rows)};  // element stride to get to next block
				hsize_t mblock[2] = {1, static_cast<hsize_t>(rows)};  // block size  rows x 1
				// hsize_t mcount[2] = {1, static_cast<hsize_t>(rows)}; // # of blocks
				// hsize_t mstride[2] = {static_cast<hsize_t>(cols), 1};  // element stride to get to next block
				// hsize_t mblock[2] = {static_cast<hsize_t>(cols), 1};  // block size  rows x 1
				H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, mstart, mstride, mcount, mblock);
			}

			// float type
			utils::hdf5::datatype<T> type;
			hid_t type_id = type.value;
			// hid_t file_type_id = H5Dget_type(dataset_id);
			// FMT_PRINT("mem type {} file type {}.  same? {}\n", type_id, file_type_id, H5Tequal(type_id, file_type_id) ? "yes" : "no");				
			
			// read data.  use mem_type is variable length string.  memspace is same as file space. 
			//auto status = 
			H5Dread(dataset_id, type_id, memspace_id, filespace_id, H5P_DEFAULT, vectors);
			// H5Tclose(file_type_id);
			H5Sclose(memspace_id);
			H5Sclose(filespace_id);
			H5Dclose(dataset_id);

			return true;
		}



	public:

		HDF5Reader(std::string const & _filename) : filename(_filename) {}
		virtual ~HDF5Reader() {}

			// // open the file for reading only.
			// fid = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);


		/*get gene expression matrix size*/
		template <typename T>
		bool getVectorSize(std::string const & path, T& count) {
			// auto stime = getSysTime();

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

			// auto etime = getSysTime();
			// FMT_ROOT_PRINT("get vector size {} in {} sec\n", count, get_duration_s(stime, etime));
			return res;
		}
		template <typename T, typename S>
		bool getSparseMatrixSize(std::string const & path, T& rows, T& cols, S& nelements, std::string &format) {
			// auto stime = getSysTime();

			// open the file for reading only.
			hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
			if (file_id < 0) return false;
			
			bool res = readSparseMatrixSize(file_id, path.c_str(), rows, cols, nelements);

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
			readStringAttribute(group_id, "format", format);
			
			H5Gclose(group_id);

			H5Fclose(file_id);

			// auto etime = getSysTime();
			// FMT_ROOT_PRINT("get sparse matrix size {}x{}, nz {} in {} sec\n", rows, cols, nelements, get_duration_s(stime, etime));
			return res;

		}
		template <typename T>
		bool getMatrixSize(std::string const & path, T& rows, T& cols) {
			// auto stime = getSysTime();

			// open the file for reading only.
			hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
			if (file_id < 0) return false;
			
			bool res = readMatrixSize(file_id, path.c_str(), rows, cols);

			H5Fclose(file_id);

			// auto etime = getSysTime();
			// FMT_ROOT_PRINT("get matrix size {}x{} in {} sec\n", rows, cols, get_duration_s(stime, etime));
			return res;
		}

		template <typename T, typename S>
		bool loadVectorData(std::string const & path, 
			S const & count, 
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

		template <typename T, typename I, typename P, typename S, typename SS>
		bool loadSparseMatrixData(std::string const & path,
			S const & rows, S const & cols, SS const & nz, std::string const & format, 
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
			size_t nptrs = 0;
			if (format == "csc")  nptrs = cols + 1;
			else if (format == "csr")  nptrs = rows + 1;
			readVector(group_id, "p", nptrs, p);
			
			H5Gclose(group_id);
			H5Fclose(file_id);
			return true;
		}


		template <typename T, typename S>
		bool loadMatrixData(std::string const & path, 
			S const & rows, S const & cols, bool & rowMajor, 
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

			// // read the names.
			// readStrings(group_id, "axis1", genes);
			// readStrings(group_id, "axis0", samples);

			// read the data.
			std::string format;
			readStringAttribute(group_id, "order", format);
			rowMajor = (format == "C");
			
			readMatrix(group_id, "block0_values", rows, cols, rowMajor, vectors);

			H5Gclose(group_id);
			H5Fclose(file_id);
			return true;
		}


};
}