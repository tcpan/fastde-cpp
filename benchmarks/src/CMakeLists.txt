if (BUILD_BENCHMARKS)


set(fastde_benchmark_compile_defs "")
set(fastde_benchmark_compile_flags "")
set(fastde_benchmark_link_flags "")
set(fastde_benchmark_link_libs "m")


#if (USE_ZLIB)
#	set(fastde_benchmark_compile_defs "${fastde_benchmark_compile_defs};-DCOMPRESSED_INPUT")
#	set(fastde_benchmark_link_libs "z;${fastde_benchmark_link_libs}")
#endif(USE_ZLIB)

#if (USE_MKL)
#	set(fastde_benchmark_compile_defs "${fastde_benchmark_compile_defs};-DWITH_MKL;${MKL_DEFINES}")
#	set(fastde_benchmark_compile_flags "${MKL_COMPILE_FLAGS}")
#	set(fastde_benchmark_link_flags "${MKL_LINK_FLAGS}")
#	include_directories(${MKL_INCLUDE_DIRS})
#	link_directories(${MKL_LIB_DIRS})
#	set(fastde_benchmark_link_libs "${MKL_LIBS};${fastde_benchmark_link_libs}")
#endif(USE_MKL)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
	set(fastde_benchmark_compile_defs "${fastde_benchmark_compile_defs};-DHAVE_INTEL_ICC")
	if (PRINT_VEC_REPORT)
		set(fastde_benchmark_compile_flags "${fastde_benchmark_compile_flags};-qopt-report=5;-qopt-report-phase=vec;-qopt-report-file=stdout")
	endif(PRINT_VEC_REPORT)

elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	if (PRINT_VEC_REPORT)
		set(fastde_benchmark_compile_flags "${fastde_benchmark_compile_flags};-ftree-vectorizer-verbose=4")
	endif(PRINT_VEC_REPORT)

endif()


set(CORRELATION_APPS )

foreach( exec IN LISTS CORRELATION_APPS)
	add_executable(${exec} ${exec}.cpp)
	foreach (cdef IN LISTS fastde_benchmark_compile_defs)
		target_compile_definitions(${exec} PRIVATE ${cdef})
	endforeach(cdef)
	foreach (cflgs IN LISTS fastde_benchmark_compile_flags)
		target_compile_options(${exec} PRIVATE ${cflgs})
	endforeach(cflgs)
	target_link_libraries(${exec} ${fastde_benchmark_link_flags} ${fastde_benchmark_link_libs} ${EXTRA_LIBS})

	if (USE_MPI)
		target_compile_definitions(${exec} PRIVATE -DUSE_MPI)
		# MPI_LIBRARIES already in EXTRA_LIBS
	endif(USE_MPI)
endforeach(exec)




endif(BUILD_BENCHMARKS)


