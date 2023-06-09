if (BUILD_TESTS)


set(fastde_test_compile_defs "-DVERBOSE")

if (ENABLE_SANITIZER)
	set(fastde_test_compile_flags ${SANITIZER_COMPILE_FLAGS})
	set(fastde_test_link_flags ${SANITIZER_LINK_FLAGS})
else()
	set(fastde_test_compile_flags "")
	set(fastde_test_link_flags "")
endif()
set(fastde_test_link_libs "m")


#if (USE_ZLIB)
#	set(fastde_test_compile_defs "${fastde_test_compile_defs};-DCOMPRESSED_INPUT")
#	set(fastde_test_link_libs "z;${fastde_test_link_libs}")
#endif(USE_ZLIB)

#if (USE_MKL)
#	set(fastde_test_compile_defs "${fastde_test_compile_defs};-DWITH_MKL;${MKL_DEFINES}")
#	set(fastde_test_compile_flags "${MKL_COMPILE_FLAGS}")
#	set(fastde_test_link_flags "${MKL_LINK_FLAGS}")
#	include_directories(${MKL_INCLUDE_DIRS})
#	link_directories(${MKL_LIB_DIRS})
#	set(fastde_test_link_libs "${MKL_LIBS};${fastde_test_link_libs}")
#endif(USE_MKL)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
	set(fastde_test_compile_defs "${fastde_test_compile_defs};-DHAVE_INTEL_ICC")
	if (PRINT_VEC_REPORT)
		set(fastde_test_compile_flags "${fastde_test_compile_flags};-qopt-report=5;-qopt-report-phase=vec;-qopt-report-file=stdout")
	endif(PRINT_VEC_REPORT)

elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	if (PRINT_VEC_REPORT)
		set(fastde_test_compile_flags "${fastde_test_compile_flags};-ftree-vectorizer-verbose=4")
	endif(PRINT_VEC_REPORT)

endif()


set(TEST_APPS gen_data_io)

foreach( exec IN LISTS TEST_APPS)
	add_executable(${exec} ${exec}.cpp)
	foreach (cdef IN LISTS fastde_test_compile_defs)
		target_compile_definitions(${exec} PRIVATE ${cdef})
	endforeach(cdef)
	foreach (cflgs IN LISTS fastde_test_compile_flags)
		target_compile_options(${exec} PRIVATE ${cflgs})
	endforeach(cflgs)
	target_link_libraries(${exec} ${fastde_test_link_flags} ${fastde_test_link_libs} ${EXTRA_LIBS})

	target_compile_definitions(${exec} PRIVATE -DFDE_TEST_DATA_DIR="${PROJECT_SOURCE_DIR}/tests/data")

	if (USE_MPI)
		target_compile_definitions(${exec} PRIVATE -DUSE_MPI)
		# MPI_LIBRARIES already in EXTRA_LIBS
	endif(USE_MPI)
endforeach(exec)


set(CATCH2_TEST_APPS test_io test_sparsemat_transpose test_sparsemat_todense test_sparsemat_sums test_sparsemat_bind test_normalize test_foldchange test_wmwtest test_ttest)

foreach( exec IN LISTS CATCH2_TEST_APPS)
	add_executable(${exec} ${exec}.cpp)
	foreach (cdef IN LISTS fastde_test_compile_defs)
		target_compile_definitions(${exec} PRIVATE ${cdef})
	endforeach(cdef)
	foreach (cflgs IN LISTS fastde_test_compile_flags)
		target_compile_options(${exec} PRIVATE ${cflgs})
	endforeach(cflgs)
	target_link_libraries(${exec} PRIVATE Catch2::Catch2WithMain ${fastde_test_link_flags} ${fastde_test_link_libs} ${EXTRA_LIBS} )

	target_compile_definitions(${exec} PRIVATE -DFDE_TEST_DATA_DIR="${PROJECT_SOURCE_DIR}/tests/data")

	if (USE_MPI)
		target_compile_definitions(${exec} PRIVATE -DUSE_MPI)
		# MPI_LIBRARIES already in EXTRA_LIBS
	endif(USE_MPI)
endforeach(exec)





endif(BUILD_TESTS)


