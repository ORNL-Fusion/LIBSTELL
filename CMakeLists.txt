cmake_minimum_required (VERSION 3.14)

project (libstell C Fortran)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package (NetCDF MODULE REQUIRED)
find_package (SCALAPACK MODULE REQUIRED)
find_package (MPI REQUIRED)
find_package (BLAS REQUIRED)
find_package (LAPACK REQUIRED)

add_library (stell STATIC)

add_subdirectory (Sources)

OPTION (USE_PROFILER "Activate code profiling." OFF)

# Activate or Deactivate fast copies. The fast copies maybe activated by setting
# -DUSE_FAST_COPY=ON on the cmake command line.
OPTION (USE_FAST_COPY "Activate code fast copies." OFF)

#  Fortran specific settings. The first setting tells the compiler to use the C
#  preprocessor.
target_compile_options (stell PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:-cpp>)
target_compile_definitions (stell PUBLIC
                            $<$<PLATFORM_ID:Darwin>:DARWIN>
                            $<$<PLATFORM_ID:Linux>:LINUX>
                            $<$<PLATFORM_ID:Windows>:WIN32>
                            $<$<BOOL:${USE_FAST_COPY}>:FAST_COPY>
                            $<$<BOOL:${USE_PROFILER}>:PROFILE_ON>
                            $<$<BOOL:${NetCDF_FOUND}>:NETCDF>
                            $<$<BOOL:${MPI_Fortran_FOUND}>:MPI_OPT>)
target_link_libraries (stell PUBLIC
                       $<$<BOOL:${MPI_Fortran_FOUND}>:MPI::MPI_Fortran>
                       $<$<BOOL:${NetCDF_FOUND}>:NetCDF::NetCDF>
                       $<$<BOOL:${SCALAPACK_FOUND}>:SCALAPACK::SCALAPACK>
                       $<$<BOOL:${BLAS_FOUND}>:${BLAS_LIBRARIES}>
                       $<$<BOOL:${LAPACK_FOUND}>:${LAPACK_LIBRARIES}>)
#target_include_directories (stell PRIVATE
#                            ${MPI_C_INCLUDE_DIRS}
#                            ${MPI_Fortran_INCLUDE_DIRS})
if (MPI_Fortran_FOUND)
    set_target_properties (stell PROPERTIES
                           MPIEXEC_EXECUTABLE "${MPIEXEC_EXECUTABLE}"
                           MPIEXEC_NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
                           MPIEXEC_MAX_NUMPROCS "${MPIEXEC_MAX_NUMPROCS}")
endif ()

################################################################################
#  Testing                                                                     #
################################################################################

#  Build test utilities.
add_executable (xstell_unit_test_runner)
add_dependencies (xstell_unit_test_runner stell)

target_link_libraries (xstell_unit_test_runner stell)
target_include_directories (xstell_unit_test_runner PUBLIC $<TARGET_PROPERTY:stell,BINARY_DIR>)

#  Build file_opts testers
add_executable (xstell_file_opts_runner)
add_dependencies (xstell_file_opts_runner stell)

target_link_libraries (xstell_file_opts_runner stell)
target_include_directories (xstell_file_opts_runner PUBLIC $<TARGET_PROPERTY:stell,BINARY_DIR>)

add_subdirectory (Testing)
