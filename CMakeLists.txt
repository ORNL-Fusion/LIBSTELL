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

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    target_compile_definitions (stell PUBLIC DARWIN)
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux" OR 
        ${CMAKE_SYSTEM_NAME} MATCHES "CrayLinuxEnvironment")
    target_compile_definitions (stell PUBLIC LINUX)
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Windows" OR 
        ${CMAKE_SYSTEM_NAME} MATCHES "MSYS")
    target_compile_definitions (stell PUBLIC WIN32)
endif ()

OPTION (USE_PROFILER "Activate code profiling." OFF)
if (USE_PROFILER)
    target_compile_definitions (stell PUBLIC PROFILE_ON)
endif ()

# Activate or Deactivate fast copies. The fast copies maybe activated by setting
# -DUSE_FAST_COPY=ON on the cmake command line.
OPTION (USE_FAST_COPY "Activate code fast copies." OFF)
if (USE_FAST_COPY)
    target_compile_definitions (stell PUBLIC FAST_COPY)
endif ()

#  Fortran specific settings. The first setting tells the compiler to use the C
#  preprocessor.
target_compile_options (stell PUBLIC -cpp)

if (MPI_Fortran_FOUND)
    target_compile_options (stell PUBLIC -DMPI_OPT)
    target_include_directories (stell PRIVATE ${MPI_C_INCLUDE_DIRS})
    target_link_libraries (stell INTERFACE MPI::MPI_Fortran)
    set_target_properties (stell PROPERTIES
                           MPIEXEC_EXECUTABLE "${MPIEXEC_EXECUTABLE}"
                           MPIEXEC_NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
                           MPI_OVERSUBSCRIBE_FLAG "${MPI_OVERSUBSCRIBE_FLAG}")
endif ()

if (NetCDF_FOUND)
    target_compile_definitions (stell PUBLIC NETCDF)
    target_link_libraries (stell PUBLIC NetCDF::NetCDF)
endif ()

if (SCALAPACK_FOUND)
    target_link_libraries (stell PUBLIC SCALAPACK::SCALAPACK)
endif ()

if (BLAS_FOUND)
    target_link_libraries (stell PUBLIC ${BLAS_LIBRARIES})
endif ()

if (LAPACK_FOUND)
    target_link_libraries (stell PUBLIC ${LAPACK_LIBRARIES})
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
