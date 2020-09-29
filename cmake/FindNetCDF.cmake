#[==[
Provides the following variables:

  * `NetCDF_FOUND`: Whether NetCDF was found or not.
  * `NetCDF_INCLUDE_DIRS`: Include directories necessary to use NetCDF.
  * `NetCDF_LIBRARIES`: Libraries necessary to use NetCDF.
  * `NetCDF_VERSION`: The version of NetCDF found.
  * `NetCDF::NetCDF`: A target to use with `target_link_libraries`.
#]==]

# Try to find a CMake-built NetCDF.
find_package (NetCDF CONFIG QUIET)
if (NetCDF_FOUND)
# Forward the variables in a consistent way.
    set (NetCDF_FOUND "${NetCDF_FOUND}")
    set (NetCDF_INCLUDE_DIRS "${NetCDF_INCLUDE_DIR}")
    set (NetCDF_LIBRARIES "${NetCDF_LIBRARIES}")
    set (NetCDF_VERSION "${NetCDFVersion}")
    if (NOT TARGET NetCDF::NetCDF)
        add_library(NetCDF::NetCDF INTERFACE IMPORTED)
        if (TARGET "netCDF::netcdf")
# 4.7.3
           set_target_properties (NetCDF::NetCDF PROPERTIES
                                  INTERFACE_LINK_LIBRARIES "netCDF::netcdf")
        elseif (TARGET "netcdf netcdff")
            set_target_properties(NetCDF::NetCDF PROPERTIES
                                  INTERFACE_LINK_LIBRARIES "netcdf netcdff")
        else ()
            set_target_properties(NetCDF::NetCDF PROPERTIES
                                  INTERFACE_LINK_LIBRARIES "${netCDF_LIBRARIES}")
        endif ()
    endif ()
# Skip the rest of the logic in this file.
    return ()
endif ()

find_package (PkgConfig QUIET)
if (PkgConfig_FOUND)
    pkg_check_modules(_NetCDF QUIET netcdf netcdff IMPORTED_TARGET)
    if (_NetCDF_FOUND)
# Forward the variables in a consistent way.
        set (NetCDF_FOUND "${_NetCDF_FOUND}")
        set (NetCDF_INCLUDE_DIRS "${_NetCDF_INCLUDE_DIRS}")
        set (NetCDF_LIBRARIES "${_NetCDF_LIBRARIES}")
        set (NetCDF_VERSION "${_NetCDF_VERSION}")
        if (NOT TARGET NetCDF::NetCDF)
            add_library (NetCDF::NetCDF INTERFACE IMPORTED)
            set_target_properties (NetCDF::NetCDF PROPERTIES
                                   INTERFACE_LINK_LIBRARIES "PkgConfig::_NetCDF")
        endif ()
# Skip the rest of the logic in this file.
        return ()
    endif ()
endif ()

# Find it manually.
find_path (NetCDF_INCLUDE_DIR
           NAMES netcdf.h
           DOC "netcdf include directories")
mark_as_advanced (NetCDF_INCLUDE_DIR)

find_library (NetCDF_C_LIBRARY
              NAMES netcdf
              DOC "netcdf C library")
mark_as_advanced (NetCDF_LIBRARY)
find_library (NetCDF_Fortran_LIBRARY
              NAMES netcdff
              DOC "netcdf Fortran library")
mark_as_advanced (NetCDF_LIBRARY)

if (NetCDF_INCLUDE_DIR)
    file (STRINGS "${NetCDF_INCLUDE_DIR}/netcdf_meta.h" _netcdf_version_lines
          REGEX "#define[ \t]+NC_VERSION_(MAJOR|MINOR|PATCH|NOTE)")
    string (REGEX REPLACE ".*NC_VERSION_MAJOR *\([0-9]*\).*" "\\1" _netcdf_version_major "${_netcdf_version_lines}")
    string (REGEX REPLACE ".*NC_VERSION_MINOR *\([0-9]*\).*" "\\1" _netcdf_version_minor "${_netcdf_version_lines}")
    string (REGEX REPLACE ".*NC_VERSION_PATCH *\([0-9]*\).*" "\\1" _netcdf_version_patch "${_netcdf_version_lines}")
    string (REGEX REPLACE ".*NC_VERSION_NOTE *\"\([^\"]*\)\".*" "\\1" _netcdf_version_note "${_netcdf_version_lines}")
    set (NetCDF_VERSION "${_netcdf_version_major}.${_netcdf_version_minor}.${_netcdf_version_patch}${_netcdf_version_note}")
    unset (_netcdf_version_major)
    unset (_netcdf_version_minor)
    unset (_netcdf_version_patch)
    unset (_netcdf_version_note)
    unset (_netcdf_version_lines)
endif ()

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args(NetCDF
                                  REQUIRED_VARS NetCDF_C_LIBRARY NetCDF_Fortran_LIBRARY NetCDF_INCLUDE_DIR
                                  VERSION_VAR NetCDF_VERSION)

if (NetCDF_FOUND)
    set (NetCDF_INCLUDE_DIRS ${NetCDF_INCLUDE_DIR})
    set (NetCDF_LIBRARIES "${NetCDF_C_LIBRARY} ${NetCDF_Fortran_LIBRARY}")

    if (NOT TARGET NetCDF::NetCDF)
        add_library (NetCDF::NetCDF UNKNOWN IMPORTED)
        set_target_properties (NetCDF::NetCDF PROPERTIES
                               IMPORTED_LOCATION "${NetCDF_Fortran_LIBRARY}"
                               INTERFACE_INCLUDE_DIRECTORIES ${NetCDF_INCLUDE_DIR})
    endif ()
endif ()
