#[==[
Provides the following variables:

  * `SCALAPACK_FOUND`: Whether SCALAPACK was found or not.
  * `SCALAPACK_LIBRARIES`: Libraries necessary to use SCALAPACK.
  * `SCALAPACK::SCALAPACK`: A target to use with `target_link_libraries`.
#]==]

# Find it manually.
find_library (SCALAPACK_LIBRARY
              NAMES scalapack scalapack-openmpi
              DOC "scalapack library")
mark_as_advanced (SCALAPACK_LIBRARY)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (SCALAPACK
                                   REQUIRED_VARS SCALAPACK_LIBRARY)

if (SCALAPACK_FOUND)
    set (SCALAPACK_LIBRARIES ${SCALAPACK_LIBRARY})

    if (NOT TARGET SCALAPACK::SCALAPACK)
       add_library (SCALAPACK::SCALAPACK UNKNOWN IMPORTED)
       set_target_properties (SCALAPACK::SCALAPACK PROPERTIES
                              IMPORTED_LOCATION ${SCALAPACK_LIBRARY})
    endif ()
endif ()
