target_sources(stell
               PRIVATE
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/gmresr.f>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/dPackgmres.f>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/gmres_lib.f>
)
