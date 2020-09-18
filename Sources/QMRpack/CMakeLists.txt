target_sources(stell
               PRIVATE
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/dutfx.f>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/daxpby.f>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/drandn.f>
)
