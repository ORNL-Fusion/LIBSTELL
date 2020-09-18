target_sources(stell
               PRIVATE
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/smpint.f>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/dcuhre.f>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/dqag.f>
)
