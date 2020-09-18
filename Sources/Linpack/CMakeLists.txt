target_sources(stell
               PRIVATE
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/liprec.f
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/li_gesl.f
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/li_gefa.f
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/li_gbsl.f
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/li_gbfa.f 
)
