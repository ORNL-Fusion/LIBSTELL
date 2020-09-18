target_sources(stell
               PRIVATE 
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/copy_file.c> 
)
