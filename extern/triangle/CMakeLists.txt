cmake_minimum_required(VERSION 2.8.12)
project(triangle)

set(TRIANGLE_LIB_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include CACHE INTERNAL "") 
set(TRIANGLE_LIB_LIBRARIES triangle CACHE INTERNAL "") 

add_library(triangle STATIC ${CMAKE_CURRENT_SOURCE_DIR}/triangle.c)
target_compile_definitions(triangle PRIVATE -DTRILIBRARY -DANSI_DECLARATORS)
if(WIN32)
   target_compile_definitions(triangle PRIVATE -DNO_TIMER)
endif()

target_include_directories(triangle PRIVATE ${TRIANGLE_LIB_INCLUDE_DIRS})

set ( INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include )

# Generate position independent code
# set_target_properties(triangle PROPERTIES POSITION_INDEPENDENT_CODE ON)
