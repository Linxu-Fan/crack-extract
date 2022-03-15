# Works with 3.11 and tested through 3.15
cmake_minimum_required(VERSION 3.11...3.16)

# Project name and a few useful settings. Other commands can pick up the results
project(crackExtraction
    VERSION 0.1
    DESCRIPTION "Crack extraction project"
    LANGUAGES CXX)
    
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
set(EXTERN_PROJECTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern)

add_subdirectory(${EXTERN_PROJECTS_DIR})

set (INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include )
list (APPEND EXTRA_LIBS  ${MCUT_LIBRARIES} ${TRIMESH_LIBRARIES} ${OPENVDB_LIBRARIES} ${VORO_PLUS_PLUS_LIBRARIES})
list (APPEND EXTRA_INCLUDES ${EXTERN_PROJECTS_DIR}/eigen  ${MCUT_INCLUDE_DIRS}  ${TRIMESH_INCLUDE_DIRS} ${INCLUDE_DIR} ${OPENVDB_INCLUDE_DIRS}  ${VORO_PLUS_PLUS_INCLUDE_DIRS} ${EARCUT_INCLUDE_DIRS} ${LIBIGL_INCLUDE_DIR})

message(STATUS "EXTRA_LIBS = ${EXTRA_LIBS}")
message(STATUS "EXTRA_INCLUDES = ${EXTRA_INCLUDES}")

add_executable(
    crackExtraction 
    src/main.cpp 
    src/utility-functions.cpp
    src/damageGradient.cpp 
    src/extractCrack.cpp 
    src/weights.cpp
    src/polygon_triangulate.cpp)

if(WIN32)
    target_compile_options(crackExtraction PRIVATE /bigobj)
else()
    target_compile_options(crackExtraction PRIVATE -fmax-errors=1 ) # stop compilation after the first error (makes life a little easier)
endif()
target_compile_features(crackExtraction PRIVATE cxx_std_14 )
target_include_directories(crackExtraction PRIVATE ${EXTRA_INCLUDES})
target_link_libraries (crackExtraction ${EXTRA_LIBS}) 
target_compile_definitions(crackExtraction PUBLIC -DROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}" )
