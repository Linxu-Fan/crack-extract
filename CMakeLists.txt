# Works with 3.11 and tested through 3.15
cmake_minimum_required(VERSION 3.11...3.16)

# Project name and a few useful settings. Other commands can pick up the results
project(crackExtraction
    VERSION 0.1
    DESCRIPTION "MPM Fracture Project"
    LANGUAGES CXX)
    
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)


OPTION (USE_OpenMP "Use OpenMP" ON)
IF(USE_OpenMP)
  FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  ENDIF()
ENDIF()


set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)



list(APPEND EXTRA_INCLUDES ${EIGEN3_INCLUDE_DIR})

set(EXTERN_PROJECTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern)

add_subdirectory(${EXTERN_PROJECTS_DIR})



set (INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include )
list (APPEND EXTRA_LIBS ${TRIANGLE_LIB_LIBRARIES} ${MCUT_LIBRARIES} ${PYSTRING_LIBRARIES} ${TRIMESH_LIBRARIES} ${OPENVDB_LIBRARIES} ${VORO_PLUS_PLUS_LIBRARIES})
list (APPEND EXTRA_INCLUDES  ${TRIANGLE_LIB_INCLUDE_DIRS} ${MCUT_INCLUDE_DIRS}  ${PYSTRING_INCLUDE_DIRS} ${TRIMESH_INCLUDE_DIRS} ${INCLUDE_DIR} ${OPENVDB_INCLUDE_DIRS} ${MESH_DECIMATION_INCLUDE_DIRS} ${VORO_PLUS_PLUS_INCLUDE_DIRS} ${EARCUT_INCLUDE_DIRS} ${LIBIGL_INCLUDE_DIR})

message(STATUS "EXTRA_LIBS = ${EXTRA_LIBS}")
message(STATUS "EXTRA_INCLUDES = ${EXTRA_INCLUDES}")

add_executable(
    crackExtraction 
    src/main.cpp 
    src/rigid-body-object.cpp
    src/utility-functions.cpp
    src/mpm.cpp
    src/damageGradient.cpp 
    src/extractCrack.cpp 
    src/particles.cpp 
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
