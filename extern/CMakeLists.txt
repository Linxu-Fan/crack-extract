
#
# eigen
#
set (EIGEN_INCLUDE_DIR ${EXTERN_PROJECTS_DIR}/eigen CACHE STRING "eigen includes")



#
# trimesh
#
set(TRIMESH_PROJECT_DIR ${EXTERN_PROJECTS_DIR}/trimesh2 CACHE STRING "location of trimesh2" )
set (TRIMESH_INCLUDE_DIRS ${TRIMESH_PROJECT_DIR}/include CACHE INTERNAL "trimesh2 include directories")
set (TRIMESH_LIBRARIES trimesh2 CACHE INTERNAL  "trimesh2 libraries")
set ( SRC_DIR ${TRIMESH_PROJECT_DIR}/libsrc)

add_library(trimesh2 STATIC
  ${SRC_DIR}/TriMesh_bounding.cc
	${SRC_DIR}/TriMesh_connectivity.cc
  ${SRC_DIR}/TriMesh_curvature.cc
  ${SRC_DIR}/TriMesh_io.cc
  ${SRC_DIR}/TriMesh_grid.cc
  ${SRC_DIR}/TriMesh_normals.cc
  ${SRC_DIR}/TriMesh_pointareas.cc
  ${SRC_DIR}/TriMesh_stats.cc
  ${SRC_DIR}/TriMesh_tstrips.cc
  ${SRC_DIR}/GLCamera.cc
  ${SRC_DIR}/ICP.cc
  ${SRC_DIR}/KDtree.cc
  ${SRC_DIR}/conn_comps.cc
  ${SRC_DIR}/diffuse.cc
  ${SRC_DIR}/edgeflip.cc
  ${SRC_DIR}/faceflip.cc
  ${SRC_DIR}/filter.cc
  ${SRC_DIR}/overlap.cc
  ${SRC_DIR}/remove.cc
  ${SRC_DIR}/reorder_verts.cc
  ${SRC_DIR}/subdiv.cc
)

target_include_directories(trimesh2 PUBLIC ${TRIMESH_INCLUDE_DIRS})

#
# libigl
#
set(LIBIGL_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libigl)
#add_subdirectory(${LIBIGL_ROOT_DIR})

#set (LIBIGL_LIBRARIES igl::core CACHE STRING "linigl libs")
set (LIBIGL_INCLUDE_DIR ${LIBIGL_ROOT_DIR}/include CACHE STRING "libigl includes")


#
# OpenVDB
#
add_subdirectory(openvdb)

set (OPENVDB_PROJECT_DIR ${EXTERN_PROJECTS_DIR}/openvdb CACHE STRING "location of openvdb" )
set (OPENVDB_INCLUDE_DIRS ${OPENVDB_PROJECT_DIR} CACHE INTERNAL "openvdb include directories")
set (OPENVDB_LIBRARIES openvdb CACHE INTERNAL  "openvdb libraries")


#
# Voro++
#

add_subdirectory(voro)
set(VORO_PLUS_PLUS_INCLUDE_DIRS ${EXTERN_PROJECTS_DIR}/voro/src CACHE INTERNAL "") 
set(VORO_PLUS_PLUS_LIBRARIES voro++ CACHE INTERNAL "") 

#
# MCUT
#

set (MCUT_BUILD_WITH_ARBITRARY_PRECISION_NUMBERS ON CACHE INTERNAL "")
set (MCUT_BUILD_WITH_MULTITHREADING OFF CACHE INTERNAL "")

add_subdirectory(mcut)
set(MCUT_INCLUDE_DIRS ${EXTERN_PROJECTS_DIR}/mcut/include CACHE INTERNAL "") 
set(MCUT_LIBRARIES mcut CACHE INTERNAL "") 

#
# Earcut (triangulation lib)
#

set(EARCUT_INCLUDE_DIRS ${EXTERN_PROJECTS_DIR}/earcut.hpp/include CACHE INTERNAL "") 


