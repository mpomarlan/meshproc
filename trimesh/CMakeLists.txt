cmake_minimum_required(VERSION 2.8.3)
project(trimesh)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Boost REQUIRED system filesystem date_time thread iostreams)
find_package(catkin REQUIRED COMPONENTS )

catkin_package(
  INCLUDE_DIRS
    include
  LIBRARIES
  CATKIN_DEPENDS
  DEPENDS
    Boost
    )

include_directories(include
                   )

include_directories(${catkin_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})
link_directories(${catkin_LIBRARY_DIRS})

set(libsrc
    libsrc/conn_comps.cc
    libsrc/diffuse.cc
    libsrc/edgeflip.cc
    libsrc/faceflip.cc
    libsrc/filter.cc
    libsrc/GLCamera.cc
    libsrc/ICP.cc
    libsrc/KDtree.cc
    libsrc/lmsmooth.cc
    libsrc/overlap.cc
    libsrc/TriMesh_bounding.cc
    libsrc/TriMesh_connectivity.cc
    libsrc/TriMesh_curvature.cc
    libsrc/TriMesh_io.cc
    libsrc/TriMesh_grid.cc
    libsrc/TriMesh_normals.cc
    libsrc/TriMesh_pointareas.cc
    libsrc/TriMesh_stats.cc
    libsrc/TriMesh_tstrips.cc
    libsrc/remove.cc
    libsrc/reorder_verts.cc
    libsrc/shared.cc
    libsrc/subdiv.cc
    )

add_definitions("-Wno-unused-result")
add_library(trimesh ${libsrc})
target_link_libraries(trimesh ${catkin_LIBRARIES} ${Boost_LIBRARIES})

install(TARGETS trimesh
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

install(DIRECTORY include/
  DESTINATION include)

