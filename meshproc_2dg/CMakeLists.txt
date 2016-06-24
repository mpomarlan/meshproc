cmake_minimum_required(VERSION 2.8.3)
project(meshproc_2dg)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

include_directories(include)

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )
if ( NOT CGAL_FOUND )
  message(WARNING "This project requires the CGAL library, and will not be compiled.")
  return()  
endif()

# include helper file
include( ${CGAL_USE_FILE} )

find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
if (EIGEN3_FOUND)
   include( ${EIGEN3_USE_FILE} )
endif()

# Boost and its components
find_package( Boost REQUIRED )
if ( NOT Boost_FOUND )
  message(WARNING "This project requires the Boost library, and will not be compiled.")
  return()  
endif()

find_package(PkgConfig)
find_package(catkin REQUIRED COMPONENTS
  trimesh
  meshproc_csg
  meshproc_msgs
  shape_msgs
  roscpp
  roslib
)

include_directories(${catkin_INCLUDE_DIRS})
link_directories(${catkin_LIBRARY_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

catkin_package(
  INCLUDE_DIRS
    include
  LIBRARIES
    meshproc_2dg_lib
  CATKIN_DEPENDS
)

set( spec_libs
     trimesh
     GL
     GLU
     GLEW
     )

set( sources
  src/2dg.cpp)

add_library(meshproc_2dg_lib ./src/PolygonEntry.cpp)
add_dependencies(meshproc_2dg_lib ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
target_link_libraries(meshproc_2dg_lib
  ${CGAL_LIBRARIES} CGAL_Core ${CGAL_3RD_PARTY_LIBRARIES} ${EIGEN_LIBRARIES} ${catkin_LIBRARIES} ${Boost_LIBRARIES})

add_executable(${PROJECT_NAME} ${sources} ${headers})
add_to_cached_list( CGAL_EXECUTABLE_TARGETS ${PROJECT_NAME} )
add_dependencies(${PROJECT_NAME} ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
target_link_libraries(${PROJECT_NAME}
  meshproc_2dg_lib meshproc_csg_lib ${spec_libs} ${CGAL_LIBRARIES} CGAL_Core ${CGAL_3RD_PARTY_LIBRARIES} ${EIGEN_LIBRARIES} ${catkin_LIBRARIES} ${Boost_LIBRARIES})

foreach(dir
    src)
  install(DIRECTORY ${dir}
      DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
      PATTERN ".svn" EXCLUDE
      PATTERN ".git" EXCLUDE)
endforeach()

install (FILES 
    DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION})

install(TARGETS meshproc_2dg_lib
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

install(TARGETS ${PROJECT_NAME}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION})
#install(DIRECTORY include/ DESTINATION include)
install(DIRECTORY launch DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION})

