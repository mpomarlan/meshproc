cmake_minimum_required(VERSION 2.8.3)
project(meshproc_msgs)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(MSG_DEPS
  std_msgs
  shape_msgs
)

set(MSG_FILES
  Dummy.msg
)

set(SRV_FILES
  AffineTransformMesh.srv
  GetMesh.srv
  GetMeshAABB.srv
  GetNearMeshVertices.srv
  LoadMesh.srv
  CSGRequest.srv
  UnloadMesh.srv
  GetMeshProps.srv
  GetLoadedMeshNames.srv
)

set(ACT_FILES
)

find_package(catkin REQUIRED COMPONENTS
  ${MSG_DEPS}
  message_generation
)

#add_action_files(DIRECTORY action FILES ${ACT_FILES})
add_message_files(DIRECTORY msg FILES ${MSG_FILES})
add_service_files(DIRECTORY srv FILES ${SRV_FILES})

generate_messages(DEPENDENCIES ${MSG_DEPS})

catkin_package(DEPENDS ${MSG_DEPS} message_runtime)
