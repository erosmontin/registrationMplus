# This is the root ITK CMakeLists file.

if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()

cmake_minimum_required(VERSION 2.8)

project(2DFFDMANGF)

# Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(2DFFDMANGF 2DFFDMANGF.cxx )

target_link_libraries(2DFFDMANGF ${ITK_LIBRARIES}) 
