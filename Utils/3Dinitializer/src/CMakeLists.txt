# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 2.4)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()

project(PRJ_3DInitializer)

# Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(3Dinitializer 3Dinitializer.cxx )

target_link_libraries(3Dinitializer ${ITK_LIBRARIES}) 
