# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 2.4)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()


# This project is designed to be built outside the Insight source tree.
project(PRJ_ReadTXTSpoint)

# Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(ReadTXTSpoint ReadTXTSpoint.cxx )

target_link_libraries(ReadTXTSpoint ${ITK_LIBRARIES}) 
