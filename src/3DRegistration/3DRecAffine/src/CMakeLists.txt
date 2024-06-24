# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 2.4)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()


# This project is designed to be built outside the Insight source tree.
project(PRJ_3DRecBsplines_new)

# Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})


find_package(Boost REQUIRED COMPONENTS program_options)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(3DRecBsplines_new 3DRecBsplines_new.cxx)
    target_link_libraries(3DRecBsplines_new ${Boost_LIBRARIES}  ${ITK_LIBRARIES})
endif()