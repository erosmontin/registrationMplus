cmake_minimum_required(VERSION 2.4)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()

project(PRJ_3DRegCalibrate)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

find_package(Boost REQUIRED COMPONENTS program_options)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(3DRegCalibrate 3DRegCalibrate.cxx)
    target_link_libraries(3DRegCalibrate ${Boost_LIBRARIES} ${ITK_LIBRARIES})
endif()
