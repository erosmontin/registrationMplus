cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(RegistrationSuite LANGUAGES CXX)

# — ITK: grab everything you built under your ITK_DIR
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

message(STATUS ">>> ITK_FOUND       = ${ITK_FOUND}")
message(STATUS ">>> ITK_VERSION     = ${ITK_VERSION}")
message(STATUS ">>> ITK_LIBRARIES   = ${ITK_LIBRARIES}")


# — Boost: for ProgramOptions
find_package(Boost REQUIRED COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIRS})

# — Settings
set(CMAKE_CXX_STANDARD        14)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
add_compile_options(-fpermissive)

# — Sub‐folders
add_subdirectory(3DRegistration)

