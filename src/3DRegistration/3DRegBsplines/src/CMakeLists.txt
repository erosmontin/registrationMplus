cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(3DRegBsplines LANGUAGES CXX)

# — ITK (all modules you built) —
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

# — Boost.ProgramOptions —
find_package(Boost REQUIRED COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIRS})

# — Compiler settings —
set(CMAKE_CXX_STANDARD 14)
# allow non‑conforming qualifiers in Mplus (if your Bspline code shares it)
set(PERMISSIVE_FLAGS -fpermissive)

# — 3DRegBsplines executable —
add_executable(3DRegBsplines
  3DRegBsplines.cxx
)
target_compile_options(3DRegBsplines PRIVATE ${PERMISSIVE_FLAGS})
target_link_libraries(3DRegBsplines PRIVATE
  ${ITK_LIBRARIES}
  ${Boost_LIBRARIES}
)
