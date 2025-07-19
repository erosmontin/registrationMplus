cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(3DRegAffine LANGUAGES CXX)

# — Find ITK (all modules you built) —
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

# — Find Boost.ProgramOptions —
find_package(Boost REQUIRED COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIRS})

# — Compiler settings —
set(CMAKE_CXX_STANDARD 14)
# allow non‑conforming qualifiers in Mplus
set(PERMISSIVE_FLAGS -fpermissive)

# — 3DRegSimilarity executable —
add_executable(3DRegSimilarity
  3DRegSimilarity.cxx
)
target_compile_options(3DRegSimilarity PRIVATE ${PERMISSIVE_FLAGS})
target_include_directories(3DRegSimilarity PRIVATE
  ${CMAKE_SOURCE_DIR}/3DRegistration/Metrics/Mplus
  ${CMAKE_SOURCE_DIR}/3DRegistration/Metrics/NGF/NGFImageToImageMetric/Code
  ${CMAKE_SOURCE_DIR}/includes
)
target_link_libraries(3DRegSimilarity PRIVATE
  ${ITK_LIBRARIES}
  ${Boost_LIBRARIES}
)

# — 3DRegAffineMultiLevel executable —
add_executable(3DRegAffineMultiLevel
  3DRegAffineMultiLevel.cxx
)
target_compile_options(3DRegAffineMultiLevel PRIVATE ${PERMISSIVE_FLAGS})
target_include_directories(3DRegAffineMultiLevel PRIVATE
  ${CMAKE_SOURCE_DIR}/3DRegistration/Metrics/Mplus
  ${CMAKE_SOURCE_DIR}/3DRegistration/Metrics/NGF/NGFImageToImageMetric/Code
  ${CMAKE_SOURCE_DIR}/includes
)
target_link_libraries(3DRegAffineMultiLevel PRIVATE
  ${ITK_LIBRARIES}
  ${Boost_LIBRARIES}
)
