# ── 3DRegAffine executables ────────────────────────────────────────────────

find_package(Boost REQUIRED COMPONENTS program_options)

set(_COMMON_INCLUDES
  ${Boost_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/includes
  ${CMAKE_SOURCE_DIR}/Metrics/Mplus
)
set(_COMMON_LIBS
  ${Boost_LIBRARIES}
  ${ITK_LIBRARIES}
)

add_executable(3DRegAffine 3DRegAffine.cxx)
target_include_directories(3DRegAffine PRIVATE ${_COMMON_INCLUDES})
target_link_libraries(3DRegAffine PRIVATE ${_COMMON_LIBS})

add_executable(3DRegAffineMultiLevel 3DRegAffineMultiLevel.cxx)
target_include_directories(3DRegAffineMultiLevel PRIVATE ${_COMMON_INCLUDES})
target_link_libraries(3DRegAffineMultiLevel PRIVATE ${_COMMON_LIBS})

add_executable(3DRegSimilarity 3DRegSimilarity.cxx)
target_include_directories(3DRegSimilarity PRIVATE ${_COMMON_INCLUDES})
target_link_libraries(3DRegSimilarity PRIVATE ${_COMMON_LIBS})

if(USE_CUDA)
  target_link_libraries(3DRegAffine PRIVATE mplus_cuda)
  target_link_libraries(3DRegAffineMultiLevel PRIVATE mplus_cuda)
  target_link_libraries(3DRegSimilarity PRIVATE mplus_cuda)
endif()