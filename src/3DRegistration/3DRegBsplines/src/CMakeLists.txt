# ── 3DRegBsplines ──────────────────────────────────────────────────────────

find_package(Boost REQUIRED COMPONENTS program_options)

add_executable(3DRegBsplines 3DRegBsplines.cxx)
target_include_directories(3DRegBsplines PRIVATE
  ${Boost_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/includes
  ${CMAKE_SOURCE_DIR}/Metrics/Mplus
)
target_link_libraries(3DRegBsplines
  PRIVATE
    ${Boost_LIBRARIES}
    ${ITK_LIBRARIES}
)
if(USE_CUDA)
  target_link_libraries(3DRegBsplines PRIVATE mplus_cuda)
endif()