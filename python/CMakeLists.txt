# ── Python pybind11 bindings for Mplus ────────────────────────────────────
#
# Included from the root CMakeLists.txt when BUILD_PYTHON=ON.
# Also usable via scikit-build-core (see pyproject.toml).

find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_core
  mplus/_core.cpp
)

target_include_directories(_core PRIVATE
  ${CMAKE_SOURCE_DIR}/includes
  ${CMAKE_SOURCE_DIR}/Metrics/Mplus
  ${CMAKE_SOURCE_DIR}/Metrics/NGF/NGFImageMetric/NGFImageToImageMetric/Code
)

target_link_libraries(_core PRIVATE
  ${ITK_LIBRARIES}
)

if(USE_CUDA)
  target_link_libraries(_core PRIVATE mplus_cuda)
  target_compile_definitions(_core PRIVATE USE_CUDA)
endif()

# Install the extension module next to the Python package
install(TARGETS _core DESTINATION mplus)
