# ── CUDA acceleration library for Mplus ───────────────────────────────────
# This CMakeLists.txt is included only when USE_CUDA=ON.

cmake_minimum_required(VERSION 3.18)

# Configuration header
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/MplusCudaConfig.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/MplusCudaConfig.h
)

# ── CUDA library ──────────────────────────────────────────────────────────
add_library(mplus_cuda STATIC
  DistanceTransform.cu
  LabelMetricKernels.cu
  DerivativeOps.cu
  DeviceMemoryPool.cu
)

target_include_directories(mplus_cuda PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}   # for MplusCudaConfig.h
)

target_link_libraries(mplus_cuda
  PUBLIC
    CUDA::cudart
)

# Separable compilation so device code can be linked from C++ TUs
set_target_properties(mplus_cuda PROPERTIES
  CUDA_SEPARABLE_COMPILATION ON
  POSITION_INDEPENDENT_CODE ON
  CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}"
)
