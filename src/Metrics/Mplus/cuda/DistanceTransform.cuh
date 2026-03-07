/*=========================================================================
 *  DistanceTransform.cuh  –  Pure C++ interface to GPU signed distance
 *                             transform (Felzenszwalb separable).
 *
 *  This header contains NO CUDA types — safe to #include from .cxx files.
 *=========================================================================*/
#ifndef MPLUS_DISTANCE_TRANSFORM_CUH
#define MPLUS_DISTANCE_TRANSFORM_CUH

#include <cstddef>

namespace mplus { namespace cuda {

/**
 * Compute a signed Euclidean distance transform on the GPU.
 *
 * The input is a binary float image (1.0 inside, 0.0 outside).
 * The output is a signed distance field: negative inside, positive outside.
 *
 * @param binary_in   Host pointer to the binary float image (row-major).
 * @param dist_out    Host pointer to output distance field (same size).
 * @param sizeX       Image size along X.
 * @param sizeY       Image size along Y.
 * @param sizeZ       Image size along Z.
 * @param spacingX    Voxel spacing in mm along X.
 * @param spacingY    Voxel spacing in mm along Y.
 * @param spacingZ    Voxel spacing in mm along Z.
 */
void gpu_signed_distance_transform(
    const float* binary_in,
    float*       dist_out,
    int sizeX, int sizeY, int sizeZ,
    float spacingX = 1.0f,
    float spacingY = 1.0f,
    float spacingZ = 1.0f);

}} // namespace mplus::cuda

#endif /* MPLUS_DISTANCE_TRANSFORM_CUH */
