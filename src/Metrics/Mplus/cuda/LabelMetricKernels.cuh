/*=========================================================================
 *  LabelMetricKernels.cuh  –  Pure C++ interface to GPU label-metric
 *                              (kappa value, derivative, Dice tracking).
 *=========================================================================*/
#ifndef MPLUS_LABEL_METRIC_KERNELS_CUH
#define MPLUS_LABEL_METRIC_KERNELS_CUH

#include <cstddef>
#include <map>
#include <vector>

namespace mplus { namespace cuda {

/**
 * Compute per-label kappa value (MCE of signed distance maps) on GPU.
 *
 * @param fixedDist      Per-label fixed distance maps, each of size nVoxels.
 * @param movingDist     Per-label moving distance maps.
 * @param labelValues    Unique label values.
 * @param labelWeights   Weight per label (parallel to labelValues).
 * @param nLabels        Number of labels.
 * @param nVoxels        Number of voxels per map.
 * @return               Weighted kappa value.
 */
double gpu_kappa_value(
    const float* const* fixedDist,
    const float* const* movingDist,
    const short*        labelValues,
    const double*       labelWeights,
    int                 nLabels,
    size_t              nVoxels);

/**
 * Compute kappa derivative on GPU.
 *
 * @param fixedDist        Per-label fixed distance maps.
 * @param movingDist       Per-label moving distance maps.
 * @param movingDistGrad   Per-label moving distance gradients (3 * nVoxels).
 * @param transformJacobian Sparse Jacobian data for B-spline.
 * @param derivative       Output derivative array (nParameters).
 * @param nLabels          Number of labels.
 * @param nVoxels          Voxels per map.
 * @param nParameters      Number of transform parameters.
 */
void gpu_kappa_derivative(
    const float* const* fixedDist,
    const float* const* movingDist,
    const float* const* movingDistGrad,
    const double*       transformJacobian,
    double*             derivative,
    int                 nLabels,
    size_t              nVoxels,
    size_t              nParameters);

/**
 * Compute per-label Dice coefficients from label maps on GPU.
 *
 * @param fixedLabels   Fixed label map (nVoxels int16 values).
 * @param movingLabels  Moving label map.
 * @param nVoxels       Number of voxels.
 * @param labelValues   Unique label integers to compute Dice for.
 * @param nLabels       Number of labels.
 * @param diceOut       Output: Dice per label (nLabels doubles).
 */
void gpu_label_dice(
    const short* fixedLabels,
    const short* movingLabels,
    size_t       nVoxels,
    const short* labelValues,
    int          nLabels,
    double*      diceOut);

}} // namespace mplus::cuda

#endif /* MPLUS_LABEL_METRIC_KERNELS_CUH */
