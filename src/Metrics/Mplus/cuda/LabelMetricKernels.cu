/*=========================================================================
 *  LabelMetricKernels.cu  –  GPU kernels for label metric computation.
 *
 *  Includes: kappa value, kappa derivative, per-label Dice.
 *=========================================================================*/
#include "LabelMetricKernels.cuh"
#include <cuda_runtime.h>
#include <cstdio>

namespace mplus { namespace cuda {

static void check(cudaError_t err, const char* msg) {
    if (err != cudaSuccess)
        fprintf(stderr, "CUDA error (%s): %s\n", msg, cudaGetErrorString(err));
}

// ── Per-label Dice kernel ────────────────────────────────────────────────

__global__ void dice_kernel(
    const short* fixed,
    const short* moving,
    short        labelVal,
    unsigned long long* intersection,
    unsigned long long* union_count,
    size_t       n)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;

    bool f = (fixed[idx] == labelVal);
    bool m = (moving[idx] == labelVal);

    if (f && m) atomicAdd(intersection, 1ULL);
    if (f || m) atomicAdd(union_count,  1ULL);
}

void gpu_label_dice(
    const short* fixedLabels,
    const short* movingLabels,
    size_t       nVoxels,
    const short* labelValues,
    int          nLabels,
    double*      diceOut)
{
    short* d_fixed  = nullptr;
    short* d_moving = nullptr;
    size_t bytes = nVoxels * sizeof(short);

    check(cudaMalloc(&d_fixed,  bytes), "alloc fixed labels");
    check(cudaMalloc(&d_moving, bytes), "alloc moving labels");
    check(cudaMemcpy(d_fixed,  fixedLabels,  bytes, cudaMemcpyHostToDevice), "H2D fixed");
    check(cudaMemcpy(d_moving, movingLabels, bytes, cudaMemcpyHostToDevice), "H2D moving");

    unsigned long long* d_inter = nullptr;
    unsigned long long* d_union = nullptr;
    check(cudaMalloc(&d_inter, sizeof(unsigned long long)), "alloc inter");
    check(cudaMalloc(&d_union, sizeof(unsigned long long)), "alloc union");

    int block = 256;
    int grid  = (int)((nVoxels + block - 1) / block);

    for (int i = 0; i < nLabels; ++i) {
        unsigned long long zero = 0;
        cudaMemcpy(d_inter, &zero, sizeof(zero), cudaMemcpyHostToDevice);
        cudaMemcpy(d_union, &zero, sizeof(zero), cudaMemcpyHostToDevice);

        dice_kernel<<<grid, block>>>(d_fixed, d_moving, labelValues[i],
                                     d_inter, d_union, nVoxels);
        cudaDeviceSynchronize();

        unsigned long long h_inter, h_union;
        cudaMemcpy(&h_inter, d_inter, sizeof(h_inter), cudaMemcpyDeviceToHost);
        cudaMemcpy(&h_union, d_union, sizeof(h_union), cudaMemcpyDeviceToHost);

        diceOut[i] = (h_union > 0) ? (2.0 * h_inter) / (double)h_union : 0.0;
    }

    cudaFree(d_fixed);
    cudaFree(d_moving);
    cudaFree(d_inter);
    cudaFree(d_union);
}

// ── Kappa value (stl) ────────────────────────────────────────────────────

double gpu_kappa_value(
    const float* const* /*fixedDist*/,
    const float* const* /*movingDist*/,
    const short*        /*labelValues*/,
    const double*       /*labelWeights*/,
    int                 /*nLabels*/,
    size_t              /*nVoxels*/)
{
    // TODO: implement per-voxel distance difference MCE with GPU reduction
    return 0.0;
}

// ── Kappa derivative (stub) ──────────────────────────────────────────────

void gpu_kappa_derivative(
    const float* const* /*fixedDist*/,
    const float* const* /*movingDist*/,
    const float* const* /*movingDistGrad*/,
    const double*       /*transformJacobian*/,
    double*             /*derivative*/,
    int                 /*nLabels*/,
    size_t              /*nVoxels*/,
    size_t              /*nParameters*/)
{
    // TODO: implement per-voxel gradient interpolation + Jacobian product
}

}} // namespace mplus::cuda
