/*=========================================================================
 *  DerivativeOps.cu  –  GPU kernels for derivative normalise, rescale,
 *                        and weighted combination.
 *=========================================================================*/
#include "DerivativeOps.cuh"
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/extrema.h>
#include <cmath>
#include <cstdio>

namespace mplus { namespace cuda {

static void chk(cudaError_t e, const char* m) {
    if (e != cudaSuccess)
        fprintf(stderr, "CUDA (%s): %s\n", m, cudaGetErrorString(e));
}

// File-scope functor for squared norm (must not be local to a function)
struct SquareOp {
    __device__ double operator()(double x) const { return x * x; }
};

// ── Normalize ────────────────────────────────────────────────────────────

__global__ void div_kernel(double* data, size_t n, double norm) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) data[i] /= norm;
}

void gpu_normalize_derivative(double* data, size_t n) {
    double* d_data = nullptr;
    size_t bytes = n * sizeof(double);
    chk(cudaMalloc(&d_data, bytes), "alloc norm");
    chk(cudaMemcpy(d_data, data, bytes, cudaMemcpyHostToDevice), "H2D norm");

    // Compute L2 norm via Thrust
    thrust::device_ptr<double> dp(d_data);
    double sq_sum = thrust::transform_reduce(
        dp, dp + n,
        SquareOp{},
        0.0,
        thrust::plus<double>());
    double norm = std::sqrt(sq_sum);

    if (norm > 1e-10) {
        int block = 256;
        int grid  = (int)((n + block - 1) / block);
        div_kernel<<<grid, block>>>(d_data, n, norm);
        cudaDeviceSynchronize();
    }

    chk(cudaMemcpy(data, d_data, bytes, cudaMemcpyDeviceToHost), "D2H norm");
    cudaFree(d_data);
}

// ── Rescale ──────────────────────────────────────────────────────────────

__global__ void rescale_kernel(double* data, size_t n, double minVal, double range) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        data[i] = 2.0 * (data[i] - minVal) / range - 1.0;
    }
}

void gpu_rescale_derivative(double* data, size_t n) {
    double* d_data = nullptr;
    size_t bytes = n * sizeof(double);
    chk(cudaMalloc(&d_data, bytes), "alloc rescale");
    chk(cudaMemcpy(d_data, data, bytes, cudaMemcpyHostToDevice), "H2D rescale");

    thrust::device_ptr<double> dp(d_data);
    auto minmax = thrust::minmax_element(dp, dp + n);
    double minVal = *minmax.first;
    double maxVal = *minmax.second;
    double range  = maxVal - minVal;

    if (range > 1e-15) {
        int block = 256;
        int grid  = (int)((n + block - 1) / block);
        rescale_kernel<<<grid, block>>>(d_data, n, minVal, range);
        cudaDeviceSynchronize();
    }

    chk(cudaMemcpy(data, d_data, bytes, cudaMemcpyDeviceToHost), "D2H rescale");
    cudaFree(d_data);
}

// ── Combine ──────────────────────────────────────────────────────────────

__global__ void combine_kernel(
    const double* const* __restrict__ derivs,
    const double* __restrict__        weights,
    int nArrays, size_t n,
    double* __restrict__ out)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    double sum = 0.0;
    for (int k = 0; k < nArrays; ++k) {
        sum += weights[k] * derivs[k][i];
    }
    out[i] = sum;
}

void gpu_combine_derivatives(
    const double* const* derivs,
    const double*        weights,
    int                  nArrays,
    size_t               n,
    double*              out)
{
    size_t bytes = n * sizeof(double);

    // Allocate device arrays
    double** d_derivs = nullptr;
    double*  d_weights = nullptr;
    double*  d_out = nullptr;

    chk(cudaMalloc(&d_derivs, nArrays * sizeof(double*)), "alloc ptrs");
    chk(cudaMalloc(&d_weights, nArrays * sizeof(double)), "alloc wts");
    chk(cudaMalloc(&d_out, bytes), "alloc out");

    // Copy each derivative array to device
    double** h_dptrs = new double*[nArrays];
    for (int k = 0; k < nArrays; ++k) {
        chk(cudaMalloc(&h_dptrs[k], bytes), "alloc deriv k");
        chk(cudaMemcpy(h_dptrs[k], derivs[k], bytes, cudaMemcpyHostToDevice), "H2D dk");
    }
    chk(cudaMemcpy(d_derivs, h_dptrs, nArrays * sizeof(double*), cudaMemcpyHostToDevice), "H2D ptrs");
    chk(cudaMemcpy(d_weights, weights, nArrays * sizeof(double), cudaMemcpyHostToDevice), "H2D wts");

    int block = 256;
    int grid  = (int)((n + block - 1) / block);
    combine_kernel<<<grid, block>>>((const double* const*)d_derivs, d_weights, nArrays, n, d_out);
    cudaDeviceSynchronize();

    chk(cudaMemcpy(out, d_out, bytes, cudaMemcpyDeviceToHost), "D2H out");

    // Cleanup
    for (int k = 0; k < nArrays; ++k) cudaFree(h_dptrs[k]);
    delete[] h_dptrs;
    cudaFree(d_derivs);
    cudaFree(d_weights);
    cudaFree(d_out);
}

}} // namespace mplus::cuda
