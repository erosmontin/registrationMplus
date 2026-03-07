/*=========================================================================
 *  DerivativeOps.cuh  –  Pure C++ interface to GPU derivative operations.
 *
 *  - Normalize (L2)
 *  - Rescale to [-1, 1]
 *  - Weighted combination of sub-metric derivatives
 *=========================================================================*/
#ifndef MPLUS_DERIVATIVE_OPS_CUH
#define MPLUS_DERIVATIVE_OPS_CUH

#include <cstddef>

namespace mplus { namespace cuda {

/**
 * In-place L2 normalisation of a derivative vector on GPU.
 *
 * @param data  Host pointer to the derivative (double array).
 * @param n     Number of elements.
 */
void gpu_normalize_derivative(double* data, size_t n);

/**
 * In-place rescale to [-1, 1] range on GPU.
 *
 * @param data  Host pointer to the derivative.
 * @param n     Number of elements.
 */
void gpu_rescale_derivative(double* data, size_t n);

/**
 * Weighted combination of multiple derivative arrays.
 *
 * out[i] = sum_k( weights[k] * derivs[k][i] )
 *
 * @param derivs   Array of pointers to derivative arrays (host).
 * @param weights  Weight per derivative array.
 * @param nArrays  Number of derivative arrays.
 * @param n        Elements per array.
 * @param out      Output combined derivative (host pointer).
 */
void gpu_combine_derivatives(
    const double* const* derivs,
    const double*        weights,
    int                  nArrays,
    size_t               n,
    double*              out);

}} // namespace mplus::cuda

#endif /* MPLUS_DERIVATIVE_OPS_CUH */
