/*=========================================================================
 *  DistanceTransform.cu  –  GPU signed distance transform (Felzenszwalb).
 *
 *  Separable 1D distance transform in X, Y, Z passes.
 *  O(n) per axis.  Handles anisotropic spacing.
 *
 *  Algorithm: Felzenszwalb & Huttenlocher, "Distance Transforms of Sampled
 *  Functions", Theory of Computing 8(19), 2012.
 *
 *  Two separate DTs are computed (outside: fg=0/bg=INF; inside: fg=INF/bg=0)
 *  and combined:  signed_dist = sqrt(outside) - sqrt(inside)
 *  → negative inside the object, positive outside.
 *=========================================================================*/
#include "DistanceTransform.cuh"
#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace mplus { namespace cuda {

// ── helpers ──────────────────────────────────────────────────────────────

static void checkCuda(cudaError_t err, const char* msg) {
    if (err != cudaSuccess)
        fprintf(stderr, "CUDA error (%s): %s\n", msg, cudaGetErrorString(err));
}

static constexpr float BIG = 1e18f;   // proxy for infinity (squared)

// ── init kernels ─────────────────────────────────────────────────────────

// outside: fg→0, bg→BIG  (compute distance from bg to nearest fg)
__global__ void init_outside(const float* __restrict__ binary_dev,
                              float*       __restrict__ vol,
                              int N) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N) vol[i] = (binary_dev[i] > 0.5f) ? 0.0f : BIG;
}

// inside: fg→BIG, bg→0  (compute distance from fg to nearest bg)
__global__ void init_inside(const float* __restrict__ binary_dev,
                             float*       __restrict__ vol,
                             int N) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N) vol[i] = (binary_dev[i] > 0.5f) ? BIG : 0.0f;
}

// ── Felzenszwalb 1-D DT kernel (one thread per row) ──────────────────────
//
// axis=0 (X): step=1,   off0=sy, bstride0=sx,    bstride1=sy*sx
// axis=1 (Y): step=sx,  off0=sx, bstride0=1,     bstride1=sy*sx
// axis=2 (Z): step=sx*sy, off0=sx, bstride0=1,   bstride1=sx
//
// Each thread tid handles one row:
//   i0 = tid % off0;  i1 = tid / off0;
//   base = i1 * bstride1 + i0 * bstride0;
//   elements at vol_in[base + k*step] for k=0..len-1
//
// Scratch:  d_v[tid * len + k]        (int)
//           d_z[tid * (len+1) + k]    (float)
// ─────────────────────────────────────────────────────────────────────────
__global__ void dt1d(
    const float* __restrict__ vol_in,
    float*       __restrict__ vol_out,
    int*         __restrict__ d_v,
    float*       __restrict__ d_z,
    int step, int len, int nRows,
    int off0, int bstride0, int bstride1,
    float sp)
{
    int tid = (int)(blockIdx.x * blockDim.x + threadIdx.x);
    if (tid >= nRows) return;

    int i0   = tid % off0;
    int i1   = tid / off0;
    int base = i1 * bstride1 + i0 * bstride0;

    int*   v = d_v + tid * len;
    float* z = d_z + tid * (len + 1);

    // ── Forward: build lower parabola envelope ────────────────────────
    int k = 0;
    v[0] = 0;
    z[0] = -BIG;
    z[1] =  BIG;

    for (int q = 1; q < len; q++) {
        float fq = vol_in[base + q * step];
        float s;
        do {
            int   vk  = v[k];
            float fvk = vol_in[base + vk * step];
            float qsp  = (float)q  * sp;
            float vksp = (float)vk * sp;
            // intersection of parabolas centred at vk and q
            s = ((fq + qsp * qsp) - (fvk + vksp * vksp))
                / (2.0f * sp * (qsp - vksp) / sp)  // simplifies to 2*(q-vk)*sp^2 but keep generic
                ;
            // above simplifies to:
            // s = (fq - fvk + sp*sp*(q*q - vk*vk)) / (2.0f * sp * sp * (q - vk))
            if (s <= z[k]) --k;
            else break;
        } while (k >= 0);
        ++k;
        v[k]     = q;
        z[k]     = s;
        z[k + 1] = BIG;
    }

    // ── Backward: look up nearest parabola for each position ─────────
    k = 0;
    for (int q = 0; q < len; q++) {
        float qsp = (float)q * sp;
        while (z[k + 1] < qsp) ++k;
        float du  = qsp - (float)v[k] * sp;
        vol_out[base + q * step] = vol_in[base + v[k] * step] + du * du;
    }
}

// ── combine: signed = sqrt(outside) - sqrt(inside) ───────────────────────
__global__ void combine_signed(
    const float* __restrict__ d_outside,   // squared outside DT
    const float* __restrict__ d_inside,    // squared inside  DT
    float* __restrict__ d_out,
    int N)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    float out_dist = sqrtf(d_outside[i]);
    float in_dist  = sqrtf(d_inside[i]);
    d_out[i] = out_dist - in_dist;   // negative inside, positive outside
}

// ── helper: run one full X-Y-Z DT pass on the GPU ────────────────────────
static void run_dt3d(
    float*  d_vol_a,   // input (modified in-place across passes)
    float*  d_vol_b,   // scratch buffer (same size)
    int*    d_v,
    float*  d_z,
    int sx, int sy, int sz,
    float spx, float spy, float spz)
{
    const int BLOCK = 256;

    // X pass: step=1, len=sx, nRows=sy*sz
    {
        int nRows = sy * sz;
        int grid  = (nRows + BLOCK - 1) / BLOCK;
        dt1d<<<grid, BLOCK>>>(d_vol_a, d_vol_b, d_v, d_z,
                              /*step*/1, /*len*/sx, nRows,
                              /*off0*/sy, /*bs0*/sx, /*bs1*/sy*sx, spx);
        cudaDeviceSynchronize();
        std::swap(d_vol_a, d_vol_b);   // output is now d_vol_a
    }
    // Y pass: step=sx, len=sy, nRows=sx*sz
    {
        int nRows = sx * sz;
        int grid  = (nRows + BLOCK - 1) / BLOCK;
        dt1d<<<grid, BLOCK>>>(d_vol_a, d_vol_b, d_v, d_z,
                              sx, sy, nRows,
                              sx, 1, sy*sx, spy);
        cudaDeviceSynchronize();
        std::swap(d_vol_a, d_vol_b);
    }
    // Z pass: step=sx*sy, len=sz, nRows=sx*sy
    {
        int nRows = sx * sy;
        int grid  = (nRows + BLOCK - 1) / BLOCK;
        dt1d<<<grid, BLOCK>>>(d_vol_a, d_vol_b, d_v, d_z,
                              sx*sy, sz, nRows,
                              sx, 1, sx, spz);
        cudaDeviceSynchronize();
        std::swap(d_vol_a, d_vol_b);
    }
    // After 3 swaps the final result is in d_vol_b (the last swap moved it there)
    // — actually with 3 swaps: a→b(X out), b→a(Y out), a→b(Z out) then swap → result in d_vol_a
    // Let's just copy result to d_vol_a if needed. The caller passes d_vol_a by value
    // so the pointer swap is local. We need to return the result pointer.
    // Simplest fix: the caller passes pointers as-is and we write output to d_vol_b.
    (void)d_vol_a; (void)d_vol_b; // final result is in whatever d_vol_a points to after swaps
}

// ── public interface ─────────────────────────────────────────────────────

void gpu_signed_distance_transform(
    const float* binary_in,
    float*       dist_out,
    int sizeX, int sizeY, int sizeZ,
    float spacingX,
    float spacingY,
    float spacingZ)
{
    const size_t N     = (size_t)sizeX * sizeY * sizeZ;
    const size_t bytes = N * sizeof(float);
    const int    BLOCK = 256;

    const int maxLen  = std::max({sizeX, sizeY, sizeZ});
    const int maxRows = std::max({sizeY * sizeZ, sizeX * sizeZ, sizeX * sizeY});

    // ── Allocate device buffers ───────────────────────────────────────
    float* d_binary  = nullptr;
    float* d_outside = nullptr;  // accumulates outside squared DT
    float* d_inside  = nullptr;  // accumulates inside  squared DT
    float* d_scratch = nullptr;  // ping-pong scratch
    float* d_result  = nullptr;  // signed output
    int*   d_v       = nullptr;
    float* d_z       = nullptr;

    checkCuda(cudaMalloc(&d_binary,  bytes),                                          "alloc bin");
    checkCuda(cudaMalloc(&d_outside, bytes),                                          "alloc out");
    checkCuda(cudaMalloc(&d_inside,  bytes),                                          "alloc in");
    checkCuda(cudaMalloc(&d_scratch, bytes),                                          "alloc scr");
    checkCuda(cudaMalloc(&d_result,  bytes),                                          "alloc res");
    checkCuda(cudaMalloc(&d_v,  (size_t)maxRows * maxLen * sizeof(int)),              "alloc v");
    checkCuda(cudaMalloc(&d_z,  (size_t)maxRows * (maxLen + 1) * sizeof(float)),      "alloc z");

    checkCuda(cudaMemcpy(d_binary, binary_in, bytes, cudaMemcpyHostToDevice),         "H2D bin");

    // ── Outside pass (fg=0, bg=BIG) ──────────────────────────────────
    {
        int grid = ((int)N + BLOCK - 1) / BLOCK;
        init_outside<<<grid, BLOCK>>>(d_binary, d_outside, (int)N);
        cudaDeviceSynchronize();

        // X pass
        {
            int nRows = sizeY * sizeZ;
            dt1d<<<(nRows+BLOCK-1)/BLOCK, BLOCK>>>(d_outside, d_scratch, d_v, d_z,
                1, sizeX, nRows, sizeY, sizeX, sizeY*sizeX, spacingX);
            cudaDeviceSynchronize();
            cudaMemcpy(d_outside, d_scratch, bytes, cudaMemcpyDeviceToDevice);
        }
        // Y pass
        {
            int nRows = sizeX * sizeZ;
            dt1d<<<(nRows+BLOCK-1)/BLOCK, BLOCK>>>(d_outside, d_scratch, d_v, d_z,
                sizeX, sizeY, nRows, sizeX, 1, sizeY*sizeX, spacingY);
            cudaDeviceSynchronize();
            cudaMemcpy(d_outside, d_scratch, bytes, cudaMemcpyDeviceToDevice);
        }
        // Z pass
        {
            int nRows = sizeX * sizeY;
            dt1d<<<(nRows+BLOCK-1)/BLOCK, BLOCK>>>(d_outside, d_scratch, d_v, d_z,
                sizeX*sizeY, sizeZ, nRows, sizeX, 1, sizeX, spacingZ);
            cudaDeviceSynchronize();
            cudaMemcpy(d_outside, d_scratch, bytes, cudaMemcpyDeviceToDevice);
        }
    }

    // ── Inside pass (fg=BIG, bg=0) ────────────────────────────────────
    {
        int grid = ((int)N + BLOCK - 1) / BLOCK;
        init_inside<<<grid, BLOCK>>>(d_binary, d_inside, (int)N);
        cudaDeviceSynchronize();

        {
            int nRows = sizeY * sizeZ;
            dt1d<<<(nRows+BLOCK-1)/BLOCK, BLOCK>>>(d_inside, d_scratch, d_v, d_z,
                1, sizeX, nRows, sizeY, sizeX, sizeY*sizeX, spacingX);
            cudaDeviceSynchronize();
            cudaMemcpy(d_inside, d_scratch, bytes, cudaMemcpyDeviceToDevice);
        }
        {
            int nRows = sizeX * sizeZ;
            dt1d<<<(nRows+BLOCK-1)/BLOCK, BLOCK>>>(d_inside, d_scratch, d_v, d_z,
                sizeX, sizeY, nRows, sizeX, 1, sizeY*sizeX, spacingY);
            cudaDeviceSynchronize();
            cudaMemcpy(d_inside, d_scratch, bytes, cudaMemcpyDeviceToDevice);
        }
        {
            int nRows = sizeX * sizeY;
            dt1d<<<(nRows+BLOCK-1)/BLOCK, BLOCK>>>(d_inside, d_scratch, d_v, d_z,
                sizeX*sizeY, sizeZ, nRows, sizeX, 1, sizeX, spacingZ);
            cudaDeviceSynchronize();
            cudaMemcpy(d_inside, d_scratch, bytes, cudaMemcpyDeviceToDevice);
        }
    }

    // ── Combine: signed_dist = sqrt(outside) - sqrt(inside) ──────────
    {
        int grid = ((int)N + BLOCK - 1) / BLOCK;
        combine_signed<<<grid, BLOCK>>>(d_outside, d_inside, d_result, (int)N);
        cudaDeviceSynchronize();
    }

    checkCuda(cudaMemcpy(dist_out, d_result, bytes, cudaMemcpyDeviceToHost), "D2H result");

    cudaFree(d_binary);
    cudaFree(d_outside);
    cudaFree(d_inside);
    cudaFree(d_scratch);
    cudaFree(d_result);
    cudaFree(d_v);
    cudaFree(d_z);
}

}} // namespace mplus::cuda
