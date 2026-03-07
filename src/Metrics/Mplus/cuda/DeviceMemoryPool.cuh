/*=========================================================================
 *  DeviceMemoryPool.cuh  –  Simple CUDA device memory pool.
 *
 *  Pre-allocates a block of GPU memory and hands out sub-allocations
 *  to avoid repeated cudaMalloc/cudaFree overhead during iterative
 *  registration.
 *=========================================================================*/
#ifndef MPLUS_DEVICE_MEMORY_POOL_CUH
#define MPLUS_DEVICE_MEMORY_POOL_CUH

#include <cuda_runtime.h>
#include <cstddef>
#include <cstdio>

namespace mplus { namespace cuda {

class DeviceMemoryPool {
public:
    /**
     * Create a pool of the given size in bytes.
     * @param poolSizeBytes  Total pre-allocated GPU memory.
     */
    explicit DeviceMemoryPool(size_t poolSizeBytes)
        : m_poolSize(poolSizeBytes), m_offset(0), m_base(nullptr)
    {
        cudaError_t err = cudaMalloc(&m_base, poolSizeBytes);
        if (err != cudaSuccess) {
            fprintf(stderr, "DeviceMemoryPool: cudaMalloc failed (%s)\n",
                    cudaGetErrorString(err));
            m_base = nullptr;
            m_poolSize = 0;
        }
    }

    ~DeviceMemoryPool() {
        if (m_base) cudaFree(m_base);
    }

    // Non-copyable
    DeviceMemoryPool(const DeviceMemoryPool&) = delete;
    DeviceMemoryPool& operator=(const DeviceMemoryPool&) = delete;

    /**
     * Allocate *bytes* from the pool.  Returns nullptr if exhausted.
     * Allocations are 256-byte aligned.
     */
    void* allocate(size_t bytes) {
        // Align to 256 bytes
        size_t aligned = (bytes + 255) & ~size_t(255);
        if (m_offset + aligned > m_poolSize) return nullptr;
        void* ptr = static_cast<char*>(m_base) + m_offset;
        m_offset += aligned;
        return ptr;
    }

    /** Reset the pool (invalidates all previous pointers). */
    void reset() { m_offset = 0; }

    /** Bytes remaining in pool. */
    size_t available() const { return m_poolSize - m_offset; }

    /** Total pool size. */
    size_t capacity() const { return m_poolSize; }

private:
    size_t m_poolSize;
    size_t m_offset;
    void*  m_base;
};

}} // namespace mplus::cuda

#endif /* MPLUS_DEVICE_MEMORY_POOL_CUH */
