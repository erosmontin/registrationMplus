"""GPU-specific tests — skipped when no CUDA device is available.

Run with:  pytest python/tests/test_cuda.py -v -m gpu
"""

import pytest
import numpy as np

# Detect CUDA availability
_CUDA_AVAILABLE = False
try:
    from mplus import _core  # type: ignore[attr-defined]
    _CUDA_AVAILABLE = getattr(_core, "cuda_available", lambda: False)()
except ImportError:
    pass

pytestmark = pytest.mark.gpu

skip_no_cuda = pytest.mark.skipif(
    not _CUDA_AVAILABLE,
    reason="CUDA device or _core extension not available",
)


@skip_no_cuda
class TestCUDADistanceTransform:
    """Verify GPU distance transform matches CPU reference."""

    def test_binary_sphere(self):
        """Signed distance of a centred sphere should be ~radially symmetric."""
        size = 64
        img = np.zeros((size, size, size), dtype=np.float32)
        c = size // 2
        r = 15
        for z in range(size):
            for y in range(size):
                for x in range(size):
                    if (z - c) ** 2 + (y - c) ** 2 + (x - c) ** 2 <= r ** 2:
                        img[z, y, x] = 1.0

        gpu_dist = _core.gpu_distance_transform(img)
        assert gpu_dist.shape == img.shape
        # Centre should be negative (inside)
        assert gpu_dist[c, c, c] < 0
        # Corner should be positive (outside)
        assert gpu_dist[0, 0, 0] > 0

    def test_matches_cpu(self):
        """GPU distance transform should match CPU within tolerance."""
        size = 32
        np.random.seed(42)
        binary = (np.random.rand(size, size, size) > 0.7).astype(np.float32)

        cpu_dist = _core.cpu_distance_transform(binary)
        gpu_dist = _core.gpu_distance_transform(binary)

        np.testing.assert_allclose(gpu_dist, cpu_dist, atol=1e-5, rtol=1e-5)


@skip_no_cuda
class TestCUDADerivatives:
    """GPU derivative operations should match CPU reference."""

    def test_normalize(self):
        np.random.seed(123)
        deriv = np.random.randn(10000).astype(np.float64)
        deriv_gpu = deriv.copy()

        # CPU normalize
        norm = np.sqrt(np.sum(deriv ** 2))
        cpu_result = deriv / norm

        # GPU normalize
        _core.gpu_normalize_derivative(deriv_gpu)

        np.testing.assert_allclose(deriv_gpu, cpu_result, atol=1e-5)

    def test_rescale(self):
        np.random.seed(456)
        deriv = np.random.randn(10000).astype(np.float64)
        deriv_gpu = deriv.copy()

        # CPU rescale to [-1, 1]
        mn, mx = deriv.min(), deriv.max()
        cpu_result = 2.0 * (deriv - mn) / (mx - mn) - 1.0

        # GPU rescale
        _core.gpu_rescale_derivative(deriv_gpu)

        np.testing.assert_allclose(deriv_gpu, cpu_result, atol=1e-5)


@skip_no_cuda
class TestCUDALabelMetric:
    """Label metric GPU kernel should produce correct Dice / kappa values."""

    def test_perfect_overlap(self):
        """Identical label maps → Dice = 1.0 for all labels."""
        size = 32
        labels = np.zeros((size, size, size), dtype=np.int16)
        labels[:16, :, :] = 1
        labels[16:, :, :] = 2

        dice = _core.gpu_label_dice(labels, labels)
        for label_val, score in dice.items():
            assert score == pytest.approx(1.0, abs=1e-5)

    def test_no_overlap(self):
        """Completely disjoint labels → Dice ≈ 0."""
        size = 32
        fixed_labels = np.zeros((size, size, size), dtype=np.int16)
        moving_labels = np.zeros((size, size, size), dtype=np.int16)
        fixed_labels[:16, :, :] = 1
        moving_labels[16:, :, :] = 1

        dice = _core.gpu_label_dice(fixed_labels, moving_labels)
        assert dice.get(1, 0.0) < 0.1


@skip_no_cuda
class TestCUDASpeedup:
    """Benchmark GPU vs CPU — not strict assertions, just sanity checks."""

    def test_distance_transform_speedup(self):
        """GPU should be at least faster than CPU on 128^3."""
        import time

        size = 128
        binary = (np.random.rand(size, size, size) > 0.5).astype(np.float32)

        t0 = time.perf_counter()
        _core.cpu_distance_transform(binary)
        cpu_time = time.perf_counter() - t0

        t0 = time.perf_counter()
        _core.gpu_distance_transform(binary)
        gpu_time = time.perf_counter() - t0

        print(f"\nDT speedup: {cpu_time / gpu_time:.1f}x "
              f"(CPU={cpu_time:.3f}s, GPU={gpu_time:.3f}s)")
        # At least shouldn't be slower
        assert gpu_time < cpu_time * 2
