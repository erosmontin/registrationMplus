"""Tests for mplus.visualization module."""

import pytest
import numpy as np

from mplus.visualization import (
    plot_registration_checkerboard,
    plot_label_dice,
    plot_metric_convergence,
    plot_before_after,
    plot_difference_map,
)

# Use non-interactive backend so tests don't open windows
try:
    import matplotlib
    matplotlib.use("Agg")
    _HAS_MPL = True
except ImportError:
    _HAS_MPL = False

skip_no_mpl = pytest.mark.skipif(not _HAS_MPL, reason="matplotlib not installed")


def _dummy_images(size=32):
    fixed = np.random.randn(size, size, size).astype(np.float32)
    moving = np.roll(fixed, 5, axis=0)
    return fixed, moving


@skip_no_mpl
class TestCheckerboard:
    def test_returns_figure(self):
        f, m = _dummy_images()
        fig = plot_registration_checkerboard(f, m)
        assert fig is not None
        assert len(fig.axes) == 3

    def test_custom_slice(self):
        f, m = _dummy_images()
        fig = plot_registration_checkerboard(f, m, slice_idx=10, axis=0)
        assert fig is not None


@skip_no_mpl
class TestLabelDice:
    def test_basic(self):
        scores = {1: 0.85, 2: 0.92, 3: 0.78}
        fig = plot_label_dice(scores)
        assert fig is not None
        assert len(fig.axes) == 1


@skip_no_mpl
class TestConvergence:
    def test_single_metric(self):
        history = [1.0, 0.8, 0.6, 0.4, 0.3, 0.25]
        fig = plot_metric_convergence(history)
        assert fig is not None

    def test_multi_metric(self):
        history = {
            "MI": [1.0, 0.8, 0.6],
            "NGF": [0.5, 0.4, 0.3],
        }
        fig = plot_metric_convergence(history)
        assert fig is not None


@skip_no_mpl
class TestBeforeAfter:
    def test_returns_figure(self):
        f, m = _dummy_images()
        warped = f  # fake "perfect" warp
        fig = plot_before_after(f, m, warped)
        assert fig is not None
        assert len(fig.axes) == 3


@skip_no_mpl
class TestDifferenceMap:
    def test_returns_figure(self):
        f, m = _dummy_images()
        fig = plot_difference_map(f, m)
        assert fig is not None
