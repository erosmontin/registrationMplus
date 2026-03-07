"""Tests for the Registration class and RegistrationResult."""

import pytest
import numpy as np

from mplus.config import MetricWeights, RegistrationConfig
from mplus.registration import Registration, RegistrationResult


# ── RegistrationResult ────────────────────────────────────────────────────

class TestRegistrationResult:
    def test_defaults(self):
        r = RegistrationResult()
        assert r.warped_image.size == 0
        assert r.warped_labels is None
        assert r.dice_scores is None
        assert r.execution_time_sec == 0.0
        assert r.convergence_history is None

    def test_with_data(self):
        warped = np.random.randn(64, 64, 64).astype(np.float32)
        r = RegistrationResult(
            warped_image=warped,
            dice_scores={1: 0.85, 2: 0.92},
            metric_values={"MI": -1.2, "NGF": 0.8},
            execution_time_sec=12.5,
        )
        assert r.warped_image.shape == (64, 64, 64)
        assert r.dice_scores[1] == 0.85
        assert r.metric_values["MI"] == -1.2


# ── Registration init & validation ────────────────────────────────────────

class TestRegistrationInit:
    def test_default_config(self):
        reg = Registration()
        assert reg.config.grid_spacing == 40.0

    def test_custom_config(self):
        cfg = RegistrationConfig(grid_spacing=25.0, num_levels=5)
        reg = Registration(cfg)
        assert reg.config.grid_spacing == 25.0
        assert reg.config.num_levels == 5

    def test_invalid_grid_spacing(self):
        with pytest.raises(ValueError, match="grid_spacing"):
            Registration(RegistrationConfig(grid_spacing=-1.0))

    def test_invalid_num_levels(self):
        with pytest.raises(ValueError, match="num_levels"):
            Registration(RegistrationConfig(num_levels=0))

    def test_invalid_iterations(self):
        with pytest.raises(ValueError, match="iterations_per_level"):
            Registration(RegistrationConfig(iterations_per_level=0))


# ── Input validation ──────────────────────────────────────────────────────

class TestRegistrationInputValidation:
    def setup_method(self):
        self.reg = Registration()

    def test_rejects_2d_images(self):
        f = np.zeros((64, 64), dtype=np.float32)
        m = np.zeros((64, 64), dtype=np.float32)
        with pytest.raises(ValueError, match="3-D"):
            self.reg.run(f, m)

    def test_rejects_shape_mismatch(self):
        f = np.zeros((64, 64, 64), dtype=np.float32)
        m = np.zeros((64, 64, 32), dtype=np.float32)
        with pytest.raises(ValueError, match="Shape mismatch"):
            self.reg.run(f, m)

    def test_rejects_partial_labels(self):
        f = np.zeros((32, 32, 32), dtype=np.float32)
        m = np.zeros((32, 32, 32), dtype=np.float32)
        fl = np.zeros((32, 32, 32), dtype=np.int16)
        with pytest.raises(ValueError, match="Both"):
            self.reg.run(f, m, fixed_labels=fl)


# ── dry_run ───────────────────────────────────────────────────────────────

class TestDryRun:
    def test_basic(self):
        reg = Registration(RegistrationConfig(
            weights=MetricWeights(mi=1.0, ngf=0.5),
            grid_spacing=32.0,
        ))
        info = reg.dry_run((128, 128, 128))
        assert "n_parameters" in info
        assert "grid_size" in info
        assert "active_metrics" in info
        assert "estimated_memory_mb" in info
        assert info["n_parameters"] > 0
        assert "MI" in info["active_metrics"]
        assert "NGF" in info["active_metrics"]

    def test_label_memory_higher(self):
        reg_no_label = Registration(RegistrationConfig(
            weights=MetricWeights(mi=1.0),
            grid_spacing=32.0,
        ))
        reg_label = Registration(RegistrationConfig(
            weights=MetricWeights(mi=1.0, label=0.5),
            grid_spacing=32.0,
        ))
        info1 = reg_no_label.dry_run((128, 128, 128))
        info2 = reg_label.dry_run((128, 128, 128))
        assert info2["estimated_memory_mb"] > info1["estimated_memory_mb"]

    def test_grid_size_scales(self):
        reg = Registration(RegistrationConfig(grid_spacing=64.0))
        small = reg.dry_run((128, 128, 128))
        large = reg.dry_run((256, 256, 256))
        assert large["n_parameters"] > small["n_parameters"]
