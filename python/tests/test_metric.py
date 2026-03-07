"""Tests for MetricWeights, RegistrationConfig and serialization."""

import json
import tempfile
import os
import pytest
import numpy as np

from mplus.config import (
    MetricWeights,
    RegistrationConfig,
    mionly_config,
    mi_label_config,
    multimodal_config,
    label_driven_config,
)


# ── MetricWeights ─────────────────────────────────────────────────────────

class TestMetricWeights:
    def test_defaults(self):
        w = MetricWeights()
        assert w.mi == 1.0
        assert w.ngf == 0.0
        assert w.mse == 0.0
        assert w.nc == 0.0
        assert w.label == 0.0
        assert w.mi_deriv is None

    def test_to_dict_omits_none(self):
        w = MetricWeights(mi=1.0, ngf=0.5)
        d = w.to_dict()
        # None values should be excluded
        assert "mi_deriv" not in d
        assert d["mi"] == 1.0
        assert d["ngf"] == 0.5

    def test_roundtrip(self):
        w = MetricWeights(mi=0.8, ngf=0.3, mse=0.1, mi_deriv=0.9)
        d = w.to_dict()
        w2 = MetricWeights.from_dict(d)
        assert w2.mi == 0.8
        assert w2.ngf == 0.3
        assert w2.mi_deriv == 0.9

    def test_custom_values(self):
        w = MetricWeights(mi=2.0, ngf=1.5, mse=0.5, nc=0.3, label=0.8)
        assert w.mi == 2.0
        assert w.label == 0.8


# ── RegistrationConfig ────────────────────────────────────────────────────

class TestRegistrationConfig:
    def test_defaults(self):
        cfg = RegistrationConfig()
        assert cfg.grid_spacing == 40.0
        assert cfg.num_levels == 3
        assert cfg.iterations_per_level == 200
        assert cfg.derivative_mode == 0
        assert cfg.auto_estimate_eta is True
        assert cfg.use_cuda is False
        assert cfg.label_weights is None

    def test_to_dict(self):
        cfg = RegistrationConfig(grid_spacing=30.0, num_levels=4)
        d = cfg.to_dict()
        assert d["grid_spacing"] == 30.0
        assert d["num_levels"] == 4
        assert "weights" in d
        assert "label_weights" not in d  # None stripped

    def test_roundtrip(self):
        cfg = RegistrationConfig(
            weights=MetricWeights(mi=1.0, ngf=0.5),
            grid_spacing=25.0,
            num_levels=5,
            iterations_per_level=300,
            use_cuda=True,
        )
        d = cfg.to_dict()
        cfg2 = RegistrationConfig.from_dict(d)
        assert cfg2.grid_spacing == 25.0
        assert cfg2.num_levels == 5
        assert cfg2.weights.ngf == 0.5
        assert cfg2.use_cuda is True

    def test_save_load(self, tmp_path):
        cfg = RegistrationConfig(
            weights=MetricWeights(mi=0.8, label=0.3),
            grid_spacing=35.0,
            label_weights={1: 0.5, 2: 0.8},
        )
        path = str(tmp_path / "cfg.json")
        cfg.save(path)

        cfg2 = RegistrationConfig.load(path)
        assert cfg2.weights.mi == 0.8
        assert cfg2.weights.label == 0.3
        assert cfg2.grid_spacing == 35.0

    def test_json_format(self):
        cfg = RegistrationConfig()
        d = cfg.to_dict()
        # Must be JSON-serializable
        s = json.dumps(d)
        assert isinstance(s, str)


# ── Preset configs ────────────────────────────────────────────────────────

class TestPresetConfigs:
    def test_mionly(self):
        cfg = mionly_config()
        assert cfg.weights.mi == 1.0
        assert cfg.weights.ngf == 0.0

    def test_mi_label(self):
        cfg = mi_label_config()
        assert cfg.weights.label == 0.5

    def test_multimodal(self):
        cfg = multimodal_config()
        assert cfg.weights.ngf == 0.5
        assert cfg.weights.mse == 0.2

    def test_label_driven(self):
        cfg = label_driven_config()
        assert cfg.weights.label == 1.0
        assert cfg.weights.mi == 0.3
