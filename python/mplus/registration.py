"""High-level registration API for Mplus combined-metric registration.

This module provides a Pythonic wrapper around the C++ Mplus metric
and ITK registration pipeline.  When the compiled ``_core`` extension
is available it delegates to C++; otherwise a pure-Python stub is
provided for configuration validation and dry-run workflows.
"""

from __future__ import annotations

import time
import subprocess
import tempfile
import os
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

from .config import MetricWeights, RegistrationConfig

# Try to import compiled pybind11 module (built via CMake / scikit-build)
try:
    from . import _core  # type: ignore[attr-defined]
    _HAS_CORE = True
except ImportError:
    _HAS_CORE = False


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------

@dataclass
class RegistrationResult:
    """Container for registration output.

    Attributes
    ----------
    transform_parameters : np.ndarray
        B-spline coefficients (or affine matrix elements).
    warped_image : np.ndarray
        Moving image resampled onto fixed-image grid.
    warped_labels : np.ndarray, optional
        Moving label map resampled (nearest-neighbour) when label
        images are provided.
    dice_scores : dict[int, float], optional
        Per-label Dice coefficients after registration.
    metric_values : dict[str, float], optional
        Final value of each active sub-metric (MI, NGF, MSE …).
    execution_time_sec : float
        Wall-clock time of the registration in seconds.
    convergence_history : list[float], optional
        Combined metric value at each optimizer iteration.
    """

    transform_parameters: np.ndarray = field(default_factory=lambda: np.array([]))
    warped_image: np.ndarray = field(default_factory=lambda: np.array([]))
    warped_labels: Optional[np.ndarray] = None
    dice_scores: Optional[Dict[int, float]] = None
    metric_values: Optional[Dict[str, float]] = None
    execution_time_sec: float = 0.0
    convergence_history: Optional[List[float]] = None


# ---------------------------------------------------------------------------
# Registration class
# ---------------------------------------------------------------------------

class Registration:
    """Multi-metric B-spline deformable image registration.

    Parameters
    ----------
    config : RegistrationConfig, optional
        Registration hyper-parameters.  Defaults to MI-only with
        standard multi-resolution settings.

    Examples
    --------
    >>> from mplus import Registration, RegistrationConfig, MetricWeights
    >>> config = RegistrationConfig(
    ...     weights=MetricWeights(mi=1.0, ngf=0.5, label=0.3),
    ...     grid_spacing=30.0,
    ...     num_levels=4,
    ... )
    >>> reg = Registration(config)
    >>> result = reg.run(fixed_img, moving_img,
    ...                  fixed_labels=fixed_seg,
    ...                  moving_labels=moving_seg)
    >>> print(result.warped_image.shape)
    >>> print(result.dice_scores)
    """

    def __init__(self, config: Optional[RegistrationConfig] = None) -> None:
        self.config = config or RegistrationConfig()
        self._validate_config()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(
        self,
        fixed: np.ndarray,
        moving: np.ndarray,
        spacing: Optional[np.ndarray] = None,
        origin: Optional[np.ndarray] = None,
        fixed_labels: Optional[np.ndarray] = None,
        moving_labels: Optional[np.ndarray] = None,
    ) -> RegistrationResult:
        """Execute registration.

        Parameters
        ----------
        fixed : np.ndarray
            Fixed (reference) image, 3-D float32 array.
        moving : np.ndarray
            Moving (source) image, same shape as *fixed*.
        spacing : np.ndarray, optional
            Voxel size in mm for each axis.  Defaults to isotropic 1 mm.
        origin : np.ndarray, optional
            World-coordinate origin.  Defaults to (0, 0, 0).
        fixed_labels : np.ndarray, optional
            Integer label map aligned with the fixed image.
        moving_labels : np.ndarray, optional
            Integer label map aligned with the moving image.

        Returns
        -------
        RegistrationResult
            Warped image, transform parameters, Dice scores, etc.
        """
        fixed = np.asarray(fixed, dtype=np.float32)
        moving = np.asarray(moving, dtype=np.float32)

        if fixed.ndim != 3 or moving.ndim != 3:
            raise ValueError("Only 3-D images are supported.")
        if fixed.shape != moving.shape:
            raise ValueError(
                f"Shape mismatch: fixed {fixed.shape} vs moving {moving.shape}"
            )

        if spacing is None:
            spacing = np.ones(3, dtype=np.float64)
        if origin is None:
            origin = np.zeros(3, dtype=np.float64)

        spacing = np.asarray(spacing, dtype=np.float64)
        origin = np.asarray(origin, dtype=np.float64)

        # Label validation
        if (fixed_labels is not None) != (moving_labels is not None):
            raise ValueError(
                "Both fixed_labels and moving_labels must be provided, or neither."
            )
        has_labels = fixed_labels is not None

        t0 = time.perf_counter()

        if _HAS_CORE:
            result = self._run_native(
                fixed, moving, spacing, origin,
                fixed_labels, moving_labels,
            )
        else:
            result = self._run_cli_fallback(
                fixed, moving, spacing, origin,
                fixed_labels, moving_labels,
            )

        result.execution_time_sec = time.perf_counter() - t0
        return result

    def dry_run(self, fixed_shape: Tuple[int, ...]) -> Dict[str, object]:
        """Report what *would* happen without running the optimizer.

        Useful for checking parameter counts, memory estimates, etc.

        Parameters
        ----------
        fixed_shape : tuple of int
            Shape of the fixed image (e.g. ``(256, 256, 256)``).

        Returns
        -------
        dict
            Keys: ``n_parameters``, ``grid_size``, ``active_metrics``,
            ``estimated_memory_mb``.
        """
        cfg = self.config
        gs = cfg.grid_spacing
        grid_size = tuple(int(s / gs) + 3 for s in fixed_shape)
        n_params = int(np.prod(grid_size)) * 3  # 3-D displacement field

        active = []
        w = cfg.weights
        if w.mi:    active.append("MI")
        if w.ngf:   active.append("NGF")
        if w.mse:   active.append("MSE")
        if w.nc:    active.append("NC")
        if w.label: active.append("Label")

        # Rough memory: images, distance maps, derivatives
        voxels = int(np.prod(fixed_shape))
        mem_mb = (voxels * 4 * 2  # fixed + moving float32
                  + n_params * 8   # derivative double
                  + (voxels * 4 * 10 if "Label" in active else 0)  # dist maps
                  ) / (1024 ** 2)

        return {
            "n_parameters": n_params,
            "grid_size": grid_size,
            "active_metrics": active,
            "estimated_memory_mb": round(mem_mb, 1),
        }

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    def _validate_config(self) -> None:
        """Sanity-check configuration values."""
        cfg = self.config
        if cfg.grid_spacing <= 0:
            raise ValueError("grid_spacing must be positive")
        if cfg.num_levels < 1:
            raise ValueError("num_levels must be >= 1")
        if cfg.iterations_per_level < 1:
            raise ValueError("iterations_per_level must be >= 1")

    def _run_native(
        self,
        fixed: np.ndarray,
        moving: np.ndarray,
        spacing: np.ndarray,
        origin: np.ndarray,
        fixed_labels: Optional[np.ndarray],
        moving_labels: Optional[np.ndarray],
    ) -> RegistrationResult:
        """Delegate to compiled pybind11 ``_core`` extension."""
        params_dict = self._config_to_core_params()

        raw = _core.register_images(  # type: ignore[union-attr]
            fixed, moving, spacing, origin, params_dict,
            fixed_labels=fixed_labels,
            moving_labels=moving_labels,
        )
        return RegistrationResult(
            transform_parameters=np.asarray(raw.get("transform", [])),
            warped_image=np.asarray(raw.get("warped", np.array([]))),
            warped_labels=raw.get("warped_labels"),
            dice_scores=raw.get("dice"),
            metric_values=raw.get("metrics"),
            convergence_history=raw.get("history"),
        )

    def _run_cli_fallback(
        self,
        fixed: np.ndarray,
        moving: np.ndarray,
        spacing: np.ndarray,
        origin: np.ndarray,
        fixed_labels: Optional[np.ndarray],
        moving_labels: Optional[np.ndarray],
    ) -> RegistrationResult:
        """Fall back to calling the CLI executables via subprocess.

        This allows the Python API to work even without a compiled
        pybind11 extension — it writes NIfTI temp files, invokes the
        C++ binary, and reads back the output.
        """
        try:
            import SimpleITK as sitk
        except ImportError:
            raise RuntimeError(
                "Neither the compiled _core extension nor SimpleITK is "
                "available.  Install mplus with C++ bindings or "
                "'pip install SimpleITK' for CLI fallback."
            )

        cfg = self.config

        with tempfile.TemporaryDirectory(prefix="mplus_") as tmpdir:
            f_path = os.path.join(tmpdir, "fixed.nii.gz")
            m_path = os.path.join(tmpdir, "moving.nii.gz")
            o_path = os.path.join(tmpdir, "warped.nii.gz")
            vf_path = os.path.join(tmpdir, "vf.mhd")

            # Write images via SimpleITK
            f_img = sitk.GetImageFromArray(fixed)
            f_img.SetSpacing(spacing.tolist())
            f_img.SetOrigin(origin.tolist())
            sitk.WriteImage(f_img, f_path)

            m_img = sitk.GetImageFromArray(moving)
            m_img.SetSpacing(spacing.tolist())
            m_img.SetOrigin(origin.tolist())
            sitk.WriteImage(m_img, m_path)

            # Build CLI command
            cmd = self._build_cli_command(f_path, m_path, o_path, vf_path)
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
            if proc.returncode != 0:
                raise RuntimeError(
                    f"Registration CLI failed (rc={proc.returncode}):\n{proc.stderr}"
                )

            # Read warped result
            warped = sitk.GetArrayFromImage(sitk.ReadImage(o_path))

            return RegistrationResult(
                warped_image=warped.astype(np.float32),
            )

    def _build_cli_command(
        self,
        fixed_path: str,
        moving_path: str,
        output_path: str,
        vf_path: str,
    ) -> List[str]:
        """Construct a CLI invocation for 3DRegBsplines."""
        cfg = self.config
        w = cfg.weights

        cmd = [
            "3DRegBsplines",
            "-f", fixed_path,
            "-m", moving_path,
            "-o", output_path,
            "-v", vf_path,
            "-g", str(cfg.grid_spacing),
            "-L", str(cfg.num_levels),
            "-l", str(w.mi),            # lambda (MI)
            "-n", str(w.ngf),           # alpha  (NGF)
            "-N", str(w.mse),           # nu     (MSE)
            "-y", str(w.nc),            # yota   (NC)
            "-p", str(cfg.iterations_per_level),
            f"--numberofthreads={os.cpu_count() or 4}",
        ]

        if cfg.auto_estimate_eta:
            cmd += ["--autoeta"]
        else:
            cmd += ["--fixedeta", str(cfg.fixed_eta),
                    "--movingeta", str(cfg.moving_eta)]

        return cmd

    def _config_to_core_params(self) -> dict:
        """Convert RegistrationConfig → dict for the C++ pybind11 interface."""
        cfg = self.config
        w = cfg.weights
        return {
            "lambda": w.mi,
            "lambda_deriv": w.mi_deriv if w.mi_deriv is not None else w.mi,
            "alpha": w.ngf,
            "alpha_deriv": w.ngf_deriv if w.ngf_deriv is not None else w.ngf,
            "nu": w.mse,
            "nu_deriv": w.mse_deriv if w.mse_deriv is not None else w.mse,
            "yota": w.nc,
            "yota_deriv": w.label_deriv if w.label_deriv is not None else w.nc,
            "label_kappa": w.label,
            "label_kappa_deriv": w.label_deriv if w.label_deriv is not None else w.label,
            "grid_spacing": cfg.grid_spacing,
            "num_levels": cfg.num_levels,
            "iterations": cfg.iterations_per_level,
            "derivative_mode": cfg.derivative_mode,
            "fixed_eta": cfg.fixed_eta,
            "moving_eta": cfg.moving_eta,
            "auto_eta": cfg.auto_estimate_eta,
            "num_samples": cfg.num_samples,
            "use_cuda": cfg.use_cuda,
            "label_weights": cfg.label_weights or {},
        }
