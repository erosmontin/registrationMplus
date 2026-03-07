"""Configuration management for Mplus registration.

Supports loading/saving registration configs from JSON files.
"""

import json
from dataclasses import dataclass, field, asdict, is_dataclass
from pathlib import Path
from typing import Optional, Dict, Any


@dataclass
class MetricWeights:
    """Weights for the Mplus combined metric."""
    mi: float = 1.0              # Mutual Information (lambda)
    ngf: float = 0.0             # Normalized Gradient Field (alpha)
    mse: float = 0.0             # Mean Squared Error (nu)
    nc: float = 0.0              # Normalized Correlation (yota)
    label: float = 0.0           # Label/Kappa metric

    # Optional per-metric derivative weights (None = same as value weight)
    mi_deriv: Optional[float] = None
    ngf_deriv: Optional[float] = None
    mse_deriv: Optional[float] = None
    label_deriv: Optional[float] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary, omitting None values."""
        return {k: v for k, v in asdict(self).items() if v is not None}

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "MetricWeights":
        """Create from dictionary."""
        return cls(**data)


@dataclass
class RegistrationConfig:
    """Configuration for B-spline registration.

    Attributes
    ----------
    weights : MetricWeights
        Combined metric weights (MI, NGF, MSE, NC, label)
    grid_spacing : float
        B-spline control point spacing in mm (default 40.0)
    num_levels : int
        Number of multi-resolution levels (default 3)
    iterations_per_level : int
        Gradient descent iterations per level (default 200)
    derivative_mode : int
        How to combine metric derivatives:
        - 0: consistent weighted sum (default, safe for LBFGS-B)
        - 1: normalize + rescale (RSGD only)
        - 2: main-metric adaptive scaling (LBFGS-B safe)
    fixed_eta : float
        NGF noise parameter for fixed image (default 5.0)
    moving_eta : float
        NGF noise parameter for moving image (default 5.0)
    auto_estimate_eta : bool
        Automatically estimate eta from image gradients (default True)
    num_samples : int
        Number of points sampled per iteration (default 50000)
    use_cuda : bool
        Enable CUDA acceleration if available (default False)
    label_weights : Dict[int, float], optional
        Per-label override weights [0, 1]. Labels not in map get 1/nLabels.
    """
    weights: MetricWeights = field(default_factory=MetricWeights)
    grid_spacing: float = 40.0
    num_levels: int = 3
    iterations_per_level: int = 200
    derivative_mode: int = 0
    fixed_eta: float = 5.0
    moving_eta: float = 5.0
    auto_estimate_eta: bool = True
    num_samples: int = 50000
    use_cuda: bool = False
    label_weights: Optional[Dict[int, float]] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        d = asdict(self)
        # Convert weights to nested dict
        d["weights"] = self.weights.to_dict()
        # Remove None label_weights
        if d["label_weights"] is None:
            del d["label_weights"]
        return d

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "RegistrationConfig":
        """Create from dictionary."""
        data = data.copy()
        # Extract and create weights
        weights_data = data.pop("weights", {})
        weights = MetricWeights.from_dict(weights_data)
        return cls(weights=weights, **data)

    def save(self, path: str) -> None:
        """Save configuration to JSON file.

        Parameters
        ----------
        path : str
            Path to output JSON file
        """
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def load(cls, path: str) -> "RegistrationConfig":
        """Load configuration from JSON file.

        Parameters
        ----------
        path : str
            Path to JSON configuration file

        Returns
        -------
        RegistrationConfig
            Loaded configuration
        """
        with open(path, "r") as f:
            data = json.load(f)
        return cls.from_dict(data)


# Example configurations (can be serialized to JSON)

def mionly_config() -> RegistrationConfig:
    """Mutual Information only (traditional registration)."""
    return RegistrationConfig(
        weights=MetricWeights(mi=1.0),
        grid_spacing=40.0,
        num_levels=3,
    )


def mi_label_config() -> RegistrationConfig:
    """MI-driven with label overlap constraint."""
    return RegistrationConfig(
        weights=MetricWeights(mi=1.0, label=0.5),
        grid_spacing=35.0,
        num_levels=4,
    )


def multimodal_config() -> RegistrationConfig:
    """Multi-metric for challenging alignment (MI + NGF + MSE)."""
    return RegistrationConfig(
        weights=MetricWeights(mi=1.0, ngf=0.5, mse=0.2),
        grid_spacing=30.0,
        num_levels=4,
        iterations_per_level=300,
    )


def label_driven_config() -> RegistrationConfig:
    """Primary focus: label overlap; secondary: intensity."""
    return RegistrationConfig(
        weights=MetricWeights(mi=0.3, label=1.0),
        grid_spacing=25.0,
        num_levels=5,
        iterations_per_level=250,
    )
