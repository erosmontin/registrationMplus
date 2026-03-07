"""mplus-registration: Multi-metric image registration with optional CUDA support.

High-level Python API for the Mplus combined-metric registration framework.

Quick Start
-----------
>>> from mplus import Registration, RegistrationConfig, MetricWeights
>>> config = RegistrationConfig(weights=MetricWeights(mi=1.0, ngf=0.5))
>>> reg = Registration(config)
>>> result = reg.run(fixed_image, moving_image)
"""

__version__ = "2.0.0"

from .config import MetricWeights, RegistrationConfig
from .registration import Registration, RegistrationResult

__all__ = [
    "MetricWeights",
    "RegistrationConfig",
    "Registration",
    "RegistrationResult",
    "__version__",
]
