# mplus-registration

Multi-metric B-spline image registration with optional CUDA acceleration.

## Installation

```bash
pip install mplus-registration
```

### From source (development)

```bash
git clone https://github.com/yourorg/registrationSuite.git
cd registrationSuite/python
pip install -e ".[dev]"
```

### With CUDA support

```bash
CMAKE_ARGS="-DUSE_CUDA=ON" pip install -e .
```

## Quick Start

```python
import numpy as np
from mplus import Registration, RegistrationConfig, MetricWeights

# Configure
config = RegistrationConfig(
    weights=MetricWeights(mi=1.0, ngf=0.5, label=0.3),
    grid_spacing=30.0,
    num_levels=4,
)

# Run
reg = Registration(config)
result = reg.run(fixed_image, moving_image,
                 spacing=np.array([1.0, 1.0, 1.0]))

print(f"Execution time: {result.execution_time_sec:.1f}s")
```

## Features

- **Multi-metric**: Combine MI, NGF, MSE, NC, and label-overlap metrics
- **B-spline deformable registration**: Control point spacing in mm
- **Label-guided**: Dice-driven registration with per-label weights
- **CUDA acceleration**: Optional GPU kernels for 10-100x speedup
- **JSON configs**: Save/load registration parameters
- **Visualization**: Checkerboard, Dice charts, convergence plots

## License

BSD-3-Clause
