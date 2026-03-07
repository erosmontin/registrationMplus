# mPlus v2 Roadmap

## Status Update

This section tracks what has been implemented so far and the recommended next steps.

### Implemented

- Modernised CMake across the repo (`cmake_minimum_required(VERSION 3.18)`, `CXX_STANDARD 17`), added `USE_CUDA`, `BUILD_PYTHON`, `BUILD_TESTS` options.
- Added `src/Metrics/Mplus/itkMplusCompat.h` to handle ITK 4.x / 5.x API differences.
- Added CUDA directory `src/Metrics/Mplus/cuda/` with initial kernel implementations and utilities (DerivativeOps, LabelMetric kernels, DistanceTransform scaffolding, DeviceMemoryPool). Some kernels (derivative ops and Dice forward pass) implemented; others need completion.
- Created Python package scaffold under `python/` with high-level `mplus.registration` API, tests, and a `python/scripts/chain_registration.py` helper for staged pipelines.
- CI workflows and a Dockerfile updated to support v2 packaging; Sphinx docs scaffold added.

### Next steps (short-term priorities)

1. Complete `DistanceTransform` CUDA kernel (finish 1-D passes) and add unit tests verifying numerical equivalence with CPU path.
2. Implement kappa-derivative CUDA kernel and finalise derivative assembly on GPU.
3. Wire CUDA runtime calls into `src/Metrics/Mplus/itkMplus.hxx` behind `#ifdef USE_CUDA` guards so the metric hot-paths can use GPU implementations when available.
4. Complete pybind11 native bindings (`mplus/_core.cpp`) to expose the full pipeline to Python and add packaging steps to produce wheels for common platforms.
5. Add GPU benchmark suite and integrate GPU tests into CI (requires a GPU-enabled runner or cloud GPUs).

If you want I can start on any of these items — say which one to prioritise and I'll implement it next.
# Mplus v2 Modernization Roadmap

This roadmap outlines the modernization of the Mplus registration metric library with three major initiatives:
1. **CUDA C++ Acceleration** – GPU kernels for distance transforms, label metrics, and derivatives
2. **Python Package** – High-level Python bindings via pybind11 with NumPy integration
3. **ITK 5.x Migration** – Modernize to ISO C++17 and latest ITK API

---

## Phase 1: ITK 5.x Migration (2–3 weeks)

### Objectives
- Update project to build against ITK 5.x
- Modernize C++ codebase to C++17 standard
- Maintain backward compatibility via compatibility header
- All existing tests pass

### Key Changes

#### CMake Updates
```cmake
# CMakeLists.txt
cmake_minimum_required(VERSION 3.18)  # Required for proper ITK 5 + CUDA support
find_package(ITK 5.0 REQUIRED)
set(CXX_STANDARD 17)
set(CXX_STANDARD_REQUIRED ON)
```

#### API Migrations

| ITK 4 (current) | ITK 5.x | Notes |
|---|---|---|
| `itkBSplineDeformableTransform<>` | `itkBSplineTransform<>` | Template parameters unchanged |
| `#include "itkMacro.h"` macros | Modern equivalents | Threading → TBB backend |
| `itkMultiThreader` (explicit manage) | `itk::MultiThreaderBase` | Auto-managed thread pool |
| Typedefs with `typedef` | Using aliases | `using Pointer = SmartPointer<Self>;` |
| `itkTypeMacro(A, B)` | `itkOverrideGetNameOfClassMacro(A)` | Name changes in 5.x |
| Virtual destructors (manual) | `~Class() override = default;` | C++17 explicit syntax |

#### Compatibility Header
Create `src/Metrics/Mplus/itkMplusCompat.h`:
- Provides `BSplineTransformCompat<>`
- Maps `ThreaderType` appropriately
- Conditional `#include` for version-specific headers

### Deliverables
- [ ] CMakeLists.txt updated to ITK 5.x
- [ ] `itkMplusCompat.h` created and tested
- [ ] All `.hxx` implementations use modern C++17 syntax
- [ ] All unit tests pass on ITK 5.x

---

## Phase 2: CUDA C++ Acceleration (3–4 weeks)

### Objectives
- Implement GPU kernels for compute-intensive operations
- Maintain CPU fallback with OpenMP
- Achievable speedups: 5–100× depending on kernel
- Compile-time toggling via `USE_CUDA` CMake option

### Architecture

#### Directory Structure
```
src/Metrics/Mplus/
├── itkMplus.h                      ← interface (unchanged)
├── itkMplus.hxx                    ← implementation with #ifdef USE_CUDA
├── CMakeLists.txt
└── cuda/
    ├── CMakeLists.txt              ← CUDA library config
    ├── MplusCudaConfig.h.in         ← config template
    ├── DistanceTransform.cuh        ← header with wrapped C++ interface
    ├── DistanceTransform.cu         ← kernel implementations
    ├── LabelMetricKernels.cuh
    ├── LabelMetricKernels.cu        ← GPU kappa value + derivative
    ├── DerivativeOps.cuh
    ├── DerivativeOps.cu             ← normalize, rescale, combine
    └── DeviceMemoryPool.cuh         ← memory management utilities
```

#### Priority Kernels (by expected impact)

| Kernel | Computational Cost | Expected Speedup | Impl. Difficulty |
|--------|-------------------|-----------------|-----------------|
| **Signed Distance Transform** | Per-label init × 2 | 10–50× | Medium |
| **Kappa Value/Derivative** | Per-voxel distance interpolation | 20–100× | High |
| **Derivative Normalize** | Full vector L2 norm + divide | 5–10× | Low |
| **Derivative Rescale** | Min/max tracking + linear remap | 5–10× | Low |
| **MI/NGF Derivative Assembly** | Jacobian × gradient products | 10–30× | High |

### Key Implementation Details

#### Signed Distance Transform (Felzenszwalb)
- Separable 1D distance transform in X, Y, Z passes
- O(n) time complexity per axis
- Handles anisotropic spacing via proper weighting
- Output: signed distance field (negative inside, positive outside)

**Key files:**
- `src/Metrics/Mplus/cuda/DistanceTransform.cu`
- `src/Metrics/Mplus/cuda/DistanceTransform.cuh`

#### Label Metric GPU Kernels
- Load fixed & moving distance maps to GPU
- Per-voxel distance difference MCE with reduction
- Kappa (weighted label MSE) computation
- Gradient interpolation for derivative evaluation
- Per-label Dice coefficient tracking

**Key files:**
- `src/Metrics/Mplus/cuda/LabelMetricKernels.cu`
- `src/Metrics/Mplus/cuda/LabelMetricKernels.cuh`

#### Derivative Operations
- `gpu_normalize_derivative()` – L2 norm via Thrust + kernel
- `gpu_rescale_derivative()` – Min/max with atomic ops + linear transform
- `gpu_combine_derivatives()` – Weighted sum of MI, NGF, MSE, label derivatives

**Key files:**
- `src/Metrics/Mplus/cuda/DerivativeOps.cu`
- `src/Metrics/Mplus/cuda/DerivativeOps.cuh`

### Usage in C++ Code

All GPU functions wrapped in pure C++ interfaces (no CUDA in header):

```cpp
// In itkMplus.hxx
#ifdef USE_CUDA
  gpu_normalize_derivative(derivative.data_block(), derivative.size());
#else
  // OpenMP fallback (existing code)
  #pragma omp parallel for reduction(+:norm)
  for (unsigned int i = 0; i < derivative.size(); ++i) {
      norm += derivative[i] * derivative[i];
  }
  norm = std::sqrt(norm);
  #pragma omp parallel for
  for (unsigned int i = 0; i < derivative.size(); ++i) {
      derivative[i] /= norm;
  }
#endif
```

### CMake Configuration

```cmake
# Root CMakeLists.txt
cmake_minimum_required(VERSION 3.18)
project(RegistrationSuite LANGUAGES CXX)

option(USE_CUDA "Enable CUDA acceleration" OFF)

if(USE_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    
    # Support modern GPU architectures
    set(CMAKE_CUDA_ARCHITECTURES "60;70;75;80;86;89;90")
    
    # Separable compilation for linking with C++ code
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -fPIC")
    
    add_subdirectory(src/Metrics/Mplus/cuda)
endif()

add_subdirectory(src)
# ...rest of config...
```

### Testing Strategy
- **Benchmark suite**: Compare CPU vs GPU on 3D medical images (256³, 512³)
- **Accuracy tests**: Verify GPU derivative matches CPU within tolerance (1e-5)
- **Memory stress**: Test on large datasets; verify no OOM on realistic data
- **CI/CD**: GitHub Actions with CUDA runner (e.g., `nvidia/cuda:11.8-runtime-ubuntu22.04`)

### Deliverables
- [ ] Distance transform kernel + wrappers
- [ ] Label metric value + derivative kernels
- [ ] Derivative ops (normalize, rescale, combine)
- [ ] Device memory pool utility
- [ ] All kernels wrapped with pure C++ interfaces
- [ ] `#ifdef USE_CUDA` guards in `.hxx`
- [ ] CMake `USE_CUDA` option functional
- [ ] Benchmark/accuracy test suite
- [ ] GPU CI/CD integration

---

## Phase 3: Python Bindings & Package (2–3 weeks)

### Objectives
- Create `mplus-registration` PyPI package
- High-level Pythonic API for registration workflows
- Direct NumPy ↔ ITK image integration
- Optional CUDA acceleration transparent to users
- scikit-build for easy compilation on all platforms

### Architecture

#### Directory Structure
```
python/
├── CMakeLists.txt
├── setup.py                         ← backwards compatibility
├── pyproject.toml                   ← PEP 517 build spec
├── mplus/
│   ├── __init__.py
│   ├── _core.cpp                    ← pybind11 bindings
│   ├── registration.py              ← high-level API
│   └── visualization.py             ← plotting helpers
├── tests/
│   ├── test_metric.py               ← MplusMetric class tests
│   ├── test_registration.py         ← Registration workflow tests
│   └── test_cuda.py                 ← GPU-specific tests
├── examples/
│   ├── basic_registration.py
│   ├── multi_metric_example.py
│   └── label_guided_registration.py
└── docs/
    ├── conf.py                      ← Sphinx config
    └── source/
        ├── index.rst
        ├── api.rst
        └── tutorials/
```

### Core Components

#### 1. pybind11 Bindings (`mplus/_core.cpp`)

**MplusMetric class exposure:**
```cpp
py::class_<MplusType>(m, "MplusMetric")
    .def_property("lambda_val", &MplusType::GetLambda, &MplusType::SetLambda)
    .def_property("alpha", &MplusType::GetAlpha, &MplusType::SetAlpha)
    .def_property("nu", &MplusType::GetNu, &MplusType::SetNu)
    .def_property("yota", &MplusType::GetYota, &MplusType::SetYota)
    .def_property("label_kappa", &MplusType::GetLabelKappa, &MplusType::SetLabelKappa)
    // ... (and derivative weights, eta params, etc.)
    .def("initialize", &MplusType::Initialize)
    .def("print_info", &MplusType::print);
```

**NumPy ↔ ITK conversion:**
- `numpy_to_itk<Dim>()` – Handles shape, spacing, origin
- `itk_to_numpy<Dim>()` – Reverse conversion for output

**High-level function:**
```cpp
py::def("register_images", 
    [](py::array_t<float> fixed, py::array_t<float> moving,
       py::array_t<double> spacing, py::array_t<double> origin,
       py::dict params) -> py::array_t<double> {
        // Convert inputs, run registration, return parameters
    });
```

#### 2. High-Level Registration API (`mplus/registration.py`)

**MetricWeights dataclass:**
```python
@dataclass
class MetricWeights:
    mi: float = 1.0         # Mutual Information
    ngf: float = 0.0        # Normalized Gradient Field
    mse: float = 0.0        # Mean Squared Error
    nc: float = 0.0         # Normalized Correlation
    label: float = 0.0      # Label/Kappa metric
    
    # Optional per-metric derivative weights
    mi_deriv: Optional[float] = None
    ngf_deriv: Optional[float] = None
    mse_deriv: Optional[float] = None
    label_deriv: Optional[float] = None
```

**RegistrationConfig dataclass:**
```python
@dataclass
class RegistrationConfig:
    weights: MetricWeights = field(default_factory=MetricWeights)
    grid_spacing: float = 40.0          # mm, B-spline control points
    num_levels: int = 3                 # Multi-resolution
    iterations_per_level: int = 200
    derivative_mode: int = 0            # 0=weighted, 1=normalized, 2=adaptive
    fixed_eta: float = 5.0              # NGF noise parameter
    moving_eta: float = 5.0
    auto_estimate_eta: bool = True
    num_samples: int = 50000
    use_cuda: bool = False              # Transparent to user
    label_weights: Optional[Dict[int, float]] = None
```

**Registration class:**
```python
class Registration:
    """Multi-metric B-spline deformable registration.
    
    Example
    -------
    >>> config = RegistrationConfig(
    ...     weights=MetricWeights(mi=1.0, ngf=0.5, label=0.3),
    ...     grid_spacing=30.0, num_levels=4
    ... )
    >>> reg = Registration(config)
    >>> result = reg.run(fixed_img, moving_img,
    ...                  fixed_labels=fixed_seg,
    ...                  moving_labels=moving_seg)
    >>> warped = result.warped_image
    >>> print(result.dice_scores)
    """
    
    def __init__(self, config: Optional[RegistrationConfig] = None): ...
    def run(self, fixed, moving, spacing=None, origin=None,
            fixed_labels=None, moving_labels=None) -> RegistrationResult: ...
```

**RegistrationResult dataclass:**
```python
@dataclass
class RegistrationResult:
    transform_parameters: np.ndarray     # B-spline coefficients
    warped_image: np.ndarray
    warped_labels: Optional[np.ndarray] = None
    dice_scores: Optional[Dict[int, float]] = None
    metric_values: Optional[Dict[str, float]] = None
    execution_time_sec: float = 0.0
```

#### 3. Visualization Helpers (`mplus/visualization.py`)

```python
def plot_registration_checkerboard(fixed, moving_warped, slice_idx=None):
    """Display checkerboard overlay of fixed and warped moving."""
    ...

def plot_label_dice(dice_scores: Dict[int, float]):
    """Bar chart of per-label Dice coefficients."""
    ...

def plot_metric_convergence(metric_history):
    """Line plot of metric values across optimizer iterations."""
    ...
```

#### 4. Package Metadata (`pyproject.toml`)

```toml
[build-system]
requires = ["scikit-build-core[pyproject]>=0.5", "pybind11>=2.11"]
build-backend = "scikit_build_core.build"

[project]
name = "mplus-registration"
version = "2.0.0"
description = "Multi-metric image registration with CUDA support"
readme = "README.md"
requires-python = ">=3.9"
authors = [{name = "Registration Suite Contributors"}]
dependencies = [
    "numpy>=1.21",
    "itk>=5.2",
]

[project.optional-dependencies]
visualization = ["matplotlib", "plotly"]
dev = ["pytest>=6.0", "pytest-cov"]
docs = ["sphinx", "sphinx-rtd-theme"]

[tool.scikit-build]
cmake.args = ["-DBUILD_PYTHON=ON"]
cmake.version = ">=3.18"
wheel.packages = ["mplus"]
```

### Installation Methods

#### From PyPI (after release)
```bash
pip install mplus-registration
# With GPU support (requires CUDA toolkit)
pip install mplus-registration[cuda]
```

#### From Source (development)
```bash
git clone https://github.com/yourorg/registrationSuite.git
cd registrationSuite/python
pip install -e .

# With CUDA
CMAKE_ARGS="-DUSE_CUDA=ON" pip install -e .
```

#### Docker (reproducible environment)
```dockerfile
FROM nvidia/cuda:11.8-runtime-ubuntu22.04
RUN apt-get update && apt-get install -y python3-pip
RUN pip install mplus-registration[cuda]
```

### Example Usage

```python
import numpy as np
from mplus import Registration, RegistrationConfig, MetricWeights

# Load images (e.g., using SimpleITK or nibabel)
fixed_img = np.random.randn(256, 256, 256).astype(np.float32)
moving_img = np.random.randn(256, 256, 256).astype(np.float32)
fixed_labels = np.random.randint(0, 10, (256, 256, 256), dtype=np.int16)
moving_labels = np.random.randint(0, 10, (256, 256, 256), dtype=np.int16)

# Configure registration
config = RegistrationConfig(
    weights=MetricWeights(
        mi=1.0,              # Primary metric
        ngf=0.5,             # Gradient field alignment
        label=0.3            # Label overlap (Dice)
    ),
    grid_spacing=35.0,       # mm
    num_levels=4,
    iterations_per_level=250,
    use_cuda=True            # Automatic GPU acceleration
)

# Run registration
reg = Registration(config)
result = reg.run(
    fixed_img, moving_img,
    spacing=np.array([1.0, 1.0, 1.0]),  # mm per voxel
    origin=np.array([0.0, 0.0, 0.0]),
    fixed_labels=fixed_labels,
    moving_labels=moving_labels
)

# Inspect results
print(f"Warped image shape: {result.warped_image.shape}")
print(f"Dice scores: {result.dice_scores}")
print(f"Execution time: {result.execution_time_sec:.2f} sec")

# Visualization
from mplus.visualization import plot_registration_checkerboard, plot_label_dice
plot_registration_checkerboard(fixed_img, result.warped_image, slice_idx=128)
plot_label_dice(result.dice_scores)
```

### Testing Strategy

**Unit tests** (`python/tests/test_metric.py`):
- MplusMetric property getters/setters
- Numerical derivative checking (finite differences)
- GPU ↔ CPU agreement (if CUDA enabled)

**Integration tests** (`python/tests/test_registration.py`):
- End-to-end registration on small synthetic images
- Label overlap computation correctness
- Result shape/type validation

**GPU tests** (`python/tests/test_cuda.py`, skip if no GPU):
- Distance transform vs CPU reference
- Label metric GPU speedup measurement
- Memory limits on large images

Sample pytest invocation:
```bash
pytest python/tests -v --cov=mplus
pytest python/tests/test_cuda.py -v -m "gpu"  # GPU tests only
```

### Documentation

**Sphinx docs** (`python/docs/`):
- API reference (auto-generated from docstrings)
- Tutorial: Basic registration workflow
- Tutorial: Multi-metric configuration
- Tutorial: Label-guided registration
- Performance benchmarks
- Troubleshooting (GPU memory, numerical stability, etc.)

Build docs:
```bash
cd python/docs && make html
```

### Deliverables
- [ ] pybind11 core bindings (`_core.cpp`)
- [ ] NumPy ↔ ITK conversion utilities
- [ ] `MetricWeights` and `RegistrationConfig` dataclasses
- [ ] `Registration` high-level class
- [ ] `RegistrationResult` dataclass
- [ ] Visualization module
- [ ] `pyproject.toml` and `setup.py`
- [ ] CMake Python build integration
- [ ] Unit + integration + GPU tests
- [ ] Sphinx documentation
- [ ] Example notebooks (Jupyter)
- [ ] PyPI release (test environment first)

---

## Phase 4: Polish & Release (1 week)

### CI/CD Pipeline

#### GitHub Actions Workflows

**Build & test on CPU (Linux, macOS, Windows):**
- Runs on every push to `main` and PR
- Tests Python 3.9, 3.10, 3.11, 3.12
- Builds wheels for release

**GPU testing (NVIDIA runner):**
- Runs on tag or manual trigger
- CUDA 11.8 + compute capabilities 60–90
- Benchmarks + accuracy validation

**Documentation:**
- Build Sphinx docs on every commit
- Deploy to GitHub Pages on release

#### Example `.github/workflows/test.yml`:
```yaml
name: Tests

on: [push, pull_request]

jobs:
  cpu-tests:
    runs-on: ${{ matrix.os }}
    strategy:
      matrix:
        os: [ubuntu-latest, macos-latest, windows-latest]
        python-version: ['3.9', '3.10', '3.11']
    steps:
      - uses: actions/checkout@v3
      - uses: actions/setup-python@v4
        with:
          python-version: ${{ matrix.python-version }}
      - run: pip install -e ".[dev]"
      - run: pytest python/tests -v

  gpu-tests:
    runs-on: [self-hosted, gpu]
    if: github.event_name == 'push'
    steps:
      - uses: actions/checkout@v3
      - run: pip install ".[dev,cuda]"
      - run: pytest python/tests/test_cuda.py -v -m gpu
```

### Versioning & Release

**Semantic Versioning:**
- `2.0.0` – First major release (CUDA + Python bindings)
- `2.0.1` – Patch releases for bug fixes
- `2.1.0` – Minor releases for new features
- `3.0.0` – Future breaking changes

**Release process:**
1. Update version in `python/pyproject.toml`, `src/CMakeLists.txt`
2. Tag commit: `git tag v2.0.0`
3. Push tag: `git push origin v2.0.0`
4. GitHub Actions builds wheels + pushes to PyPI
5. Create Release notes on GitHub

### Documentation

**README.md updates:**
- Installation instructions (pip, conda, from source)
- Quick-start example
- Feature matrix (ITK5, CUDA, etc.)
- Citation info

**API Documentation:**
- Sphinx auto-doc from docstrings
- Tutorial notebooks (Jupyter, nbsphinx)
- Gallery of registration examples

**Benchmarks:**
- GPU vs CPU runtime tables (weak scaling)
- Memory usage profiles
- Multi-GPU scaling (if supported)

### Final Deliverables
- [ ] All three phases completed and merged
- [ ] GitHub Actions CI/CD fully functional
- [ ] PyPI package released
- [ ] Complete Sphinx documentation
- [ ] Example Jupyter notebooks
- [ ] Release notes + migration guide (for ITK4 → ITK5)
- [ ] Contributor guide (CONTRIBUTING.md)

---

## High-Level Timeline

| Phase | Duration | Key Milestones |
|-------|----------|-----------------|
| **1: ITK5 Migration** | 2–3 weeks | CMake updated, all tests green, C++17 syntax |
| **2: CUDA Kernels** | 3–4 weeks | Distance transform, label metric, derivative ops; GPU benchmarks |
| **3: Python Package** | 2–3 weeks | pybind11 bindings, Registration class, PyPI ready |
| **4: Polish & Release** | 1 week | CI/CD, documentation, PyPI publish |
| **Total** | ~9–11 weeks | Production-ready v2.0 |

---

## Success Criteria

### Phase 1
- ✓ Builds cleanly against ITK 5.x (Linux, macOS, Windows)
- ✓ All unit tests pass
- ✓ No functional changes to registration output

### Phase 2
- ✓ Signed distance transform GPU kernel correct (tested vs reference)
- ✓ Label metric GPU speedup ≥ 10× on typical 3D images
- ✓ Derivatives numerically match CPU within 1e-5 relative tolerance
- ✓ Zero functional changes to registration results

### Phase 3
- ✓ Python package installs cleanly via pip
- ✓ NumPy arrays convert correctly to/from ITK
- ✓ Registration API is intuitive (docstring examples run)
- ✓ Test suite covers ≥ 80% code

### Phase 4
- ✓ v2.0.0 released on PyPI (test + production)
- ✓ GitHub Pages documentation live
- ✓ CI/CD passes on all platforms
- ✓ Example notebooks work end-to-end

---

## Risk Mitigation

| Risk | Mitigation |
|------|-----------|
| ITK 5.x API changes break existing code | Build compatibility header early; use in-house fork tests |
| GPU kernels unstable on edge devices | Support architectures 60+; extensive CI/CD testing |
| pybind11 NumPy conversion inefficient | Benchmark; use `pybind11/stl.h` + DLPack if needed |
| Documentation lag | Auto-doc; notebook examples as tests (nbval) |
| Dependency hell (numpy, itk versions) | Specify version ranges carefully; use `python_requires` |

---

## References & Resources

- **ITK Migration Guide:** https://itk.org/migrating-to-itk-5/
- **pybind11 Docs:** https://pybind11.readthedocs.io/
- **CUDA C++ Programming:** https://docs.nvidia.com/cuda/cuda-c-programming-guide/
- **scikit-build-core:** https://scikit-build-core.readthedocs.io/
- **Sphinx Documentation:** https://www.sphinx-doc.org/

---

**Last Updated:** March 6, 2026  
**Status:** Planning Phase  
**Owner:** Registration Suite Development Team
