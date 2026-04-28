# RegistrationSuite v1.2 — Documentation

> **Authors:** Dr. Eros Montin Ph.D. (eros.montin@gmail.com) and contributors.  
> **Primary citation:** Montin E., Belfatto A., Bologna M., Meroni S., Cavatorta C., Pecori E., Diletto B., Massimino M., Oprandi M.C., Poggi G., Arrigoni F., Peruzzo D., Pignoli E., Gandola L., Cerveri P., & Mainardi L. (2020). *A multi-metric registration strategy for the alignment of longitudinal brain images in pediatric oncology.* Medical & Biological Engineering & Computing, 58(4), 843–855. https://doi.org/10.1007/s11517-019-02109-4

---

## Part 1 — Introduction

RegistrationSuite is a modular, ITK-based 3-D image registration toolbox designed for medical imaging research. It provides parametric rigid, similarity, affine, and deformable (B-spline) registration as stand-alone command-line executables and as a C++ library with Python bindings. The central design choice is the **Mplus composite metric** (`itk::Mplus`), a custom `ImageToImageMetric` that linearly combines up to seven independent sub-metrics in a single optimization pass, enabling multi-metric and multi-modal registration without running sequential pipelines.

### Intended use cases

| Scenario | Recommended tool | Metric preset |
|---|---|---|
| Monomodal longitudinal alignment (e.g., MRI–MRI) | `3DRegAffine` → `3DRegBsplines` | `--modality singlemodal` |
| Multimodal alignment (e.g., CT–MR, CT–PET) | `3DRegAffine` → `3DRegBsplines` | `--modality multimodal` |
| Rigid/similarity pre-alignment | `3DRegSimilarity` | custom weights |
| Coarse-to-fine parametric alignment | `3DRegAffineMultiLevel` | custom weights |

### Repository layout (v1.2)

```
src/
  3DRegistration/
    Version.h                    # version constants shared by all executables
    3DRegAffine/src/             # 3DRegAffine.cxx, 3DRegAffineMultiLevel.cxx,
    |                            #   3DRegSimilarity.cxx
    3DRegBsplines/src/           # 3DRegBsplines.cxx
  includes/
    RegistrationCommon.h         # shared option printing / utilities
    imageUtils.h                 # image I/O wrappers
    registrationUtils.h          # reusable registration helpers
  Metrics/
    Mplus/                       # itkMplus.h, itkMplus.hxx  <- composite metric
    NGF/NGFImageMetric/          # Normalized Gradient Field implementation
bld2/bin/                        # compiled executables (build tree)
examples/                        # Python usage examples and JSON config
python/mplus/                    # Python bindings
```

---

## Part 2 — Registration Modules and ITK Transforms

### 2.1 Common pipeline

Every registration executable in the suite follows the same four-stage pipeline:

```
1. Image I/O          Read fixed image, moving image; optionally read masks /
                       label maps / an initial transform.
2. Initialization     Centre-of-mass or identity transform initialisation
                       (itkCenteredTransformInitializer).
3. Optimization       itk::ImageRegistrationMethod drives the optimizer over
                       the Mplus composite metric.
4. Output             Write resampled moving image, optional transform file,
                       optional deformation-vector-field image.
```

All tools share the same set of **Mplus metric flags** (alpha, lambda, nu, yota, rho, sigma, kappa — see Part 3) so the metric configuration is fully portable between registration stages.

---

### 2.2 `3DRegAffine` — Affine registration

**Source:** `src/3DRegistration/3DRegAffine/src/3DRegAffine.cxx`  
**Binary:** `bld2/bin/3DRegAffine`

Performs global linear registration using `itk::AffineTransform<double, 3>` optimised with a **Regular Step Gradient Descent (RSGD)** optimizer.

#### Transform

`itk::AffineTransform<double, 3>` encodes a 12-parameter affine map (3×3 matrix + 3D translation). The `--subtype` flag restricts the degrees of freedom:

| `--subtype` | DoF | Description |
|---|---|---|
| `translation` | 3 | Translation only |
| `rotation` | ≤ 6 | Rotation about image centre |
| `scaling` | ≤ 6 | Isotropic or anisotropic scaling |
| `affine` | 12 | Full affine (default) |

#### Optimizer — `itk::RegularStepGradientDescentOptimizer`

| Flag | Default | Description |
|---|---|---|
| `--maxnumberofiterations,-I` | 1000 | Maximum iterations |
| `--minimumsteplength,-S` | 0.1 | Stop when step falls below this value |
| `--maximumsteplength,-X` | 1.0 | Initial step size |
| `--relaxationfactor,-R` | 0.5 | Step reduction factor on direction change |
| `--gradientmagnitudetolerance,-G` | 1e-4 | Gradient convergence threshold |

#### Quick start

```bash
3DRegAffine -f fixed.nii.gz -m moving.nii.gz -o out.nii.gz \
    --modality multimodal \
    -T affine.tfm
```

---

### 2.3 `3DRegSimilarity` — Similarity registration

**Source:** `src/3DRegistration/3DRegAffine/src/3DRegSimilarity.cxx`  
**Binary:** `bld2/bin/3DRegSimilarity`

Same pipeline as `3DRegAffine` but uses `itk::Similarity3DTransform<double>` — a 7-parameter transform (rotation quaternion, uniform scale, 3D translation). Suitable when isotropic scaling correction is needed between acquisitions (e.g., scanner calibration drift).

Same Mplus metric, same RSGD optimizer, and same flag set as `3DRegAffine`.

---

### 2.4 `3DRegBsplines` — Deformable B-spline registration

**Source:** `src/3DRegistration/3DRegBsplines/src/3DRegBsplines.cxx`  
**Binary:** `bld2/bin/3DRegBsplines`

Performs deformable registration with a **cubic B-spline Free-Form Deformation (FFD)** transform optimised by the **L-BFGS-B** bounded-memory quasi-Newton optimizer.

#### Transform — `itk::BSplineTransform<double, 3, 3>`

The deformation field is parameterised as a regular control-point lattice with cubic B-spline basis functions, defined over the fixed-image domain.

| Flag | Default | Description |
|---|---|---|
| `--gridresolution,-g` | 50 | Control-point spacing in mm |
| `--gridposition,-G` | N | Read initial grid from an external file |
| `--transformin,-W` | N | Initialise with a pre-computed rigid/affine transform |
| `--bound` | 0 | Bound type: 0 = unbounded, 1 = lower, 2 = both, 3 = upper |
| `--lbound` / `--ubound` | 0 / 0 | Displacement bounds (mm) when bound ≠ 0 |
| `--overlappadding` | 1 | Control points outside image border (min = spline order = 3) |
| `--meshmarginsize` | 0.0 | Extra domain margin (mm) around the fixed image |

The transform is saved as an `itk::CompositeTransform` (initial affine + B-spline) when both `--transformin` and `--transformout` are used together, allowing downstream tools to apply the full chain with a single `itk::ResampleImageFilter`.

#### Optimizer — `itk::LBFGSBOptimizer`

L-BFGS-B is a memory-limited quasi-Newton method that supports per-variable box constraints, making it well suited to B-spline control-point optimization.

| Flag | Default | Description |
|---|---|---|
| `--maxnumberofiterations,-I` | 1000 | Maximum iterations |
| `--numberofevaluations,-E` | 500 | Maximum function evaluations |
| `--numberofcorrections,-C` | 5 | L-BFGS-B history length (5–20 recommended) |
| `--costfunctionconvergencefactor,-F` | 1e12 | Convergence: 1e12 = low, 1e7 = moderate, 1e1 = very high accuracy |
| `--projectedgradienttolerance,-P` | 1e-5 | Gradient convergence threshold |

#### Modality presets

The `--modality` flag provides ready-made metric weight combinations:

| Preset | Active metrics | Weights |
|---|---|---|
| `multimodal` | MI + NGF | alpha=1.0, lambda=0.5; MSE, NC, GD, NMI = 0 |
| `singlemodal` | MSE + NC | nu=1.0, yota=0.5; MI, NGF, GD, NMI = 0 |
| `custom` (default) | User-defined | All weights follow individual CLI flags |

Any individual weight flag supplied on the command line overrides the preset.

#### Quick start

```bash
# Step 1 — affine pre-alignment
3DRegAffine -f fixed.nii.gz -m moving.nii.gz -o pre.nii.gz \
    --modality multimodal -T pre_affine.tfm

# Step 2 — deformable refinement
3DRegBsplines -f fixed.nii.gz -m moving.nii.gz -o deformable.nii.gz \
    --modality multimodal \
    --transformin pre_affine.tfm \
    --gridresolution 20 \
    --maxnumberofiterations 500 \
    -T bspline.tfm
```

---

### 2.5 `3DRegAffineMultiLevel` — Multi-level affine registration

**Source:** `src/3DRegistration/3DRegAffine/src/3DRegAffineMultiLevel.cxx`  
**Binary:** `bld2/bin/3DRegAffineMultiLevel`

Performs coarse-to-fine affine registration using an image pyramid. Successive pyramid levels use progressively finer resolution. Shares the full Mplus metric and RSGD optimizer; see the `3DRegAffine` parameter table for flags.

---

### 2.6 ITK Transform I/O

All executables use `itk::TransformFileWriter` / `itk::TransformFileReader` to serialise transforms as human-readable text files. The header `ITKFactoryRegistration/itkTransformIOFactoryRegisterManager.h` ensures all ITK image and transform I/O factories are registered at program start, guaranteeing that NIfTI, MetaImage, and ITK transform formats are available without manual factory registration in user code.

Transforms can be chained via `itk::CompositeTransform<double, 3>` and applied with `itk::ResampleImageFilter` using the same physical-space grid as the fixed image.

---

## Part 3 — Metrics

### 3.1 The Mplus composite metric — design rationale

All four registration tools in v1.2 use the custom ITK metric class `itk::Mplus` (`src/Metrics/Mplus/itkMplus.h`), which inherits from `itk::ImageToImageMetric<FixedImageType, MovingImageType>`.

`Mplus` evaluates up to **seven independent sub-metrics** simultaneously and returns a single combined value and derivative:

$$M_\text{total} = \alpha \cdot M_\text{MI} + \lambda \cdot M_\text{NGF} + \nu \cdot M_\text{MSE} + \iota \cdot M_\text{NC} + \rho \cdot M_\text{GD} + \sigma \cdot M_\text{NMI} + \kappa \cdot M_\text{Label}$$

Setting any weight to `0` disables that sub-metric entirely with no computation overhead. The design intention — described in the primary citation — is to allow complementary metrics to guide registration jointly, for example combining a multimodal metric (MI) with a structural metric (NGF) to handle both intensity and gradient alignment simultaneously.

Each sub-metric has its own **independent sampling percentage** so that expensive metrics (e.g., MI) can use fewer samples than cheaper ones (e.g., MSE) within the same iteration.

---

### 3.2 Sub-metrics reference

#### 3.2.1 Mattes Mutual Information (MI) — multimodal

**ITK class:** `itk::MattesMutualInformationImageToImageMetric`  
**Weight flags:** `--alpha` (value), `--alphaderivative` (derivative)

Mutual Information measures the statistical dependence between the intensity distributions of the two images. The Mattes implementation estimates the joint probability density using a Parzen-windowed histogram with stochastic pixel sampling, making it smooth and differentiable.

**When to use:** Multimodal registration (CT–MR, CT–PET, PET–MR). Primary metric for the `multimodal` preset.

| Flag | Default | Description |
|---|---|---|
| `--alpha,-a` | 1.0 | Value weight of MI in the combined cost |
| `--alphaderivative,-A` | 1.0 | Derivative weight (set to 0 to exclude MI from the gradient) |
| `--mattesnumberofbins,-b` | 64 | Number of histogram bins. Typical range: 32–128 |
| `--mattespercentage,-p` | 0.1 | Fraction of fixed-image pixels sampled (0.1 = 10 %) |
| `--explicitPDFderivatives` | false | Explicit (memory-intensive) vs implicit (default) PDF derivatives |
| `--bsplinecaching,-B` | true | Cache B-spline weights across MI evaluations (faster) |

**Recommended values:** `--mattesnumberofbins 64 --mattespercentage 0.1`.  
Increase `--mattespercentage` toward 0.3–0.5 for small images or when convergence is noisy.

**Primary citation:**  
Maes F., Collignon A., Vandermeulen D., Marchal G., Suetens P. (1997). *Multimodality image registration by maximization of mutual information.* IEEE Transactions on Medical Imaging, 16(2), 187–198.

---

#### 3.2.2 Normalized Gradient Fields (NGF) — multimodal / structural

**ITK class:** `itkNormalizedGradientFieldImageToImageMetric` (custom, `src/Metrics/NGF/`)  
**Weight flags:** `--lambda` (value), `--lambdaderivative` (derivative)

NGF compares the orientation of image gradients rather than intensity values. Two images are considered aligned when their normalised gradient vectors are parallel (or anti-parallel), making the metric insensitive to intensity scale and offset. This is particularly powerful for multimodal data where structures that appear bright in one modality appear dark in another.

The **eta** parameter ($\eta$) controls normalisation: gradients with magnitude below $\eta$ are suppressed as noise. It can be set manually or auto-estimated from the image.

**When to use:** Complement to MI in multimodal registration; also useful independently when structural boundaries are the primary registration target.

| Flag | Default | Description |
|---|---|---|
| `--lambda,-l` | 0.0 | Value weight of NGF |
| `--lambdaderivative,-L` | 0.0 | Derivative weight |
| `--etavaluefixed,-r` | -1 | η for fixed image (−1 = auto-estimate) |
| `--etavaluemoving,-s` | -1 | η for moving image (−1 = auto-estimate) |
| `--NGFevaluator` | 0 | Evaluator kernel (see table below) |
| `--ngfprecompute` | false | Pre-compute moving NGF once and resample (faster, approximate) |
| `--ngfspacing` | `4,4,4` | NGF sub-sampling spacing per dimension (x,y,z in mm) |
| `--ngfpercentage` | 0.1 | Fraction of pixels sampled |

**NGF evaluator kernels:**

| Value | Name | Description |
|---|---|---|
| 0 | scalar | $\langle \hat{n}_F, \hat{n}_M \rangle$ — dot product of normalised gradients |
| 1 | cross | $\|\hat{n}_F \times \hat{n}_M\|$ — cross product magnitude (penalises misalignment) |
| 2 | scdelta | Scalar evaluator with delta-based normalisation |
| 3 | Delta | Delta-function based evaluator |
| 4 | Delta2 | Squared delta evaluator |

**Primary citation:**  
Haber E., Modersitzki J. (2006). *Intensity gradient based registration and fusion of multi-modal images.* Methods of Information in Medicine, 45(1), 153–161.

---

#### 3.2.3 Mean Squared Error (MSE) — monomodal

**ITK class:** `itk::MeanSquaresImageToImageMetric`  
**Weight flags:** `--nu` (value), `--nuderivative` (derivative)

Minimises the average squared difference in intensity between the fixed and transformed moving image. The simplest and fastest metric. Requires a linear (or near-linear) intensity relationship between images.

**When to use:** Monomodal registration (MRI–MRI, CT–CT) with consistent intensity scales; primary metric for the `singlemodal` preset.

| Flag | Default | Description |
|---|---|---|
| `--nu,-n` | 0.0 | Value weight |
| `--nuderivative,-N` | 0.0 | Derivative weight |
| `--msepercentage` | 0.1 | Fraction of pixels sampled |

**Tip:** Normalise image intensities to [0, 1] or zero-mean/unit-variance before using MSE to ensure the weight scale is comparable to other sub-metrics.

---

#### 3.2.4 Normalized Correlation (NC) — monomodal

**ITK class:** `itk::NormalizedCorrelationImageToImageMetric`  
**Weight flags:** `--yota` (value), `--yotaderivative` (derivative)

Measures the normalised linear correlation (Pearson's $r$) between intensity values. Unlike MSE, NC is invariant to global intensity offset and scale, making it more robust when scanner gain or bias field introduces a multiplicative/additive relationship between images.

**When to use:** Monomodal registration with global intensity offset or scale difference (e.g., different echo times, mild protocol changes). Combined with MSE in the `singlemodal` preset.

| Flag | Default | Description |
|---|---|---|
| `--yota,-y` | 0.0 | Value weight |
| `--yotaderivative,-Y` | 0.0 | Derivative weight |
| `--ncpercentage` | 0.1 | Fraction of pixels sampled |

---

#### 3.2.5 Gradient Difference (GD) — monomodal / edge-driven

**ITK class:** `itk::GradientDifferenceImageToImageMetric`  
**Weight flags:** `--rho` (value), `--rhoderivative` (derivative)

Computes the difference between image gradients (spatial derivatives) rather than intensities, making it more sensitive to edge positions and less sensitive to smooth intensity variations in homogeneous background regions.

**When to use:** Monomodal cases where boundary accuracy is critical; complements MSE by emphasising structural edges.

| Flag | Default | Description |
|---|---|---|
| `--rho` | 0.0 | Value weight |
| `--rhoderivative` | 0.0 | Derivative weight |
| `--gdpercentage` | 0.1 | Fraction of pixels sampled |

---

#### 3.2.6 Normalized Mutual Information (NMI) — multimodal

**ITK class:** `itk::NormalizedMutualInformationHistogramImageToImageMetric`  
**Weight flags:** `--sigma` (value), `--sigmaderivative` (derivative)

NMI is defined as $\text{NMI} = (H(F) + H(M)) / H(F,M)$, where $H$ denotes Shannon entropy. Normalisation makes NMI less sensitive to the image overlap area, advantageous when the overlap changes significantly during optimization (e.g., large initial misalignments, limited FOV).

**When to use:** Multimodal registration as an alternative or supplement to Mattes MI; particularly stable under varying overlap.

| Flag | Default | Description |
|---|---|---|
| `--sigma` | 0.0 | Value weight |
| `--sigmaderivative` | 0.0 | Derivative weight |
| `--nmibins` | 64 | Number of histogram bins |
| `--nmipercentage` | 0.1 | Fraction of pixels sampled |

**Primary citation:**  
Studholme C., Hill D.L.G., Hawkes D.J. (1999). *An overlap invariant entropy measure of 3D medical image alignment.* Pattern Recognition, 32(1), 71–86.

---

#### 3.2.7 Label Map Distance (Label / kappa) — anatomy-guided

**Implementation:** Custom distance-map MSE inside `itkMplus.hxx`  
**Weight flags:** `--labelkappa` (value), `--labelkappaderiv` (derivative)

When fixed and moving label maps (segmentations) are available, this sub-metric computes the MSE of signed distance maps for each label, penalising misalignment of anatomical structures. Per-label weights allow fine control, downweighting uncertain or irrelevant structures.

To improve numerical stability and make weights portable across datasets, signed distances are clamped and normalized before the loss is applied. Optional narrow-band sampling can focus the metric near boundaries, and an optional Huber loss can reduce outlier influence.

**When to use:** Any modality when anatomical segmentations are available; particularly useful for brain registration where major structures (ventricles, tumour, cortex) can anchor the deformation.

| Flag | Default | Description |
|---|---|---|
| `--fixedlabelmap` | N | Path to fixed-image label map (NIfTI short integer) |
| `--movinglabelmap` | N | Path to moving-image label map |
| `--labelkappa` | 0.0 | Global weight for all labels (0 = disabled) |
| `--labelkappaderiv` | 0.0 | Global derivative weight |
| `--labelkappavec` | "" | Per-label value weights: `L1:w1,L2:w2,...` |
| `--labelkappaderivvec` | "" | Per-label derivative weights |
| `--labelsamples` | 0.1 | Voxel samples per label. Fraction `(0,1]` **or** absolute count `>1` |
| `--labeldistmax` | 20.0 | Clamp signed distances to ±this value (mm) before loss |
| `--labelnarrowband` | false | Use only samples within a distance band from either boundary |
| `--labelbandwidth` | 5.0 | Narrow-band half-width (mm) |
| `--labelhuber` | false | Use Huber loss on normalized residuals |
| `--labelhuberdelta` | 0.25 | Huber threshold in normalized units |
| `--labelreport` | 1 | Report Dice coefficient every N iterations (0 = off) |

Label loss details:

- Signed distances are computed from fixed and moving label maps at the internal working resolution.
- Distances are clamped to `[-labeldistmax,+labeldistmax]`, residual is normalized by `labeldistmax`, then squared loss (default) or Huber loss is applied.
- If `labelkappaderiv=0`, moving distance-map gradients are not precomputed, substantially reducing initialization time.

---

### 3.3 Derivative combination modes

The way per-metric gradients are merged into the single optimizer gradient is controlled by `--derivativemode`:

| Mode | Name | Formula | Compatible optimizers |
|---|---|---|---|
| 0 | Consistent weighted sum | $\nabla M = \sum_i w_i \nabla M_i$ | All (default) |
| 1 | Normalized + rescale | Each $\nabla M_i$ unit-normalised before summing | RSGD only |
| 2 | Main-metric adaptive | $\nabla M_i$ scaled to match gradient magnitude of the main metric | All (L-BFGS-B safe) |

> **Note:** `--derivativemode 1` is incompatible with L-BFGS-B. `3DRegBsplines` will reject it at startup with an explicit error message.

The main metric for mode 2 is selected with `--mainmetric`:  
0 = MI · 1 = NGF · 2 = MSE · 3 = NC · 4 = Label · 5 = GD · 6 = NMI

---

### 3.4 Global metric options

| Flag | Default | Description |
|---|---|---|
| `--metricoverlap` | true | Restrict metric to the spatial overlap between fixed and transformed moving image |
| `--metricpadding` | 0 | Additional voxel padding around the overlap region |
| `--fixedimagethreshold,-t` | −∞ | Mask out fixed-image pixels below this intensity |
| `--numberofthreads` | 2 | ITK multi-threading for metric evaluation |
| `--workingresolution` | `0,0,0` | Internal registration spacing in mm (`x,y,z`). Non-zero values resample fixed/moving images, and any label maps used by the metric, before optimisation while preserving original fixed-image geometry for the final written result. |

---

### 3.5 Recommended starting parameters by use case

| Use case | `--modality` | Key adjustments |
|---|---|---|
| CT–MR brain | `multimodal` | `--alpha 1.0 --lambda 0.5 --mattesnumberofbins 64` |
| MRI–MRI longitudinal | `singlemodal` | `--nu 1.0 --yota 0.5 --msepercentage 0.2` |
| CT–PET | `multimodal` | `--mattespercentage 0.2 --sigma 0.3` (add NMI for stability) |
| Anatomy-guided (with segmentations) | `custom` | `--labelkappa 0.5 --labeldistmax 20`; optionally `--labelnarrowband true --labelbandwidth 5` |
| Edge-sensitive monomodal | `custom` | `--nu 1.0 --rho 0.5 --yota 0.3` |

---

### 3.6 Monitoring per-metric contributions

During registration, `Mplus` tracks the weighted contribution of each active sub-metric (`LastValMI`, `LastValNGF`, `LastValMSE`, `LastValNC`, `LastValGD`, `LastValNMI`, `LastValLabel`, `LastValTotal`). Enable `--verbose true` to print these values at each iteration — useful for diagnosing dominance by one metric and for tuning weights.

---

## References

1. Montin E., et al. (2020). *A multi-metric registration strategy for the alignment of longitudinal brain images in pediatric oncology.* Medical & Biological Engineering & Computing, 58(4), 843–855. https://doi.org/10.1007/s11517-019-02109-4
2. Maes F., Collignon A., Vandermeulen D., Marchal G., Suetens P. (1997). *Multimodality image registration by maximization of mutual information.* IEEE Transactions on Medical Imaging, 16(2), 187–198.
3. Haber E., Modersitzki J. (2006). *Intensity gradient based registration and fusion of multi-modal images.* Methods of Information in Medicine, 45(1), 153–161.
4. Studholme C., Hill D.L.G., Hawkes D.J. (1999). *An overlap invariant entropy measure of 3D medical image alignment.* Pattern Recognition, 32(1), 71–86.
5. Insight Segmentation and Registration Toolkit (ITK). https://itk.org

---

## Part 4 — Python API, Configuration, and Customisation

### 4.1 Overview

The `mplus` Python package (`python/mplus/`) provides a high-level interface to the registration executables. The two central classes are:

| Class | File | Purpose |
|---|---|---|
| `MetricWeights` | `config.py` | Declares per-metric weights (value and derivative) |
| `RegistrationConfig` | `config.py` | Groups all registration parameters; loads/saves JSON |
| `Registration` | *(bindings)* | Wraps the C++ pipeline; accepts `RegistrationConfig` |

---

### 4.2 `MetricWeights`

```python
from mplus import MetricWeights

weights = MetricWeights(
    mi=1.0,          # Mattes MI     (--alpha / --alphaderivative)
    ngf=0.5,         # NGF           (--lambda / --lambdaderivative)
    mse=0.0,         # MSE           (--nu / --nuderivative)
    nc=0.0,          # Normalized Correlation  (--yota / --yotaderivative)
    label=0.0,       # Label-map distance      (--labelkappa / --labelkappaderiv)
    # Optional: independent derivative weights (None = same as value weight)
    mi_deriv=1.0,
    ngf_deriv=0.5,
)
```

Setting any weight to `0.0` disables that sub-metric entirely, matching the C++ behaviour.

---

### 4.3 `RegistrationConfig`

```python
from mplus import RegistrationConfig, MetricWeights

config = RegistrationConfig(
    weights=MetricWeights(mi=1.0, ngf=0.5),
    grid_spacing=30.0,          # B-spline control-point spacing (mm)
    num_levels=4,               # Pyramid levels
    iterations_per_level=200,   # Optimizer iterations per level
    derivative_mode=0,          # 0=consistent, 2=main-metric adaptive
    fixed_eta=5.0,              # NGF η for fixed image
    moving_eta=5.0,             # NGF η for moving image
    auto_estimate_eta=True,     # Auto-estimate η from image gradients
    num_samples=50000,          # Spatial samples for metric evaluation
    use_cuda=False,             # Enable CUDA if available
    label_weights=None,         # Optional Dict[int, float] per-label overrides
)
```

#### Save and load as JSON

```python
# Save
config.save("my_config.json")

# Load
config = RegistrationConfig.load("my_config.json")

# Modify and re-save
config.weights.ngf = 0.3
config.grid_spacing = 25.0
config.save("my_config_v2.json")
```

The JSON format mirrors the Python dataclass structure exactly, making configs version-controllable and shareable. The `examples/configs/` directory contains ready-to-use presets described in §4.6.

---

### 4.4 Running registration

```python
import numpy as np
from mplus import Registration, RegistrationConfig

# Load images as float32 NumPy arrays (e.g. via SimpleITK or nibabel)
fixed  = np.load("fixed.npy").astype(np.float32)
moving = np.load("moving.npy").astype(np.float32)
spacing = np.array([1.0, 1.0, 1.5])  # voxel spacing in mm (x, y, z)

config = RegistrationConfig.load("examples/configs/multimodal.json")
reg = Registration(config)
result = reg.run(fixed, moving, spacing=spacing)

print(f"Done in {result.execution_time_sec:.1f} s")
warped = result.warped_image          # float32 ndarray, same shape as fixed
params = result.transform_parameters  # 1-D array of transform parameters
```

`result.metric_values` is a dict keyed by sub-metric name (e.g., `"combined"`, `"MI"`, `"NGF"`) and is populated after the final iteration, useful for quality control.

---

### 4.5 Recommended two-stage workflow

```python
from mplus import Registration, RegistrationConfig

# Stage 1 — coarse affine alignment
affine_cfg = RegistrationConfig.load("examples/configs/affine_multilevel.json")
affine_reg = Registration(affine_cfg)
affine_result = affine_reg.run(fixed, moving, spacing=spacing)

# Stage 2 — deformable B-spline refinement on the affine-warped image
bspline_cfg = RegistrationConfig.load("examples/configs/bsplines.json")
bspline_cfg.grid_spacing = 20.0   # finer grid for this dataset
bspline_reg = Registration(bspline_cfg)
final_result = bspline_reg.run(fixed, affine_result.warped_image, spacing=spacing)
```

This mirrors the command-line workflow in §2.4 and is the pattern used throughout `examples/algorithm_examples.py`.

---

### 4.6 Bundled JSON config presets

All presets live in `examples/configs/` and can be loaded with `RegistrationConfig.load()`.

| File | Metrics active | Typical use |
|---|---|---|
| `mi_only.json` | MI | Multimodal, fast baseline |
| `multimodal.json` | MI + NGF + MSE | Challenging multimodal alignment |
| `mi_label.json` | MI + Label | Multimodal with segmentation anchors |
| `label_driven.json` | MI (0.3) + Label (1.0) | Segmentation-first alignment |
| `similarity.json` | MI | Rigid with isotropic scale (7 DOF) |
| `affine.json` | MI | Single-level affine |
| `affine_multilevel.json` | MI | Coarse-to-fine affine |
| `bsplines.json` | MI | Deformable B-spline |
| `gd_singlemodal.json` | MSE + GD | Monomodal, edge-driven |
| `nmi_multimodal.json` | MI + NMI | Multimodal, overlap-robust |

The full JSON schema is documented in `examples/configs/config_schema.json`.

---

### 4.7 Creating and extending configurations

#### Programmatic construction

```python
from mplus import MetricWeights, RegistrationConfig

# Anatomy-guided multimodal preset with per-label weights
config = RegistrationConfig(
    weights=MetricWeights(
        mi=1.0,
        ngf=0.3,
        label=0.5,
        label_deriv=0.5,
    ),
    grid_spacing=28.0,
    num_levels=4,
    iterations_per_level=200,
    label_weights={
        1: 0.8,   # ventricles — high weight
        2: 1.0,   # tumour — highest weight
        3: 0.4,   # cortex — lower weight
    },
    use_cuda=True,
)
config.save("my_research_config.json")
```

#### Derivative mode selection

`derivative_mode` controls how the per-metric gradients are merged (see §3.3):

```python
# Mode 0 — default, compatible with all optimizers
config.derivative_mode = 0

# Mode 2 — scale all sub-metric gradients to the MI gradient magnitude
config.derivative_mode = 2
# (also set weights.mi_deriv or ngf_deriv independently if needed)
```

---

### 4.8 Algorithm selection guide

| Initial misalignment | Anatomy | Recommended sequence |
|---|---|---|
| Small (< 5 mm) | Rigid (skull, bone) | `similarity` only |
| Moderate | Any | `affine_multilevel` → `bsplines` |
| Large | Any | `affine_multilevel` (4 levels) → `bsplines` |
| Any | Segmentations available | add `label` weight at B-spline stage |
| Multimodal | Any | `--modality multimodal` or `multimodal.json` |
| Monomodal | Any | `--modality singlemodal` or `gd_singlemodal.json` |

---
