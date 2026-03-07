# Mplus Configuration Files

This directory contains example JSON configuration files for the Mplus multi-metric image registration library.

## Quick Start

```python
from mplus import Registration, RegistrationConfig

# Load a predefined configuration
config = RegistrationConfig.load("configs/mi_only.json")

# Run registration
reg = Registration(config)
result = reg.run(fixed_image, moving_image)
```

## Available Configurations

### `mi_only.json`
**Purpose:** Traditional intensity-based registration using Mutual Information alone.

**When to use:**
- Accurate baseline intensity alignment
- When segmentation labels are not available
- Unimodal or near-unimodal image pairs

**Key parameters:**
- `mi: 1.0` – Full MI weighting
- `grid_spacing: 40.0` – Coarse control points (faster)
- `num_levels: 3` – Standard multi-resolution
- `use_cuda: false` – CPU-only for compatibility

---

### `multimodal.json`
**Purpose:** Robust multi-metric alignment combining intensity and gradient information.

**When to use:**
- Multi-modal registration (e.g., CT to MR)
- Low-contrast regions with missing intensity information
- When MI alone produces poor results

**Key parameters:**
- `mi: 1.0` – Primary intensity metric
- `ngf: 0.5` – Normalized Gradient Field (50% weight)
- `mse: 0.2` – Mean Squared Error (regularization)
- `grid_spacing: 30.0` – Finer control points
- `num_levels: 4` – Extra refinement level
- `use_cuda: true` – GPU acceleration recommended

**Expected speedup:** 20–50× on GPU vs CPU

---

### `mi_label.json`
**Purpose:** Label-guided registration with segmentation consistency.

**When to use:**
- Both images have anatomical labels/segmentations
- Want to balance intensity + label overlap
- Optimizing Dice coefficient matters

**Key parameters:**
- `mi: 1.0` – Intensity alignment
- `label: 0.5` – Label overlap (50% weight)
- `num_levels: 4` – More refinement
- `label_weights` – Per-label customization (optional)

**Example:** In the JSON, `"label_weights": {"1": 1.0, "2": 0.8}` emphasizes label 1 over label 2.

---

### `label_driven.json`
**Purpose:** Segmentation-primary registration; intensity secondary.

**When to use:**
- Label alignment critical (e.g., surgical planning)
- Intensity variations unimportant
- Low signal-to-noise intensity images

**Key parameters:**
- `mi: 0.3` – Weak intensity component (regularization)
- `label: 1.0` – Full label weighting
- `num_levels: 5` – Maximum refinement
- `grid_spacing: 25.0` – Finest control points

---

## Configuration File Format

All configs are JSON with this structure:

```json
{
  "weights": {
    "mi": 1.0,
    "ngf": 0.5,
    "mse": 0.2,
    "nc": 0.0,
    "label": 0.0,
    "mi_deriv": null,
    "ngf_deriv": null,
    "mse_deriv": null,
    "label_deriv": null
  },
  "grid_spacing": 30.0,
  "num_levels": 4,
  "iterations_per_level": 300,
  "derivative_mode": 0,
  "fixed_eta": 5.0,
  "moving_eta": 5.0,
  "auto_estimate_eta": true,
  "num_samples": 75000,
  "use_cuda": true,
  "label_weights": {
    "1": 1.0,
    "2": 0.9,
    "3": 0.95
  }
}
```

### Field Reference

| Field | Type | Range | Notes |
|-------|------|-------|-------|
| `weights.mi` | float | [0, ∞) | Mutual Information (lambda) |
| `weights.ngf` | float | [0, ∞) | Normalized Gradient Field (alpha) |
| `weights.mse` | float | [0, ∞) | Mean Squared Error (nu) |
| `weights.nc` | float | [0, ∞) | Normalized Correlation (yota) |
| `weights.label` | float | [0, ∞) | Label/Kappa metric |
| `weights.*_deriv` | float | [0, ∞) | Transform derivative weights (optional) |
| `grid_spacing` | float | (0, ∞) | Control point spacing in mm |
| `num_levels` | int | [1, 8] | Multi-resolution levels |
| `iterations_per_level` | int | [10, ∞) | Optimizer steps per level |
| `derivative_mode` | int | {0, 1, 2} | 0=weighted, 1=normalized, 2=adaptive |
| `fixed_eta` | float | [0.1, ∞) | NGF fixed-image noise param |
| `moving_eta` | float | [0.1, ∞) | NGF moving-image noise param |
| `auto_estimate_eta` | bool | true/false | Auto-compute eta from gradients |
| `num_samples` | int | [1000, ∞) | Voxels sampled per iteration |
| `use_cuda` | bool | true/false | GPU acceleration |
| `label_weights` | dict | {label: [0,1]} | Per-label Dice weights |

---

## Creating Custom Configurations

### Method 1: Copy and Edit

```bash
cp mi_only.json my_config.json
# Edit my_config.json in your favorite editor
# Then load it:
```

```python
config = RegistrationConfig.load("my_config.json")
```

### Method 2: Programmatic Creation

```python
from mplus import RegistrationConfig, MetricWeights

weights = MetricWeights(
    mi=1.0,
    ngf=0.3,
    label=0.5
)

config = RegistrationConfig(
    weights=weights,
    grid_spacing=25.0,
    num_levels=4,
    iterations_per_level=250,
    use_cuda=True
)

# Save for reuse
config.save("my_research_config.json")
```

### Method 3: Load, Modify, Save

```python
# Start from a template
config = RegistrationConfig.load("multimodal.json")

# Tweak for your data
config.grid_spacing = 20.0  # Finer control points
config.weights.ngf = 0.7    # More gradient weighting
config.use_cuda = True      # Enable GPU

# Save the modified config
config.save("my_modified_config.json")
```

---

## Parameter Tuning Guide

### Choosing Metric Weights

**Rule of thumb:** Weights don't need to sum to 1; optimizer scales them internally.

- **MI-only:** Set `mi=1.0`, others to `0.0`
- **Balanced combination:** Use `1.0` for primary, `0.3–0.7` for secondary
- **Label-critical:** Set `label > 0.5 × mi`
- **Gradient-important:** `ngf = 0.3–0.7 × mi`

### Grid Spacing

| Spacing | Speed | Accuracy | Use Case |
|---------|-------|----------|----------|
| 60+ mm | Very fast | Coarse (for quick tests) | Preprocessing tests |
| 40 mm | Fast | Good | Default, most cases |
| 25–30 mm | Moderate | Better | Multi-metric or labels |
| 15–20 mm | Slow | Excellent | Final accurate registration |
| <15 mm | Very slow | Maximum | Not recommended; use parametric post-registration |

### Multi-Resolution Levels

- **1 level:** Fast but prone to local minima (not recommended)
- **3 levels:** Standard (coarse → medium → fine)
- **4–5 levels:** More refinement, helps with difficult alignments
- **>5 levels:** Diminishing returns; may not converge faster overall

### Iterations Per Level

- **100:** Fast (but may not converge)
- **200:** Standard default
- **300–500:** Better convergence for difficult cases
- **>500:** Diminishing returns; check optimizer convergence

### Number of Samples

- **~5% of volume voxels:** `num_samples = 0.05 × sizeX × sizeY × sizeZ`
- **Small images (128³):** 100K samples OK
- **Typical images (256³):** 50K–100K samples
- **Large images (512³+):** 200K+ samples; GPU recommended

### NGF Parameters (fixed_eta, moving_eta)

- **Default (5.0):** Works for most medical images
- **Noisy images:** Increase to 10–20 (more blur tolerance)
- **Very clean images:** Decrease to 1–3 (stricter gradients)
- **auto_estimate_eta: true:** Usually best; overrides if set manually

---

## Validation

Check JSON syntax:
```bash
# Using Python
python -m json.tool mi_only.json
# Or using jq (if installed)
jq . mi_only.json
```

Validate against schema:
```python
import json
import jsonschema

with open("config_schema.json") as schema_file:
    schema = json.load(schema_file)

with open("my_config.json") as config_file:
    config = json.load(config_file)

jsonschema.validate(config, schema)
print("✓ Config is valid!")
```

---

## Common Issues

### "Metric value = NaN"
- **Cause:** Extreme weight imbalances (e.g., `ngf=100, mi=0.001`)
- **Fix:** Keep weights in comparable ranges: `1.0, 0.1–1.0, 0.1–0.5`

### "Registration diverges"
- **Cause:** Too aggressive (fine grid, many samples, high iterations)
- **Fix:** Start with coarser grid (50+ mm), fewer samples (10K), default iterations

### "Segmentation fault when using -W"
- **Cause:** (Fixed) Previously, passing a transform of a different type (e.g. `Similarity3DTransform` from `3DRegSimilarity` into `3DRegAffine`) would crash
- **Fix:** Now auto-converted. Update your executables. Any linear ITK transform (Similarity, Rigid, Euler, etc.) can be passed to `-W` on affine executables.

### "Slow convergence"
- **Cause:** Grid spacing too fine for initial alignment
- **Fix:** Use multi-resolution (num_levels ≥ 3); start coarse, refine

### "GPU out of memory"
- **Cause:** Too many samples on large image
- **Fix:** Reduce `num_samples` or use smaller voxel region of interest

---

## Example Workflows

### Workflow 1: Quick Baseline (5 min)
```json
{
  "weights": {"mi": 1.0},
  "grid_spacing": 50.0,
  "num_levels": 2,
  "iterations_per_level": 100,
  "use_cuda": true
}
```

### Workflow 2: Production (30 min)
```json
{
  "weights": {"mi": 1.0, "ngf": 0.3},
  "grid_spacing": 30.0,
  "num_levels": 4,
  "iterations_per_level": 250,
  "num_samples": 100000,
  "use_cuda": true
}
```

### Workflow 3: Research Accuracy (60+ min)
```json
{
  "weights": {"mi": 1.0, "ngf": 0.5, "label": 0.4},
  "grid_spacing": 20.0,
  "num_levels": 5,
  "iterations_per_level": 400,
  "num_samples": 200000,
  "use_cuda": true
}
```

---

## See Also
- [JSON Config Examples](../json_config_examples.py) – Code examples
- [Config Schema](config_schema.json) – Formal specification
- Full API docs: `help(RegistrationConfig)`
