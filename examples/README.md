# Registration Examples & Documentation

Complete guide to all registration algorithms, configurations, and usage patterns for the Mplus multi-metric image registration library.

---

## 📚 Quick Navigation

### 🎯 If you need...

| Need | File | Best For |
|------|------|----------|
| **Algorithm comparison** | [QUICK_REFERENCE.md](QUICK_REFERENCE.md) | 1-minute overview |
| **How to choose algorithm** | [ALGORITHM_GUIDE.md](#algorithm-selection-decision-tree) (Decision Tree section) | Deciding which algorithm to use |
| **Detailed algorithm info** | [ALGORITHM_GUIDE.md](ALGORITHM_GUIDE.md) | In-depth reading |
| **Code examples** | [algorithm_examples.py](algorithm_examples.py) | Running examples in Python |
| **JSON configurations** | [configs/](configs/) | Ready-to-use config files |
| **Config parameters explained** | [configs/README.md](configs/README.md) | Understanding each parameter |
| **Parameter tuning guide** | [configs/README.md](#parameter-tuning-guide) (Parameter Tuning section) | Fine-tuning for your data |

---

## 🚀 Getting Started (30 seconds)

### Option A: Use a Preset Config
```python
from mplus import Registration, RegistrationConfig
import numpy as np

# Load a predefined configuration
config = RegistrationConfig.load("configs/affine_multilevel.json")

# Run registration
reg = Registration(config)
result = reg.run(fixed_image, moving_image, spacing=np.array([1, 1, 1]))

print(f"Done in {result.execution_time_sec:.1f} sec")
```

### Option B: Customize a Config
```python
# Start from template
config = RegistrationConfig.load("configs/bsplines.json")

# Customize
config.grid_spacing = 30.0      # Finer control points
config.num_levels = 4           # Better convergence
config.use_cuda = True          # GPU acceleration

# Run
reg = Registration(config)
result = reg.run(fixed_image, moving_image)
```

### Option C: Command Line (C++)
```bash
cd ../bld2/bin
./3DRegAffineMultiLevel fixed.nii moving.nii warped.nii
```

---

## 📋 Available Registration Algorithms

### 1. **Affine** (Single-Level)
- **File:** `configs/affine.json`
- **Executable:** `bld2/bin/3DRegAffine`
- **Transform:** Rotation + Scaling + Shearing (12 DOF)
- **Speed:** ⚡⚡⚡ Very Fast (~30 sec)
- **Accuracy:** ⭐⭐ Moderate (~10 mm TRE)
- **When to use:** Quick initialization, global alignment
- **Python:** See [algorithm_examples.py](algorithm_examples.py) → `example_affine()`

### 2. **Affine Multi-Level** ⭐ Recommended Global
- **File:** `configs/affine_multilevel.json`
- **Executable:** `bld2/bin/3DRegAffineMultiLevel`
- **Transform:** Affine (12 DOF) + Coarse-to-fine pyramid
- **Speed:** ⚡⚡ Fast (~1-5 min)
- **Accuracy:** ⭐⭐⭐ Good (~5-6 mm TRE)
- **When to use:** Robust global alignment, pre-registration for B-splines
- **Python:** See [algorithm_examples.py](algorithm_examples.py) → `example_affine_multilevel()`

### 3. **B-Spline** ⭐ Recommended Non-Rigid
- **File:** `configs/bsplines.json`
- **Executable:** `bld2/bin/3DRegBsplines`
- **Transform:** Free-Form Deformation (1000s DOF)
- **Speed:** 🐢 Slow CPU / ⚡ Fast GPU (~5-30 min with GPU)
- **Accuracy:** ⭐⭐⭐⭐⭐ Excellent (~1-2 mm TRE)
- **When to use:** Non-rigid deformations, anatomy with local shape variation
- **GPU Note:** **9× faster on GPU!** Highly recommended.
- **Workflow:** Use Affine MultiLevel as pre-registration, then B-Spline
- **Python:** See [algorithm_examples.py](algorithm_examples.py) → `example_bsplines()`

### 4. **Similarity**
- **File:** `configs/similarity.json`
- **Executable:** `bld2/bin/3DRegSimilarity`
- **Transform:** Rotation + Isotropic Scaling (7 DOF)
- **Speed:** ⚡⚡⚡ Very Fast (~20 sec)
- **Accuracy:** ⭐⭐ Moderate (~5-8 mm TRE)
- **When to use:** Rigid bodies (bones, skull), when fewer DOF is better
- **Python:** See [algorithm_examples.py](algorithm_examples.py) → `example_similarity()`

---

## 📊 Algorithm Comparison Table

```
┌────────────────┬──────────┬──────────┬──────────┬──────────┬─────────┐
│ Algorithm      │ DOF      │ Speed    │ GPU Help │ Accuracy │ Robust  │
├────────────────┼──────────┼──────────┼──────────┼──────────┼─────────┤
│ Affine         │ 12       │ ⚡⚡⚡    │ 1.5×     │ ⭐⭐    │ ✓       │
│ Affine ML      │ 12       │ ⚡⚡     │ 2×       │ ⭐⭐⭐  │ ✓✓      │
│ B-Spline       │ 1000s    │ 🐢⚡     │ 9×       │ ⭐⭐⭐⭐⭐│ ✓       │
│ Similarity     │ 7        │ ⚡⚡⚡    │ 1.3×     │ ⭐⭐    │ ✓✓✓     │
└────────────────┴──────────┴──────────┴──────────┴──────────┴─────────┘

⚡ = CPU speed level
🐢⚡ = Slow on CPU, fast on GPU
✓✓✓ = Very robust
```

---

## 🎯 Decision Tree: Choosing an Algorithm

```
START
  │
  ├─→ Is it a rigid body (bone, skull)?
  │   YES → Use SIMILARITY (7 DOF, most robust)
  │   NO  ↓
  │
  ├─→ Are images roughly aligned already?
  │   YES → Use AFFINE (single-level, fast)
  │   NO  ↓
  │
  ├─→ Large initial misalignment?
  │   YES → Use AFFINE_MULTILEVEL (coarse-to-fine)
  │   NO  ↓
  │
  └─→ Need sub-mm accuracy or non-rigid deformations?
      YES → Use BSPLINES (enable GPU!)
      NO  → Use AFFINE_MULTILEVEL (good balance)
```

---

## 💡 Recommended Workflows

### Workflow 1: Speed (1-2 min)
```python
config = RegistrationConfig.load("configs/affine.json")
result = Registration(config).run(fixed, moving)
# Time: ~30 sec | TRE: ~10 mm
```

### Workflow 2: Balance (5-10 min)
```python
config = RegistrationConfig.load("configs/affine_multilevel.json")
result = Registration(config).run(fixed, moving)
# Time: ~3 min | TRE: ~5-6 mm
```

### Workflow 3: High Precision (15-30 min)
```python
# Stage 1: Global alignment (similarity)
config_sim = RegistrationConfig.load("configs/similarity.json")
result_sim = Registration(config_sim).run(fixed, moving)

# Stage 2: Refine with affine, warm-started from similarity transform
config_aff = RegistrationConfig.load("configs/affine_multilevel.json")
config_aff.transform_in = result_sim.transform_path  # auto-converts Similarity3D → Affine
result_aff = Registration(config_aff).run(fixed, moving)

# Stage 3: Local deformation (GPU)
config_bs = RegistrationConfig.load("configs/bsplines.json")
config_bs.use_cuda = True
result_bs = Registration(config_bs).run(fixed, result_aff.warped_image)
# Time: ~15 min | TRE: <2 mm
```

> **Note:** The `-W` / `transform_in` flag on affine executables accepts any linear ITK transform file
> (e.g. `Similarity3DTransform` from `3DRegSimilarity`) and auto-converts it to `AffineTransform`.

### Workflow 4: Multi-Metric (15-20 min)
```python
# Combine multiple information sources
config = RegistrationConfig.load("configs/bsplines.json")
config.weights.mi = 1.0         # Intensity (Mattes MI)
config.weights.ngf = 0.5        # Gradient (NGF)
config.weights.label = 0.3      # Segmentation (label distance)
config.weights.gd = 0.0         # Gradient Difference: set > 0 for singlemodal edge-based alignment
config.weights.nmi = 0.0        # Normalized MI: set > 0 for multimodal histogram-based alignment
config.use_cuda = True

result = Registration(config).run(
    fixed, moving,
    fixed_labels=fixed_seg,
    moving_labels=moving_seg
)
```

---

## 📂 File Structure

```
examples/
│
├── QUICK_REFERENCE.md              ← Start here! One-page cheat sheet
├── ALGORITHM_GUIDE.md              ← Comprehensive algorithm details
├── algorithm_examples.py           ← Python code examples + tutorial
├── json_config_examples.py         ← How to load/save/customize configs
│
└── configs/
    ├── README.md                   ← Config file reference
    ├── config_schema.json          ← JSON schema (for validation)
    │
    ├── affine.json                 ← Single-level affine
    ├── affine_multilevel.json      ← Multi-level affine ⭐
    ├── bsplines.json               ← B-spline deformable ⭐
    ├── similarity.json             ← Similarity transform
    │
    ├── mi_only.json                ← MI metric alone
    ├── multimodal.json             ← MI + NGF + MSE (advanced)
    ├── mi_label.json               ← MI + label overlap
    ├── label_driven.json           ← Label-primary (advanced)
    ├── gd_singlemodal.json         ← GD only (Gradient Difference, singlemodal structural)
    └── nmi_multimodal.json         ← NMI + MI (Normalized MI, multimodal)
```

---

## 🔧 How to Use Configuration Files

### Load a Config
```python
from mplus import RegistrationConfig

config = RegistrationConfig.load("configs/affine_multilevel.json")
```

### Modify and Save
```python
config.num_levels = 4
config.grid_spacing = 30.0
config.use_cuda = True
config.save("my_custom_config.json")
```

### Validate JSON
```bash
python -m json.tool configs/affine.json
# or
jq . configs/affine.json
```

---

## 📖 Important Sections

### In ALGORITHM_GUIDE.md:

1. **Quick Reference** – Algorithm comparison table
2. **Algorithm Details** – Deep dive into each algorithm
3. **Decision Tree** – How to choose
4. **Complete Workflows** – Real-world examples
5. **Performance Benchmarks** – Speed/accuracy data

### In configs/README.md:

1. **Configuration Overview** – Each config explained
2. **Field Reference** – Parameter meanings & ranges
3. **Parameter Tuning** – How to adjust for your data
4. **Common Issues** – Troubleshooting
5. **Example Workflows** – Quick, production, research setups

### In QUICK_REFERENCE.md:

- One-page summary
- Performance table
- Quick start code
- Decision tree
- Common issues

---

## 🎓 Learning Path

**5 minutes:** Read [QUICK_REFERENCE.md](QUICK_REFERENCE.md)

**30 minutes:** Read [ALGORITHM_GUIDE.md](ALGORITHM_GUIDE.md) → Algorithm Details section

**1 hour:** Run examples from [algorithm_examples.py](algorithm_examples.py)

**2 hours:** Read config parameters in [configs/README.md](configs/README.md) and customize

**Full mastery:** Experiment with different configs on your own data

---

## ⚡ Performance Tips

1. **Use GPU for B-Splines:** 5–20× faster than CPU
2. **Use Affine MLi for initialization:** Quick + robust
3. **Reduce `num_samples` for speed:** Trade-off metric smoothness
4. **Increase `num_levels` for robustness:** Coarse-to-fine convergence
5. **Two-stage workflow:** Affine → B-Spline = best results

---

## 🆘 Troubleshooting

| Issue | Likely Cause | Solution |
|-------|--------------|----------|
| "Metric = NaN" | Weights too extreme | Use balanced weights (1.0, 0.3, 0.1) |
| "Registration diverges" | Grid too fine | Increase `grid_spacing` or add levels |
| "Slow convergence" | Too few iterations | Increase `iterations_per_level` |
| "GPU out of memory" | Too many samples | Reduce `num_samples` or use ROI |
| "Folding/checkerboard" | Over-deformation | Increase `grid_spacing` |

See [configs/README.md](configs/README.md) → **Common Issues** for more.

---

## 📚 See Also

- Python API: `help(Registration)`, `help(RegistrationConfig)`
- C++ Docs: `src/Metrics/Mplus/itkMplus.h`
- ITK Documentation: https://itk.org/
- Mplus v2 Roadmap: `Mplusv2roadmap.md`

---

## 🎉 Quick Commands

```bash
# View a config
cat configs/affine_multilevel.json

# Validate a config
python -m json.tool configs/bsplines.json out.json

# Run Python examples
python algorithm_examples.py

# Run C++ examples
cd ../bld2/bin
./3DRegBsplines fixed.nii moving.nii warped.nii
```

---

**Last Updated:** March 6, 2026  
**Version:** Mplus v2 Documentation  
**Status:** Complete & Ready to Use ✅
