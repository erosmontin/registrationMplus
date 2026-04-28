# Registration Suite - Simplified CLI v2.0

> **New in v2.0:** Simplified metric and label weight submission via arrays and presets. Full backward compatibility maintained.

---

## Quick Start

### **Simplest (Preset)**
```bash
./3DRegAffine \
  --fixed fixed_image.nii \
  --moving moving_image.nii \
  --output transform.txt \
  --preset multimodal
```

### **Explicit Weights (Array Format)**
```bash
./3DRegAffine \
  --fixed fixed_image.nii --moving moving_image.nii --output T.txt \
  --metrics "1.0,0.5,0,0,0,0,0"  # [MI, NGF, MSE, GD, NC, NMI, Label]
```

### **Traditional (Individual Parameters)**
```bash
./3DRegAffine \
  --fixed fixed_image.nii --moving moving_image.nii --output T.txt \
  --alpha 1.0 --lambda 0.5 --nu 0  # Old style still works!
```

---

## Presets

> Use `--preset` for instant metric configuration.

| Preset | Best For | Weights |
|--------|----------|---------|
| `multimodal` | Different imaging modalities (CT→MRI) | MI=1.0, NGF=0.5 |
| `singlemodal` | Same modality over time (longitudinal MRI) | MSE=1.0, NC=0.5 |
| `rigid` | Balanced structural + soft-tissue alignment | MI=1.0, NGF=0.5, MSE=0.25 |
| (empty) | Custom weights via `--metrics` | User-specified |

**Example:**
```bash
# Brain registration template-to-patient
./3DRegAffine \
  --fixed template.nii --moving patient.nii \
  --output patient_to_template.txt \
  --preset multimodal \
  --threads 8
```

---

## Metrics Format

### **Array Format: `--metrics`**

Single comma-separated array specifies all metric weights:

```
--metrics "alpha,lambda,nu,rho,yota,sigma"
```

**Layout:**
- **alpha** (0–1.5): Mattes Mutual Information (MI) weight
- **lambda** (0–1.5): Normalized Gradient Field (NGF) weight (structure)
- **nu** (0–1.5): Mean Squared Error (MSE) weight (intensity)
- **rho** (0–1.5): Gradient Difference (GD) weight
- **yota** (0–1.5): Normalized Correlation (NC) weight
- **sigma** (0–1.5): Normalized Mutual Information (NMI) weight

> **Note:** Label penalty weight (`kappa`) is separate — use `--labelkappa` or `--label-weights` instead.

**Examples:**
```bash
# MI + NGF (multimodal)
--metrics "1.0,0.5,0,0,0,0"

# MSE + NC (single-modal)
--metrics "0,0,1.0,0,0.5,0"

# All four main metrics
--metrics "1.0,0.5,0.2,0.1,0,0"

# Disable all (custom config)
--metrics "0,0,0,0,0,0"
```

### **Derivatives: `--metric-derivatives`**

Optionally specify separate derivatives (default: auto-derived from weights):

```bash
--metric-derivatives "1.0,0.5,0,0,0,0"
```

If not provided, derivatives = weights (or 0 if weight is 0).

### **Sampling Percentages: `--metric-sampling`**

Override sampling percentages per-metric (default: all 0.1 = 10%):

```bash
--metric-sampling "0.15,0.1,0.1,0.1,0.1,0.1"
```

Layout: `ma%,ngf%,mse%,gd%,nc%,nmi%`

> **Note:** Label metric sampling is separate — use `--labelsamples`.

**Example:** Use more samples for MI, less for MSE:
```bash
./3DRegAffine --fixed A.nii --moving B.nii --output T.txt \
  --preset multimodal \
  --metric-sampling "0.2,0.15,0.05,0.1,0.1,0.1"
```

---

## Label Weights

### **Scalar (Simple)**
One weight applied to ALL structures:

```bash
--labelkappa 0.5
```

**Use when:** All structures equally important, or unknown structure count

### **Vector (Advanced)**
Different weight per structure:

```bash
--label-weights "0.8,0.5,0.2"
```

For a brain with 3 structures:
- Structure 1 (hippocampus): 0.8 weight
- Structure 2 (ventricles): 0.5 weight
- Structure 3 (other): 0.2 weight

**Example:**
```bash
./3DRegAffine \
  --fixed brain.nii --fixed-labelmap brain_rois.nii \
  --moving patient.nii --moving-labelmap patient_rois.nii \
  --output transform.txt \
  --preset multimodal \
  --label-weights "0.8,0.5,0.2"
```

---

## Command Reference

### **Input/Output** (Required)
```bash
--fixed <file>              Fixed image (reference)
--moving <file>             Moving image (to be registered)
--output <file>             Output transform file
```

### **Image Processing**
```bash
--fixed-labelmap <file>     Label map for fixed image (N = none)
--moving-labelmap <file>    Label map for moving image (N = none)
```

### **Metrics (NEW - Simplified)**
```bash
--preset {multimodal|singlemodal|rigid}  Auto-configure weights
--metrics "a,b,c,d,e,f"                  Explicit weight array (alpha,lambda,nu,rho,yota,sigma)
--metric-derivatives "a,b,c,d,e,f"       Optional explicit derivatives
--metric-sampling "a,b,c,d,e,f"          Optional per-metric sampling %
```

### **Metrics (Traditional - Still Works)**
```bash
--alpha <val>                Mattes MI weight (default 1.0)
--alphaderivative <val>      MI derivative weight (default 1.0)
--lambda <val>               NGF weight (default 1.0)
--lambdaderivative <val>     NGF derivative weight (default 1.0)
--nu <val>                   MSE weight (default 1.0)
--nuderivative <val>         MSE derivative weight (default 1.0)
--normalizemse {true|false}  Normalize MSE by intensity^2 (default false)
--rho <val>                  GD weight (default 0)
--rhoderivative <val>        GD derivative weight (default 0)
--yota <val>                 NC weight (default 0)
--yotaderivative <val>       NC derivative weight (default 0)
--sigma <val>                NMI weight (default 0)
--sigmaderivative <val>      NMI derivative weight (default 0)
```

### **Label Metric (NEW)**
```bash
--label-weights "w1,w2,w3"     Per-label weights (array)
--label-derivatives "d1,d2,d3" Per-label derivatives (optional)
```

### **Label Metric (Traditional)**
```bash
--labelkappa <val>           Single weight for all labels
--labelkappaderiv <val>      Derivative (default: same as kappa)
--labelsamples <val>         Label samples per label. Fraction (0,1] OR absolute count >1 (default 0.1)
```

### **Transform & Optimizer**
```bash
--iterations <int>              Max iterations (default 1000)
--step-length <val>             Minimum step length (default 0.1)
--max-step-length <val>         Maximum step length (default 1.0)
--tolerance <val>               Gradient magnitude tolerance (default 1e-4)
--relaxation-factor <val>       Step size factor (default 0.5)
```

### **Advanced**
```bash
--threads <int>            Number of threads (default 1)
--verbose {true|false}     Verbose output (default false)
--derivativemode {0|1|2}   0=consistent, 1=normalized, 2=adaptive (default 0)
--normalizemse {true|false} Normalize MSE to [0,1] (default false)
```

---

## Examples

### **Example 1: Multi-Modal Brain Registration**
```bash
./3DRegAffine \
  --fixed CT_brain.nii \
  --moving MRI_brain.nii \
  --output CT_to_MRI.txt \
  --preset multimodal \
  --iterations 1500 \
  --threads 8
```

### **Example 2: Longitudinal MRI with ROI Constraints**
```bash
./3DRegAffine \
  --fixed baseline_mri.nii --fixed-labelmap baseline_rois.nii \
  --moving followup_mri.nii --moving-labelmap followup_rois.nii \
  --output baseline_to_followup.txt \
  --preset singlemodal \
  --label-weights "0.7,0.4,0.3" \
  --threads 4
```

### **Example 3: Custom Multi-Metric with Fine-Tuning**
```bash
./3DRegAffine \
  --fixed A.nii --moving B.nii --output T.txt \
  --metrics "1.0,0.4,0.1,0,0,0,0" \
  --metric-sampling "0.2,0.15,0.08,0.1,0.1,0.1,0.1" \
  --tolerance 1e-5 \
  --threads 12
```

### **Example 4: Old Style (Backward Compatible)**
```bash
./3DRegAffine \
  --fixed A.nii --moving B.nii --output T.txt \
  --alpha 1.0 --lambda 0.5 --nu 0 \
  --mapercentage 0.15 --ngfpercentage 0.1 --msepercentage 0.1
```

---

## Conflict Resolution

If you specify **both** array and individual parameters:

```bash
./3DRegAffine \
  --fixed A.nii --moving B.nii --output T.txt \
  --metrics "1.0,0.5,0,0,0,0,0" \
  --alpha 2.0  # ⚠️ Conflicts with array!
```

**Output:**
```
WARNING: Both --metrics array and individual parameters provided: --alpha. 
         Individual parameters take precedence.

Main Metrics:
  MI (alpha): w=2.0, d=2.0, samp=0.1         ← Overridden by --alpha
  NGF (lambda): w=0.5, d=0.5, samp=0.1      ← From array (unchanged)
  ...
```

**Resolution:** Individual parameters override array values (for selective tweaking).
To avoid warnings, use **only one format**:
- ✅ Array-only: `--metrics "1.0,0.5,..." `
- ✅ Individual-only: `--alpha 1.0 --lambda 0.5 ...`
- ⚠️ Mixed: `--metrics "..." --alpha 2.0` (works but warned)

---

## Executables

All registration tools support the new CLI:

```bash
./3DRegAffine           # Rigid affine registration
./3DRegAffineMultiLevel # Multi-resolution affine
./3DRegBsplines         # B-spline deformable registration
./3DRegSimilarity       # Metric evaluation only
```

**All accept the same CLI format.**

---

## Backward Compatibility

✅ **Fully backward compatible** — old individual parameter format still works:

```bash
# Old style (still works)
./3DRegAffine --fixed A.nii --moving B.nii --output T.txt \
  --alpha 1.0 --lambda 0.5 --nu 0 \
  --mapercentage 0.1 --ngfpercentage 0.1

# New style (cleaner)
./3DRegAffine --fixed A.nii --moving B.nii --output T.txt \
  --metrics "1.0,0.5,0,0,0,0,0"
```

Both produce identical results.

---

## FAQ

**Q: How do I know if I have more than one label in my label map?**

A: Use an image viewer (ITK-SNAP, 3D Slicer) or:
```bash
# Python
import nibabel as nib
import numpy as np
labels = np.unique(nib.load('labels.nii.gz').get_fdata())
print("Labels:", labels[labels > 0])
```

**Q: What's the difference between `--preset` and `--metrics`?**

A: 
- `--preset`: Auto-configures weights (fast, predefined)
- `--metrics`: Explicit array (full control)

Use preset for standard cases, metrics for custom tuning.

**Q: Can I mix presets and individual parameters?**

A: Yes, but later parameters override preset:
```bash
./3DRegAffine ... --preset multimodal --alpha 2.0
# Result: Uses multimodal preset, but MI weight becomes 2.0
```

**Q: Should I use new or old CLI?**

A: **New** (array format) is recommended:
- ✅ Cleaner, fewer parameters
- ✅ Less error-prone
- ✅ Auto-derives derivatives
- ✅ Same backward compatibility

Old still works perfectly if you prefer it.

**Q: Why does `--normalizemse` exist?**

A: MSE can dominate other metrics if intensity ranges are large. Normalize divides by (intensity_range)² to keep MSE in [0, 1] scale like MI/NGF/NC.

```bash
--normalizemse true   # Recommended if blending MSE with MI/NGF
```

---

## Version History

### v2.0 (Current)
- ✨ **NEW:** Array-based metric weights (`--metrics`)
- ✨ **NEW:** Metric presets (`--preset`)
- ✨ **NEW:** Unified label weight format (`--label-weights`)
- ✨ **NEW:** Conflict detection & warnings
- ✅ Full backward compatibility

### v1.0
- Individual metric parameters

---

## Support

For issues or questions:
- Check examples above
- Read `CLI_INTEGRATION_GUIDE.md` (developer reference)
- Review `MetricsConfig.h` and `LabelWeightsParser.h` (source code)
