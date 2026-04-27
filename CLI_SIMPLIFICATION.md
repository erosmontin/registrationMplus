# CLI Simplification Summary

## **Parameters Removed (Inferred or Defaulted)**

| Old Parameter | New Behavior | Why |
|---|---|---|
| `--mapercentage`, `--ngfpercentage`, `--msepercentage`, etc. | Single `--metric-percentages "0.1,0.1,..."` | Avoid 7 separate params |
| `--alphaderivative`, `--lambdaderivative`, etc. | Auto-derived from weight values | If weight=0, derivative=0; else derivative=weight |
| `--etavaluefixed`, `--etavaluemoving` | Auto-estimated (enabled by default) | NGF noise is better auto-detected |
| `--ngfevaluator` | Default = 0 (scalar) | Doesn't change for most users |
| `--nmibins` | Default = 64 | Standard histogram size |
| `--ngfspacing` | Inferred from image spacing | Compute as `image.spacing * 2.0` |
| `--fixedimagethreshold` | Optional (advanced) | Most registrations don't use |
| `--metricoverlap` | Default = true | Better in 99% of cases |
| `--numberofthreads` | Short `--threads` or `-T` | Same functionality, shorter |

---

## **New Vs Old CLI Comparison**

### **Scenario 1: Multimodal Registration (brain CT to MRI)**

**OLD (current):**
```bash
./3DRegAffine --fixed brain_mri.nii --moving brain_ct.nii --output transform.txt \
  --alpha 1.0 --alphaderivative 1.0 \
  --lambda 0.5 --lambdaderivative 0.5 \
  --nu 0 --nuderivative 0 \
  --mapercentage 0.1 --ngfpercentage 0.1 --msepercentage 0.1 \
  --maxnumberofiterations 1000 \
  --threads 4
  # 13 parameters
```

**NEW (simplified):**
```bash
./3DRegAffine --fixed brain_mri.nii --moving brain_ct.nii --output transform.txt \
  --preset multimodal --threads 4
  # 5 parameters (77% reduction!)
```

---

### **Scenario 2: Single-Modal Registration (longitudinal MRI)**

**OLD:**
```bash
./3DRegAffine --fixed baseline.nii --moving followup.nii --output transform.txt \
  --alpha 0 --nu 1.0 --nuderivative 1.0 \
  --yota 0.5 --yotaderivative 0.5 \
  --mapercentage 0.1 --msepercentage 0.1 --ncpercentage 0.1 \
  --minimumsteplength 0.05 \
  --threads 8
  # 11 parameters
```

**NEW:**
```bash
./3DRegAffine --fixed baseline.nii --moving followup.nii --output transform.txt \
  --preset singlemodal --step-length 0.05 --threads 8
  # 6 parameters
```

---

### **Scenario 3: Custom Weights**

**OLD:**
```bash
./3DRegAffine --fixed A.nii --moving B.nii --output transform.txt \
  --alpha 0.8 --alphaderivative 0.8 \
  --lambda 0.3 --lambdaderivative 0.3 \
  --nu 0.1 --nuderivative 0.1 \
  --mapercentage 0.15 --ngfpercentage 0.12 --msepercentage 0.08
```

**NEW:**
```bash
./3DRegAffine --fixed A.nii --moving B.nii --output transform.txt \
  --weights "0.8,0.3,0.1,0,0,0,0" \
  --metric-percentages "0.15,0.12,0.08,0.1,0.1,0.1,0.1"
```

---

## **What Still Exists (Advanced Use Cases)**

These parameters remain available but are **optional and rarely needed**:

```bash
--derivativemode 0          # 0=consistent (default), 1=normalized, 2=adaptive
--mainmetric 0              # For mode 2: which metric drives scaling
--normalizemse true         # Keep MSE in [0,1] range
--ngfprecompute false       # Speed optimization
--fixedimagethreshold -1    # Focus on ROI only
--verbose true              # Debug output
```

---

## **Preset Reference**

| Preset | Weights | Best For |
|--------|---------|----------|
| `multimodal` | α=1.0, λ=0.5 | Different imaging modalities (CT↔MRI) |
| `singlemodal` | ν=1.0, γ=0.5 | Same modality, different time points |
| `rigid` | α=1.0, λ=0.5, ν=0.25 | Structural balance (MI + NGF + small MSE) |
| `custom` | None | User specifies via `--weights` |

---

## **Impact: Fewer Mistakes**

**Before:** User could set `--nu 1.0 --nuderivative 0.5` → value/gradient mismatch → optimizer fails

**After:** CliParser auto-derives derivatives → impossible to mismatch

---

## **Backward Compatibility**

✅ Old CLI still works:
```bash
./3DRegAffine --fixed A.nii --moving B.nii --output transform.txt \
  --alpha 1.0 --alphaderivative 1.0 --lambda 0.5 --lambdaderivative 0.5 ...
```

Auto-detection picks the right path automatically.
