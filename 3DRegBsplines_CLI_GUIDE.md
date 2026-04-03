# 3DRegBsplines CLI Usage Guide

## NEW DEFAULTS (After Fixes)

| Parameter | Default | What It Does |
|-----------|---------|--------------|
| `--lambda` | **0.5** | NGF weight (structure preservation) |
| `--alpha` | **1.0** | MI weight (intensity matching) |
| `--nu` | **0** | MSE weight (off - use NGF+MI instead) |
| `--costfunctionconvergencefactor` | **1.e7** | Convergence precision (moderate accuracy) |
| `--overlappadding` | **5** | Control point boundary support |
| `--gridresolution` | **50 mm** | B-spline mesh spacing |
| `--maxnumberofiterations` | **1000** | Max optimization steps |

---

## 1. BASIC USAGE (Uses ALL Fixed Defaults)

```bash
./bin/3DRegBsplines \
  --fixedimage    fixed.nii.gz \
  --movingimage   moving.nii.gz \
  --outputimage   result.nii.gz
```

**What happens:**
- ✅ MI (1.0) + NGF (0.5) metric enabled
- ✅ Moderate convergence criterion (1e7)
- ✅ Safe boundary control points (padding=5)
- ✅ 50mm mesh resolution
- ✅ Up to 1000 iterations

---

## 2. QUICK DEFORMATION (Fast, Lower Quality)

```bash
./bin/3DRegBsplines \
  --fixedimage    fixed.nii.gz \
  --movingimage   moving.nii.gz \
  --outputimage   result.nii.gz \
  --gridresolution 100 \
  --maxnumberofiterations 200 \
  --costfunctionconvergencefactor 1.e12
```

**Use when:** Speed is priority, rough alignment OK
- Coarser mesh (100mm vs 50mm)
- Fewer iterations (200 vs 1000)
- Looser convergence (1e12 vs 1e7)
- **Time:** ~1-2 minutes

---

## 3. HIGH PRECISION REGISTRATION (Slow, Best Quality)

```bash
./bin/3DRegBsplines \
  --fixedimage    fixed.nii.gz \
  --movingimage   moving.nii.gz \
  --outputimage   result.nii.gz \
  --gridresolution 25 \
  --maxnumberofiterations 2000 \
  --costfunctionconvergencefactor 1.e1 \
  --overlappadding 7
```

**Use when:** Maximum alignment accuracy needed
- Fine mesh (25mm vs 50mm)
- Many iterations (2000 vs 1000)
- Tight convergence (1e1 vs 1e7)
- Extra boundary support (7 vs 5)
- **Time:** ~10-30 minutes

---

## 4. MULTIMODAL REGISTRATION (MR + CT, PET, etc.)

```bash
./bin/3DRegBsplines \
  --fixedimage    ct.nii.gz \
  --movingimage   mr.nii.gz \
  --outputimage   ct_to_mr.nii.gz \
  --alpha 1.0 \
  --lambda 0.5 \
  --mattespercentage 0.1 \
  --gridresolution 50
```

**Use when:** Different imaging modalities (contrast, intensity ranges vary)
- MI (1.0) handles modality differences
- NGF (0.5) preserves edge structure
- 10% sampling for speed
- Standard 50mm mesh

---

## 5. SINGLE MODALITY REGISTRATION (MR to MR, CT to CT)

```bash
./bin/3DRegBsplines \
  --fixedimage    scan1.nii.gz \
  --movingimage   scan2.nii.gz \
  --outputimage   registered.nii.gz \
  --alpha 1.0 \
  --lambda 1.0 \
  --nu 0.25 \
  --gridresolution 40
```

**Use when:** Same modality (higher confidence in intensity)
- Same metrics (1.0, 1.0, 0.25)
- Can use MSE (0.25) for intensity differences
- Slightly finer mesh (40mm)

---

## 6. WITH INTENSITY NORMALIZATION (Brain, Skull-stripped)

```bash
./bin/3DRegBsplines \
  --fixedimage    brain_fixed.nii.gz \
  --movingimage   brain_moving.nii.gz \
  --outputimage   brain_registered.nii.gz \
  --alpha 1.0 \
  --lambda 0.8 \
  --mattespercentage 0.2 \
  --gridresolution 40 \
  --maxnumberofiterations 1500
```

**Use when:** Well-preprocessed brain images
- Higher NGF weight (0.8) for structural details
- More sampling (20% vs 10%)
- Finer mesh (40mm)
- More iterations (1500)

---

## 7. WITH LABEL SEGMENTATION (Anatomy-guided)

```bash
./bin/3DRegBsplines \
  --fixedimage    fixed.nii.gz \
  --movingimage   moving.nii.gz \
  --outputimage   result.nii.gz \
  --fixedlabelmap fixed_labels.nii.gz \
  --movinglabelmap moving_labels.nii.gz \
  --labelkappa 1.0 \
  --labelsamples 0.5 \
  --alpha 1.0 \
  --lambda 0.5 \
  --gridresolution 50
```

**Use when:** Segmentation maps available (organs, ROIs)
- Adds label correspondence metric (1.0)
- Uses 50% of labeled pixels
- Helps alignment in anatomically important regions

---

## 8. WITH DEFORMATION FIELD OUTPUT (For analysis)

```bash
./bin/3DRegBsplines \
  --fixedimage    fixed.nii.gz \
  --movingimage   moving.nii.gz \
  --outputimage   result.nii.gz \
  --vfout         deformation_field.nii.gz \
  --transformout  bspline_transform.tfm \
  --gridresolution 50
```

**Use when:** Need to analyze/visualize deformation
- VF: Deformation field (3D vector image with magnitude/direction)
- Transform: B-spline knot mesh (can reapply to new images)

---

## 9. WITH SNAPSHOTS (Progressive visualization)

```bash
./bin/3DRegBsplines \
  --fixedimage    fixed.nii.gz \
  --movingimage   moving.nii.gz \
  --outputimage   result.nii.gz \
  --snapshotdir   ./snapshots \
  --snapshotevery 10 \
  --snapshots grid true \
  --gridresolution 50
```

**Use when:** Want to see registration progress
- Saves PNG every 10 iterations
- Shows deformation grid overlay
- Helps debug convergence issues
- Directory created automatically

---

## 10. FINE-TUNING EXAMPLE (Custom metric combination)

```bash
./bin/3DRegBsplines \
  --fixedimage    fixed.nii.gz \
  --movingimage   moving.nii.gz \
  --outputimage   result.nii.gz \
  --alpha 1.0 \
  --alphaderivative 1.0 \
  --lambda 0.6 \
  --lambdaderivative 0 \
  --nu 0.2 \
  --nuderivative 0 \
  --etavaluefixed -1 \
  --etavaluemoving -1 \
  --mattespercentage 0.1 \
  --ngfpercentage 0.1 \
  --msepercentage 0.1 \
  --gridresolution 40 \
  --maxnumberofiterations 1500 \
  --costfunctionconvergencefactor 1.e7 \
  --overlappadding 5
```

**Use when:** Need precise control over all parameters
- MI (1.0) with derivatives
- NGF (0.6) for structure
- MSE (0.2) for intensity differences
- Auto-compute noise levels (-1)
- Full metric combination

---

## 11. BATCH PROCESSING EXAMPLE

```bash
#!/bin/bash
# Register multiple subjects

for subject in subj001 subj002 subj003; do
  echo "Registering $subject..."
  
  ./bin/3DRegBsplines \
    --fixedimage    template.nii.gz \
    --movingimage   data/${subject}_native.nii.gz \
    --outputimage   data/${subject}_registered.nii.gz \
    --vfout         data/${subject}_deformation.nii.gz \
    --gridresolution 50 \
    --maxnumberofiterations 1000
    
  echo "✓ $subject complete"
done
```

---

## 12. PIPELINE AFTER AFFINE (Two-stage registration)

```bash
# Stage 1: Rigid+Affine
./bin/3DRegAffine \
  --fixedimage    fixed.nii.gz \
  --movingimage   moving.nii.gz \
  --outputimage   affine_registered.nii.gz \
  --transformout  affine_transform.tfm

# Stage 2: B-spline deformable (refine affine)
./bin/3DRegBsplines \
  --fixedimage    fixed.nii.gz \
  --movingimage   affine_registered.nii.gz \
  --outputimage   final_registered.nii.gz \
  --gridresolution 40 \
  --alpha 1.0 \
  --lambda 0.5
```

**This is the recommended workflow:**
- ✅ Stage 1: Get rigid/affine alignment (fast, global)
- ✅ Stage 2: Refine with B-spline deformation (slower, local)

---

## 13. TROUBLESHOOTING COMMANDS

### Poor convergence?
```bash
./bin/3DRegBsplines \
  --fixedimage fixed.nii.gz \
  --movingimage moving.nii.gz \
  --outputimage result.nii.gz \
  --costfunctionconvergencefactor 1.e1 \
  --maxnumberofiterations 2000 \
  --verbose true
```

### Boundary artifacts?
```bash
./bin/3DRegBsplines \
  --fixedimage fixed.nii.gz \
  --movingimage moving.nii.gz \
  --outputimage result.nii.gz \
  --overlappadding 8 \
  --gridresolution 45
```

### Too much deformation (over-registration)?
```bash
./bin/3DRegBsplines \
  --fixedimage fixed.nii.gz \
  --movingimage moving.nii.gz \
  --outputimage result.nii.gz \
  --lambda 0.7 \
  --nu 0.1 \
  --gridresolution 70
```

---

## PARAMETER QUICK REFERENCE

### Metric Weights (sum ~2.0 for stability)
```
--alpha                  (MI):  1.0 (default) - intensity matching
--lambda                (NGF):  0.5 (default) - structure preservation
--nu                    (MSE):  0.0 (default) - off for deformable
--rho                    (GD):  0.0 (default) - off
--yota                   (NC):  0.0 (default) - off
--sigma                 (NMI):  0.0 (default) - off
```

### Mesh & Iterations
```
--gridresolution         50 mm (default) - coarse: 100, fine: 25
--overlappadding          5 (default) - safety: 3-7, large: 10
--maxnumberofiterations 1000 (default) - fast: 200, precise: 2000
```

### Convergence
```
--costfunctionconvergencefactor 1.e7 (default)
  Fast:     1.e12
  Moderate: 1.e7 (default)
  Precise:  1.e1
```

### Sampling
```
--mattespercentage  0.1 (10%) - MI pixel sampling
--ngfpercentage     0.1 (10%) - NGF pixel sampling
--msepercentage     0.1 (10%) - MSE pixel sampling
```

---

## PRESET SHORTCUTS (if available)

```bash
# Multimodal (MI+NGF)
./bin/3DRegBsplines \
  --fixedimage fixed.nii.gz \
  --movingimage moving.nii.gz \
  --outputimage result.nii.gz \
  --preset multimodal

# Single modality (MSE+NC)
./bin/3DRegBsplines \
  --fixedimage fixed.nii.gz \
  --movingimage moving.nii.gz \
  --outputimage result.nii.gz \
  --preset singlemodal
```

---

## TYPICAL RUNTIMES

| Scenario | Grid | Iters | Time | Quality |
|----------|------|-------|------|---------|
| Quick (Example 2) | 100mm | 200 | 1-2 min | Rough |
| Standard (Example 1) | 50mm | 1000 | 5-15 min | Good |
| High precision (Example 3) | 25mm | 2000 | 20-60 min | Excellent |

*Times depend on image size, GPU/CPU, and convergence speed*

---

## VALIDATION TIPS

After registration, check:

1. **Visual inspection:**
   ```bash
   # Look at registered image overlaid on fixed
   # Use any NIfTI viewer (ITK-SNAP, Slicer, etc.)
   ```

2. **Statistics:**
   ```bash
   # Compare intensity distributions
   # Check for artifacts at boundaries
   ```

3. **Deformation magnitude:**
   ```bash
   # Compute deformation field magnitude
   # Large values (~20+ mm) might indicate over-registration
   ```

4. **Segmentation overlap (if labels available):**
   ```bash
   # Compute Dice coefficient before/after
   # Should improve significantly
   ```
