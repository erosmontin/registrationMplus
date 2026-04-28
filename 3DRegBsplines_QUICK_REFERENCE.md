# 3DRegBsplines - Command Cheat Sheet

## 📋 Copy-Paste Ready Commands

### SIMPLEST (Just Input/Output)
```bash
./bin/3DRegBsplines -f fixed.nii.gz -m moving.nii.gz -o result.nii.gz
```

### STANDARD QUALITY
```bash
./bin/3DRegBsplines \
  --fixedimage fixed.nii.gz \
  --movingimage moving.nii.gz \
  --outputimage result.nii.gz
```

### HIGH QUALITY
```bash
./bin/3DRegBsplines \
  --fixedimage fixed.nii.gz \
  --movingimage moving.nii.gz \
  --outputimage result.nii.gz \
  --gridresolution 25 \
  --maxnumberofiterations 2000 \
  --costfunctionconvergencefactor 1.e1
```

### FAST (Less Accurate)
```bash
./bin/3DRegBsplines \
  --fixedimage fixed.nii.gz \
  --movingimage moving.nii.gz \
  --outputimage result.nii.gz \
  --gridresolution 100 \
  --maxnumberofiterations 200
```

### WITH DEFORMATION FIELD
```bash
./bin/3DRegBsplines \
  --fixedimage fixed.nii.gz \
  --movingimage moving.nii.gz \
  --outputimage result.nii.gz \
  --vfout deformation_field.nii.gz \
  --transformout transform.tfm
```

### WITH LABELS
```bash
./bin/3DRegBsplines \
  --fixedimage fixed.nii.gz \
  --movingimage moving.nii.gz \
  --outputimage result.nii.gz \
  --fixedlabelmap fixed_seg.nii.gz \
  --movinglabelmap moving_seg.nii.gz \
  --labelkappa 1.0
```

### WITH PROGRESS SNAPSHOTS
```bash
./bin/3DRegBsplines \
  --fixedimage fixed.nii.gz \
  --movingimage moving.nii.gz \
  --outputimage result.nii.gz \
  --snapshotdir ./snapshots \
  --snapshotevery 10 \
  --verbose true
```

### CUSTOM METRICS
```bash
./bin/3DRegBsplines \
  --fixedimage fixed.nii.gz \
  --movingimage moving.nii.gz \
  --outputimage result.nii.gz \
  --alpha 1.0 \
  --lambda 0.6 \
  --nu 0.2
```

### TWO-STAGE (Affine + B-spline)
```bash
# Step 1: Affine
./bin/3DRegAffine \
  --fixedimage fixed.nii.gz \
  --movingimage moving.nii.gz \
  --outputimage affine_result.nii.gz

# Step 2: B-spline (refine)
./bin/3DRegBsplines \
  --fixedimage fixed.nii.gz \
  --movingimage affine_result.nii.gz \
  --outputimage final_result.nii.gz
```

---

## 🔧 Short Options

| Short | Long | Example |
|-------|------|---------|
| `-f` | `--fixedimage` | `-f fixed.nii.gz` |
| `-m` | `--movingimage` | `-m moving.nii.gz` |
| `-o` | `--outputimage` | `-o result.nii.gz` |
| `-a` | `--alpha` | `-a 1.0` |
| `-l` | `--lambda` | `-l 0.5` |
| `-n` | `--nu` | `-n 0.0` |
| `-g` | `--gridresolution` | `-g 50` |
| `-I` | `--maxnumberofiterations` | `-I 1000` |
| `-V` | `--verbose` | `-V 1` |

**Example with short options:**
```bash
./bin/3DRegBsplines -f fixed.nii.gz -m moving.nii.gz -o result.nii.gz -l 0.6 -g 40 -I 1500
```

---

## 📊 Default Values (FIXED)

```
Metric Weights:
  --alpha 1.0          (MI - intensity)
  --lambda 0.5         (NGF - structure) ← FIXED from 0
  --nu 0.0             (MSE - off)
  
Mesh & Convergence:
  --gridresolution 50  (mm)
  --overlappadding 5   (voxels) ← FIXED from 1
  --costfunctionconvergencefactor 1.e7 ← FIXED from 1.e12
  
Optimization:
  --maxnumberofiterations 1000
  --projectedgradienttolerance 1.e-5
  --numberofevaluations 500
  --numberofcorrections 5
  
Sampling:
  --mattespercentage 0.1 (10%)
  --mattesnumberofbins 64
  --ngfpercentage 0.1
  --msepercentage 0.1
```

---

## ⚡ Common Adjustments

### If alignment is poor:
```bash
--lambda 0.7 \           # More structure guidance
--gridresolution 40 \    # Finer mesh
--maxnumberofiterations 1500
```

### If too much deformation:
```bash
--lambda 0.3 \           # Less structure guidance
--gridresolution 70 \    # Coarser mesh
--nu 0.1                 # Add MSE penalty
```

### If slow/timeout:
```bash
--gridresolution 100 \   # Coarser mesh
--maxnumberofiterations 200 \
--mattespercentage 0.05  # Fewer pixels
```

### If you have labels:
```bash
--fixedlabelmap seg_fixed.nii.gz \
--movinglabelmap seg_moving.nii.gz \
--labelkappa 1.0 \
--labelsamples 0.5            # fraction (0,1]; or pass an integer >1 for absolute count
```

---

## 🧪 Test Command

```bash
./bin/3DRegBsplines --help
```

Should show:
```
lambda value NGF 0.5 (enables structural regularization for deformable registration)
costfunctionconvergencefactor ... 1e+7 for moderate accuracy ...
overlappadding ... recommended >= 5 for B-splines ...
```

✅ If you see these → **Fixes were applied!**

---

## 📁 Output Files

**From:** `--outputimage result.nii.gz`
- ✅ result.nii.gz — Registered moving image

**From:** `--vfout deformation.nii.gz`
- ✅ deformation.nii.gz — 3D displacement field (vector image)

**From:** `--transformout transform.tfm`
- ✅ transform.tfm — B-spline transform (reusable)

**From:** `--snapshotdir ./snapshots`
- ✅ snapshots/iteration_0000.png — Progress visualizations

---

## 🔗 Full Documentation

- **CLI Guide:** `3DRegBsplines_CLI_GUIDE.md` (13 detailed examples)
- **Audit Report:** `CLI_DEFAULTS_AUDIT.md` (what was broken)
- **Fixes Applied:** `FIXES_APPLIED.md` (what was fixed)
- **Verification:** `CLI_DEFAULTS_VERIFICATION.md` (complete analysis)

---

## 💡 Tips

1. **Always** use fixed image as reference (don't swap)
2. **Test** with standard parameters first
3. **Compare** before & after with viewer
4. **Check** for boundary artifacts if increasing lambda
5. **Use** two-stage (Affine + B-spline) for best results
6. **Save** deformation field for later analysis/reuse

---

## ⚙️ Rebuild After Fixes

```bash
cd /data/PROJECTS/registrationSuite
rm -rf build && mkdir build && cd build
cmake .. && make -j$(nproc)
./bin/3DRegBsplines --help
```
