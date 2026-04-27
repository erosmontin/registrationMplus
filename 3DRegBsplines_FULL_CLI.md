# 3DRegBsplines - COMPLETE CLI WITH ALL PARAMETERS

## 📋 ABSOLUTE MINIMAL

```bash
./bin/3DRegBsplines \
  --fixedimage fixed.nii.gz \
  --movingimage moving.nii.gz \
  --outputimage result.nii.gz
```

✅ Uses all NEW defaults automatically

---

## 📊 COMPLETE COMMAND (Every Parameter)

```bash
./bin/3DRegBsplines \
  \
  # ═══════════════════════════════════════════════════════════════
  # REQUIRED INPUTS (Must specify)
  # ═══════════════════════════════════════════════════════════════
  --fixedimage fixed.nii.gz \
  --movingimage moving.nii.gz \
  --outputimage result.nii.gz \
  \
  # ═══════════════════════════════════════════════════════════════
  # OPTIONAL OUTPUTS (Default: N = off)
  # ═══════════════════════════════════════════════════════════════
  --vfout deformation_field.nii.gz \
  --transformout transform.tfm \
  --transformin N \
  --gridposition N \
  --snapshotdir ./snapshots \
  \
  # ═══════════════════════════════════════════════════════════════
  # METRIC WEIGHTS (PRIMARY - NEW DEFAULTS!)
  # ═══════════════════════════════════════════════════════════════
  --alpha 1.0 \                    # Mutual Information (MI) weight
  --lambda 0.5 \                   # ⭐ NGF weight (NEW: was 0) - STRUCTURE
  --nu 0.0 \                       # Sum Squared Differences (MSE) - off for deformable
  --rho 0.0 \                      # Gradient Difference (GD) - off
  --yota 0.0 \                     # Normalized Correlation (NC) - off
  --sigma 0.0 \                    # Normalized Mutual Information (NMI) - off
  \
  # ═══════════════════════════════════════════════════════════════
  # METRIC DERIVATIVES (Analytical gradients of metrics)
  # ═══════════════════════════════════════════════════════════════
  --alphaderivative 1.0 \          # MI derivative weight
  --lambdaderivative 0 \           # NGF derivative - off (numerical ok)
  --nuderivative 0 \               # MSE derivative - off
  --rhoderivative 0.0 \            # GD derivative - off
  --yotaderivative 0 \             # NC derivative - off
  --sigmaderivative 0.0 \          # NMI derivative - off
  \
  # ═══════════════════════════════════════════════════════════════
  # MI/MATTES CONFIGURATION
  # ═══════════════════════════════════════════════════════════════
  --mattespercentage 0.1 \         # Use 10% of image pixels for MI
  --mattesnumberofbins 64 \        # MI histogram bins (64-256 typical)
  --bsplinecaching true \          # Cache MI B-spline interpolation
  --explicitPDFderivatives false \ # Use implicit PDF derivatives
  \
  # ═══════════════════════════════════════════════════════════════
  # NGF CONFIGURATION (Normalized Gradient Field)
  # ═══════════════════════════════════════════════════════════════
  --ngfspacing 4,4,4 \             # NGF gradient spacing (x,y,z in mm)
  --ngfpercentage 0.1 \            # Use 10% of pixels for NGF
  --NGFevaluator 0 \               # NGF type: 0=scalar, 1=cross, 2=scdelta, 3=Delta, 4=Delta2
  --etavaluefixed -1 \             # Fixed image noise estimate (-1 = auto)
  --etavaluemoving -1 \            # Moving image noise estimate (-1 = auto)
  --ngfprecompute false \          # Recompute NGF every iteration (vs precompute once)
  \
  # ═══════════════════════════════════════════════════════════════
  # MSE CONFIGURATION
  # ═══════════════════════════════════════════════════════════════
  --msepercentage 0.1 \            # Use 10% of pixels for MSE
  --normalizemse false \           # Normalize by intensity range (usually off)
  \
  # ═══════════════════════════════════════════════════════════════
  # GRADIENT DIFFERENCE (GD)
  # ═══════════════════════════════════════════════════════════════
  --gdpercentage 0.1 \             # Use 10% of pixels for GD
  \
  # ═══════════════════════════════════════════════════════════════
  # NORMALIZED CORRELATION (NC)
  # ═══════════════════════════════════════════════════════════════
  --ncpercentage 0.1 \             # Use 10% of pixels for NC
  \
  # ═══════════════════════════════════════════════════════════════
  # NORMALIZED MUTUAL INFORMATION (NMI)
  # ═══════════════════════════════════════════════════════════════
  --nmipercentage 0.1 \            # Use 10% of pixels for NMI
  --nmibins 64 \                   # NMI histogram bins
  \
  # ═══════════════════════════════════════════════════════════════
  # MULTI-METRIC CONFIGURATION
  # ═══════════════════════════════════════════════════════════════
  --derivativemode 0 \             # Derivative merge: 0=consistent, 2=adaptive
  --mainmetric 0 \                 # Main metric for mode 2: 0=MI, 1=NGF, 2=MSE, 3=NC
  --metricoverlap true \           # Compute overlap between fixed/moving
  \
  # ═══════════════════════════════════════════════════════════════
  # B-SPLINE MESH & OPTIMIZATION (NEW DEFAULTS HIGHLIGHTED!)
  # ═══════════════════════════════════════════════════════════════
  --gridresolution 50 \            # B-spline mesh spacing (mm) [coarse:100, fine:25]
  --overlappadding 5 \             # ⭐ Control points outside domain (NEW: was 1)
  --meshmarginsize 0.0 \           # Extra margin (mm) to extend mesh
  --maxnumberofiterations 1000 \   # Max optimization iterations
  --costfunctionconvergencefactor 1.e7 \  # ⭐ Convergence criterion (NEW: was 1.e12)
  --projectedgradienttolerance 1.e-5 \    # Gradient magnitude tolerance
  --numberofevaluations 500 \      # L-BFGS-B: max function evaluations
  --numberofcorrections 5 \        # L-BFGS-B: Hessian correction pairs
  \
  # ═══════════════════════════════════════════════════════════════
  # BOUNDARY CONDITIONS & CONSTRAINTS
  # ═══════════════════════════════════════════════════════════════
  --fixedimagethreshold -99999999 \# Ignore fixed pixels below this value
  --dfltpixelvalue 0 \             # Default for resampled pixels
  --bound 0 \                      # Boundary: 0=unbounded, 1=lower, 2=both, 3=upper
  --lbound 0 \                     # Lower bound for parameters
  --ubound 0 \                     # Upper bound for parameters
  \
  # ═══════════════════════════════════════════════════════════════
  # LABEL/SEGMENTATION SUPPORT
  # ═══════════════════════════════════════════════════════════════
  --fixedlabelmap N \              # Fixed segmentation (N=none)
  --movinglabelmap N \             # Moving segmentation (N=none)
  --labelkappa 0.0 \               # Label metric weight (0=off)
  --labelkappaderiv 0.0 \          # Label metric derivative weight
  --labelsamples 0.1 \             # Use 10% of labeled pixels
  --labelreport 1 \                # Report Dice every N iterations (0=off)
  \
  # ═══════════════════════════════════════════════════════════════
  # VISUALIZATION & DEBUGGING
  # ═══════════════════════════════════════════════════════════════
  --snapshotevery 1 \              # Save snapshot every N iterations
  --snapshotstack false \          # PNG vs full 3D .nii.gz
  --snapshotgrid false \           # Show deformation grid (vs knot mesh)
  --snapshotgridspacing 20 \       # Grid line spacing (voxels)
  --verbose false \                # Print verbose output
  \
  # ═══════════════════════════════════════════════════════════════
  # SYSTEM & PERFORMANCE
  # ═══════════════════════════════════════════════════════════════
  --numberofthreads 2              # Number of CPU threads
```

---

## 🎯 ORGANIZED BY USE CASE

### 1️⃣ MULTIMODAL (CT + MR, PET, etc.)

```bash
./bin/3DRegBsplines \
  --fixedimage ct.nii.gz \
  --movingimage mr.nii.gz \
  --outputimage ct_mr_aligned.nii.gz \
  \
  # Metrics for different modalities
  --alpha 1.0 \                    # MI handles modality differences well
  --lambda 0.5 \                   # NGF preserves edges
  --nu 0.0 \                       # MSE not useful across modalities
  \
  # Configuration
  --mattespercentage 0.1 \         # Standard MI sampling
  --gridresolution 50 \            # Standard mesh
  --maxnumberofiterations 1000
```

---

### 2️⃣ SAME MODALITY (MR to MR, CT to CT)

```bash
./bin/3DRegBsplines \
  --fixedimage baseline_mr.nii.gz \
  --movingimage followup_mr.nii.gz \
  --outputimage followup_aligned.nii.gz \
  \
  # Can use intensity difference metric
  --alpha 1.0 \                    # MI
  --lambda 0.8 \                   # Stronger structure guidance
  --nu 0.25 \                      # Enable MSE for intensity similarity
  \
  # Configuration
  --mattespercentage 0.15 \        # Slightly higher sampling
  --msepercentage 0.1 \
  --gridresolution 40 \            # Finer mesh for better detail
  --maxnumberofiterations 1500
```

---

### 3️⃣ WITH SEGMENTATION LABELS

```bash
./bin/3DRegBsplines \
  --fixedimage fixed.nii.gz \
  --movingimage moving.nii.gz \
  --outputimage result.nii.gz \
  --fixedlabelmap fixed_organs.nii.gz \
  --movinglabelmap moving_organs.nii.gz \
  \
  # Image metrics
  --alpha 1.0 \
  --lambda 0.5 \
  \
  # Label correspondence
  --labelkappa 1.0 \               # Weight for label alignment
  --labelsamples 0.5 \             # Use 50% of labeled voxels
  --labelreport 50 \               # Report Dice every 50 iterations
  \
  # Configuration
  --gridresolution 45 \
  --maxnumberofiterations 1200
```

---

### 4️⃣ WITH TRANSFORMATION OUTPUT

```bash
./bin/3DRegBsplines \
  --fixedimage fixed.nii.gz \
  --movingimage moving.nii.gz \
  --outputimage result.nii.gz \
  \
  # Save outputs for reuse
  --transformout bspline_transform.tfm \    # Reusable transform
  --vfout deformation_field.nii.gz \        # Displacement vectors
  \
  --alpha 1.0 \
  --lambda 0.5 \
  --gridresolution 50
```

---

### 5️⃣ WITH PROGRESS VISUALIZATION

```bash
./bin/3DRegBsplines \
  --fixedimage fixed.nii.gz \
  --movingimage moving.nii.gz \
  --outputimage result.nii.gz \
  \
  # Snapshots for monitoring
  --snapshotdir ./registration_progress \
  --snapshotevery 10 \             # Save PNG every 10 iterations
  --snapshotgrid true \            # Show deformation grid
  --snapshotgridspacing 15 \
  --verbose true \                 # Print iteration info
  \
  --alpha 1.0 \
  --lambda 0.5 \
  --gridresolution 50
```

---

### 6️⃣ HIGH PRECISION

```bash
./bin/3DRegBsplines \
  --fixedimage fixed.nii.gz \
  --movingimage moving.nii.gz \
  --outputimage result.nii.gz \
  \
  # Metrics
  --alpha 1.0 \
  --lambda 0.6 \                   # Slightly stronger structure
  --mattespercentage 0.2 \         # More pixels (20%)
  --ngfpercentage 0.2 \
  \
  # More refinement
  --gridresolution 25 \            # Fine mesh
  --overlappadding 7 \             # Extra boundary support
  --maxnumberofiterations 2000 \
  --costfunctionconvergencefactor 1.e1 \   # Tight convergence
  --projectedgradienttolerance 1.e-6 \     # Strict gradient tolerance
  --numberofevaluations 1000 \
  --numberofcorrections 10
```

---

### 7️⃣ FAST (ROUGH ALIGNMENT)

```bash
./bin/3DRegBsplines \
  --fixedimage fixed.nii.gz \
  --movingimage moving.nii.gz \
  --outputimage result.nii.gz \
  \
  # Metrics
  --alpha 1.0 \
  --lambda 0.3 \                   # Less structure (faster)
  --mattespercentage 0.05 \        # Fewer pixels (5%)
  \
  # Quick convergence
  --gridresolution 100 \           # Coarse mesh
  --maxnumberofiterations 200 \
  --costfunctionconvergencefactor 1.e12 \  # Loose convergence
  --numberofevaluations 100 \
  --numberofcorrections 2
```

---

### 8️⃣ CUSTOM MULTI-METRIC

```bash
./bin/3DRegBsplines \
  --fixedimage fixed.nii.gz \
  --movingimage moving.nii.gz \
  --outputimage result.nii.gz \
  \
  # Combine multiple metrics
  --alpha 1.0 \                    # MI (intensity)
  --alphaderivative 1.0 \
  --lambda 0.5 \                   # NGF (structure)
  --lambdaderivative 0 \
  --nu 0.2 \                       # MSE (intensity difference)
  --nuderivative 0 \
  --rho 0.1 \                      # GD (gradient difference)
  --rhoderivative 0 \
  \
  # Sampling percentages per metric
  --mattespercentage 0.1 \         # MI sampling
  --ngfpercentage 0.1 \            # NGF sampling
  --msepercentage 0.1 \            # MSE sampling
  --gdpercentage 0.1 \             # GD sampling
  \
  # Multi-metric configuration
  --derivativemode 0 \             # Consistent derivatives
  --mainmetric 0 \                 # Primary: MI
  \
  # Optimization
  --gridresolution 40 \
  --maxnumberofiterations 1500 \
  --costfunctionconvergencefactor 1.e7
```

---

## 🔥 PARAMETER CATEGORIES REFERENCE

### MUST SPECIFY (No defaults)
```bash
--fixedimage FILE        # Required
--movingimage FILE       # Required
--outputimage FILE       # Required
```

### ALMOST ALWAYS USEFUL
```bash
--alpha 1.0              # MI weight (use in almost all cases)
--lambda 0.5             # ⭐ NGF (NEW default - use in almost all cases)
--gridresolution 50      # Mesh fineness (main tuning knob)
--maxnumberofiterations 1000  # Convergence control
```

### OPTIONAL - USUALLY DON'T CHANGE
```bash
--nu 0.0                 # MSE (off for deformable)
--mattespercentage 0.1   # 10% pixels is good default
--mattesnumberofbins 64  # Standard histogram bins
--overlappadding 5       # ⭐ NEW safe default
--costfunctionconvergencefactor 1.e7  # ⭐ NEW good default
--projectgradienttolerance 1.e-5      # Almost never needs change
```

### ADVANCED - RARELY USED
```bash
--rho 0.0                # GD metric (off)
--yota 0.0               # NC metric (off)
--sigma 0.0              # NMI metric (off)
--etavaluefixed -1       # Auto-detect noise
--etavaluemoving -1      # Auto-detect noise
--bound 0                # Boundary conditions
--meshmarginsize 0.0     # Mesh extension
```

### OPTIONAL OUTPUTS
```bash
--transformout FILE      # Save B-spline transform for reuse
--vfout FILE             # Save deformation field (vector image)
--snapshotdir DIR        # Save progress PNGs
```

---

## 📊 DEFAULT VALUES TABLE (NEW DEFAULTS!)

| Parameter | Default | Type | Category |
|-----------|---------|------|----------|
| `--alpha` | 1.0 | double | MI weight |
| `--lambda` | **0.5** ⭐ | double | NGF weight (NEW) |
| `--nu` | 0.0 | double | MSE weight |
| `--gridresolution` | 50 | double | Mesh (mm) |
| `--overlappadding` | **5** ⭐ | int | Boundary (NEW) |
| `--costfunctionconvergencefactor` | **1.e7** ⭐ | double | Convergence (NEW) |
| `--maxnumberofiterations` | 1000 | int | Max iterations |
| `--mattespercentage` | 0.1 | double | MI sampling (10%) |
| `--mattesnumberofbins` | 64 | int | MI histogram |
| `--ngfpercentage` | 0.1 | double | NGF sampling (10%) |
| `--msepercentage` | 0.1 | double | MSE sampling (10%) |
| `--projectedgradienttolerance` | 1.e-5 | double | Gradient tol |
| `--numberofevaluations` | 500 | int | L-BFGS-B evals |
| `--numberofcorrections` | 5 | int | L-BFGS-B corr |
| `--numberofthreads` | 2 | int | CPU threads |

---

## 🎬 TYPICAL USAGE SEQUENCES

### Sequence 1: ONE-LINER
```bash
./bin/3DRegBsplines -f fixed.nii.gz -m moving.nii.gz -o result.nii.gz
```
Uses all fixed NEW defaults ✅

### Sequence 2: STANDARD
```bash
./bin/3DRegBsplines \
  --fixedimage fixed.nii.gz \
  --movingimage moving.nii.gz \
  --outputimage result.nii.gz \
  --alpha 1.0 \
  --lambda 0.5 \
  --gridresolution 50 \
  --maxnumberofiterations 1000
```
Explicit standard configuration

### Sequence 3: TWO-STAGE (Recommended)
```bash
# Stage 1: Affine alignment
./bin/3DRegAffine \
  -f fixed.nii.gz -m moving.nii.gz \
  -o affine_result.nii.gz

# Stage 2: B-spline refinement
./bin/3DRegBsplines \
  -f fixed.nii.gz -m affine_result.nii.gz \
  -o final_result.nii.gz \
  --lambda 0.5 \
  --gridresolution 40 \
  --maxnumberofiterations 1500
```

### Sequence 4: BATCH PROCESSING
```bash
for img in *.nii.gz; do
  echo "Processing $img..."
  ./bin/3DRegBsplines \
    --fixedimage template.nii.gz \
    --movingimage "$img" \
    --outputimage "aligned_$img" \
    --alpha 1.0 \
    --lambda 0.5 \
    --gridresolution 50
done
```

---

## 🔑 KEY PARAMETERS TO ADJUST

| Problem | Adjustment |
|---------|-------------|
| Poor alignment | `--lambda 0.7` (more structure), `--gridresolution 40` (finer) |
| Too much deformation | `--lambda 0.3` (less structure), `--gridresolution 70` (coarser) |
| Slow/timeout | `--gridresolution 100` (coarser), `--maxnumberofiterations 200` |
| Boundary artifacts | `--overlappadding 8` (more support) |
| Not converging | `--maxnumberofiterations 2000`, `--costfunctionconvergencefactor 1.e1` |

---

## ✅ VERIFICATION

```bash
# Check if fixes are in place
./bin/3DRegBsplines --help | grep -A1 "lambda\|convergence\|padding"
```

Should show:
```
--lambda (=0.5)               ✅
--costfunctionconvergencefactor (=1e+07)  ✅
--overlappadding (=5)         ✅
```
