# 3DRegBsplines - QUICK START CARD

## 🎯 The Simplest Command

```bash
./bin/3DRegBsplines -f fixed.nii.gz -m moving.nii.gz -o result.nii.gz
```

**That's it.** Everything else is optional.

---

## 📊 What Happens Inside

```
Your Command ──────────────────────────────────────────────────
    ↓
Load Images
    ↓
Initialize B-spline Transform (50mm mesh, padding=5)
    ↓
Initialize Optimizer (L-BFGS-B)
    ↓
Compute Initial Metrics
    ├─ MI (1.0)    ← Intensity matching
    ├─ NGF (0.5)   ← Structure preservation (NOW ENABLED!)
    └─ MSE (0.0)   ← Off
    ↓
Iterate (up to 1000 times):
    ├─ Compare fixed vs moving
    ├─ Update deformation
    ├─ Check convergence (1.e7 tolerance)
    ├─ if converged → STOP
    └─ else → continue
    ↓
Output:
    └─ result.nii.gz (Registered moving image)
```

---

## 📈 Quality vs Speed Trade-off

```
Speed/Quality Continuum:

⚡ FAST                                               🎯 PRECISE
├─────────────────────────────────────────────────┤
1-2min    5-15min          20-60min
Rough     Standard ← DEFAULT → High Precision
          (50mm)    (Fixed!)   (25mm)
```

---

## 🔧 Three Knobs to Adjust

### 1. Mesh Resolution (`--gridresolution`)
```
Coarse (100mm)  ←────────→  Fine (25mm)
Speed: ⚡⚡⚡              Speed: 🐢
Accuracy: ⭐              Accuracy: ⭐⭐⭐⭐⭐
DEFAULT: 50mm (balanced)
```

### 2. Iterations (`--maxnumberofiterations`)
```
Few (200)     ←────────→  Many (2000)
Speed: ⚡⚡⚡              Speed: 🐢
Accuracy: ⭐              Accuracy: ⭐⭐⭐⭐⭐
DEFAULT: 1000 (balanced)
```

### 3. Structure Weight (`--lambda`)
```
Low (0.3)     ←────────→  High (0.8)
Structure:    Weak         Strong
DEFAULT: 0.5 (balanced)
```

---

## 🚀 Copy-Paste Commands

### FASTEST (1-2 min)
```bash
./bin/3DRegBsplines -f fixed.nii.gz -m moving.nii.gz -o result.nii.gz \
  --gridresolution 100 --maxnumberofiterations 200
```

### BALANCED (5-15 min) ← DEFAULT
```bash
./bin/3DRegBsplines -f fixed.nii.gz -m moving.nii.gz -o result.nii.gz
```

### HIGHEST QUALITY (20-60 min)
```bash
./bin/3DRegBsplines -f fixed.nii.gz -m moving.nii.gz -o result.nii.gz \
  --gridresolution 25 --maxnumberofiterations 2000
```

---

## 🎯 Add Outputs

### Output Deformation Field (vector image)
```bash
./bin/3DRegBsplines -f fixed.nii.gz -m moving.nii.gz -o result.nii.gz \
  --vfout deformation.nii.gz
```

### Output Transform (reusable)
```bash
./bin/3DRegBsplines -f fixed.nii.gz -m moving.nii.gz -o result.nii.gz \
  --transformout transform.tfm
```

### Save Progress (PNG snapshots every 10 iterations)
```bash
./bin/3DRegBsplines -f fixed.nii.gz -m moving.nii.gz -o result.nii.gz \
  --snapshotdir ./snapshots --snapshotevery 10 --verbose 1
```

---

## ✅ Verify Installation

```bash
./bin/3DRegBsplines --help | head -20
```

Should show:
```
--lambda (=0.5)              ← ⭐ NEW DEFAULT
--costfunctionconvergencefactor (=1e+07)  ← ⭐ NEW DEFAULT
--overlappadding (=5)        ← ⭐ NEW DEFAULT
```

---

## 🆘 If Something Goes Wrong

| Problem | Solution |
|---------|----------|
| Poor alignment | `--lambda 0.7 --gridresolution 40` |
| Too much deformation | `--lambda 0.3 --gridresolution 70` |
| Too slow | `--gridresolution 100 --maxnumberofiterations 200` |
| Bad boundaries | `--overlappadding 7` |
| Not converging | `--maxnumberofiterations 2000` |

---

## 📚 For More Info

- **All examples** → `3DRegBsplines_CLI_GUIDE.md` (13 detailed examples)
- **Quick reference** → `3DRegBsplines_QUICK_REFERENCE.md` (copy-paste)
- **What changed** → `3DRegBsplines_BEFORE_AND_AFTER.md` (fixes explained)
- **Full audit** → `CLI_DEFAULTS_AUDIT.md` (complete analysis)

---

## 🎬 Typical Workflow

```
Step 1: Rough alignment (Affine)
$ ./bin/3DRegAffine -f fixed.nii.gz -m moving.nii.gz -o affine.nii.gz
(2-3 minutes)

Step 2: Fine refinement (B-spline)
$ ./bin/3DRegBsplines -f fixed.nii.gz -m affine.nii.gz -o final.nii.gz
(10 minutes)

Step 3: Check result
$ # Open final.nii.gz in viewer
$ # Compare with fixed.nii.gz
```

---

## 💡 Pro Tips

✅ **DO:**
- Use two-stage registration (Affine → B-spline)
- Start with defaults, then adjust
- Save deformation field for reuse
- Monitor convergence with verbose mode

❌ **DON'T:**
- Change all parameters at once
- Use very large lambda (>1.0)
- Forget to rebuild after code changes
- Mix different image resolutions

---

## 🎯 Default Parameters (FIXED)

| What | Value | Why |
|------|-------|-----|
| NGF weight | 0.5 | ⭐ Structure preservation |
| Convergence | 1.e7 | ⭐ Proper optimization |
| Padding | 5 | ⭐ Boundary support |
| Mesh | 50mm | Balanced speed/accuracy |
| Max iterations | 1000 | Sufficient for convergence |
| Sampling | 10% | Speed vs accuracy |

---

## ⏱️ Expected Runtimes (for 256×256×256 images)

| Config | Time | Quality |
|--------|------|---------|
| Fast (100mm, 200 iter) | 1-2 min | ⭐ |
| Default (50mm, 1000 iter) | 5-15 min | ⭐⭐⭐⭐ |
| Precise (25mm, 2000 iter) | 30-60 min | ⭐⭐⭐⭐⭐ |

*Depends on: image size, GPU/CPU, convergence speed*

---

## 🔗 Full Command Template

```bash
./bin/3DRegBsplines \
  --fixedimage FIXED.nii.gz \
  --movingimage MOVING.nii.gz \
  --outputimage RESULT.nii.gz \
  \
  # Optional: outputs
  --vfout DEFORMATION.nii.gz \
  --transformout TRANSFORM.tfm \
  --snapshotdir ./snapshots \
  \
  # Optional: tuning
  --lambda 0.5 \
  --gridresolution 50 \
  --maxnumberofiterations 1000 \
  --verbose 1
```

---

**Ready? Run this:**
```bash
./bin/3DRegBsplines -f fixed.nii.gz -m moving.nii.gz -o result.nii.gz
```

🎉 That's it!
