# CLI Default Values Audit

## 🚨 Issues Found (Especially B-splines)

### **Critical Issues**

#### 1. **3DRegBsplines: Metric Weights Are Zero**

| Parameter | 3DRegAffine | 3DRegBsplines | Issue |
|-----------|-------------|---------------|-------|
| `lambda` (NGF) | **1.0** | **0** ⚠️ | NGF disabled by default! |
| `nu` (MSE) | **1.0** | **0** ⚠️ | MSE disabled by default! |
| `alpha` (MI) | 1.0 | 1.0 | ✅ OK |

**Impact:** B-splines only uses MI by default, no structural (NGF) or intensity (MSE) metrics. This is likely unintended.

**Lines in 3DRegBsplines.cxx:**
- Line 94: `("lambda,l", po::value<double>()->default_value(0), ...)`
- Line 100: `("nu,n", po::value<double>()->default_value(0), ...)`

---

#### 2. **3DRegBsplines: Overlap Padding Too Small**

| Parameter | 3DRegAffine | 3DRegBsplines | Issue |
|-----------|-------------|---------------|-------|
| `overlappadding` | 20 | **1** ⚠️ | B-splines needs more padding! |

**Impact:** Control points at image borders may not have enough support. B-spline knot mesh needs wider domain.

**Line 153 in 3DRegBsplines.cxx:**
```cpp
("overlappadding", po::value<unsigned int>()->default_value(1), ...)
```

The comment says it needs `>= 3` (spline order), so 1 is dangerously low.

---

#### 3. **3DRegBsplines: Convergence Factor Way Too Loose**

| Parameter | Value | Issue |
|-----------|-------|-------|
| `costfunctionconvergencefactor` | **1.e12** ⚠️ | Very loose criterion! |

**Impact:** Optimizer may stop prematurely with poor alignment. L-BFGS-B won't refine enough.

**Line 104 in 3DRegBsplines.cxx:**
```cpp
("costfunctionconvergencefactor,F", po::value<double>()->default_value(1.e12), ...)
```

Suggested: `1.e7` (moderate) or `1.e1` (high precision)

---

### **Minor Issues**

#### 4. **3DRegBsplines: Missing NGF Config**

B-splines doesn't set `NGFPrecomputeGradient` or `NGFSpacing` defaults as clearly as others. When `lambda=0`, these don't matter, but if user enables NGF they should expect sensible defaults.

#### 5. **3DRegAffineMultiLevel: Inconsistent Lambda**

Actually, 3DRegAffineMultiLevel looks OK (lambda=1.0), but worth confirming it's NOT also set to 0somewhere.

---

## ✅ What's Correct

| Parameter | Default | Status |
|-----------|---------|--------|
| `alpha` (MI) | 1.0 | ✅ Good across all |
| `nu` deriv | 1.0 | ✅ Good (Affine/MultiLevel) |
| `mattespercentage` | 0.1 (10%) | ✅ Reasonable |
| `mattesnumberofbins` | 64 | ✅ Standard |
| `maxnumberofiterations` | 1000 | ✅ OK |
| `minimumsteplength` | 0.1 | ✅ Conservative |
| `maximumsteplength` | 1.0 | ✅ Reasonable |
| `relaxationfactor` | 0.5 | ✅ Standard RSGD default |
| `gradientmagnitudetolerance` | 1e-4 | ✅ Good |
| `ngfspacing` | "4,4,4" | ✅ Standard |
| `overlappadding` (Affine) | 20 | ✅ Good |
| `metricpadding` (Bsplines) | 0 | ✅ OK (separate) |
| `labelkappa` | 0.0 | ✅ Off by default (good) |
| Sampling %: MI,NGF,MSE,etc | 0.1 (all) | ✅ Consistent |

---

## Recommended Fixes

### **Fix 1: 3DRegBsplines Lambda**

**Current (Wrong):**
```cpp
("lambda,l", po::value<double>()->default_value(0), "lambda value NGF 1.0")
```

**Should be:**
```cpp
("lambda,l", po::value<double>()->default_value(0.5), "lambda value NGF 0.5 (was 0, re-enabled for deformable)")
```

---

### **Fix 2: 3DRegBsplines Nu**

**Current (Wrong):**
```cpp
("nu,n", po::value<double>()->default_value(0), "nu value MSE 1.0")
```

**Should be:**
```cpp
("nu,n", po::value<double>()->default_value(0), "nu value MSE 0 (deformable typically uses MI+NGF, not MSE)")
// OR
("nu,n", po::value<double>()->default_value(0.25), "nu value MSE 0.25 (light regularization)")
```

**Decision:** Keep at 0, but clarify comment. B-splines are typically rigid-affine + deformable, so MI+NGF makes sense. If user wants MSE they can set `--nu 0.5`.

---

### **Fix 3: 3DRegBsplines Overlap Padding**

**Current (Wrong):**
```cpp
("overlappadding", po::value<unsigned int>()->default_value(1),
    "Number of B-spline control points outside the image domain per side "
    "(min = spline order = 3). Higher values give more deformation support at image borders.")
```

**Should be:**
```cpp
("overlappadding", po::value<unsigned int>()->default_value(5),
    "Number of B-spline control points outside the image domain per side "
    "(min = spline order = 3, recommended >= 5 for B-splines). Higher values give more deformation support at image borders.")
```

---

### **Fix 4: 3DRegBsplines Cost Function Convergence Factor**

**Current (Loose):**
```cpp
("costfunctionconvergencefactor,F", po::value<double>()->default_value(1.e12), 
   "CostFunctionConvergenceFactor 1e+12 for low accuracy; 1e+7 for moderate accuracy and 1e+1 for extremely high accuracy.")
```

**Should be:**
```cpp
("costfunctionconvergencefactor,F", po::value<double>()->default_value(1.e7), 
   "CostFunctionConvergenceFactor 1e+12 for low accuracy; 1e+7 for moderate accuracy and 1e+1 for high accuracy.")
```

---

## Summary Table

| Executable | Issue | Severity | Fix |
|------------|-------|----------|-----|
| 3DRegBsplines | lambda=0 | 🔴 High | Change to 0.5 or document decision |
| 3DRegBsplines | nu=0 | 🟡 Medium | Keep as 0, but clarify comment |
| 3DRegBsplines | overlappadding=1 | 🔴 High | Change to 5 |
| 3DRegBsplines | costfunctionconvergencefactor=1e12 | 🔴 High | Change to 1e7 |
| Others | All OK | ✅ Green | No changes |

---

## Files to Modify

1. **[3DRegBsplines.cxx](src/3DRegistration/3DRegBsplines/src/3DRegBsplines.cxx)** — Lines: 94, 100, 153, 104
