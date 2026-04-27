# Metric Parameter Resolution Strategy

## **Three Ways to Handle Conflicts**

When user provides both array (`--metrics`) and individual (`--alpha`) parameters:

### **Strategy 1: Silent Override (Default, Current)**
**Precedence:** Individual > Array
```bash
./3DRegAffine --metrics "1.0,0.5,0,0,0,0" --alpha 2.0
# Result: alpha=2.0 (overrides array value of 1.0)
# User sees: array applied, then individual overrides silently
```
**Pros:** Flexible, allows selective tweaking
**Cons:** User might not realize they're overriding

---

### **Strategy 2: Warn on Conflict (Recommended)**
Uses `DetectConflicts()` to warn user but still apply override:

```cpp
// In your CLI parsing code:
bool hasMetricsArray = (!vm["metrics"].as<std::string>().empty() || 
                        !vm["metric-derivatives"].as<std::string>().empty());

// Detect conflicts
std::string warning = MetricsConfig::DetectConflicts(
    hasMetricsArray,
    vm["alpha"].as<double>(),      // -1 if not provided
    vm["lambda"].as<double>(),     // -1 if not provided
    // ... etc
    true  // verbose=true: print warning
);

// Then merge (individual still takes precedence)
MainMetricsConfig cfg = MetricsConfig::MergeIndividual(baseConfig, ...);
```

**Output:**
```
WARNING: Both --metrics array and individual parameters provided: --alpha, --lambda. 
         Individual parameters take precedence.

Main Metrics:
  MI (alpha): w=2.0, d=2.0, samp=0.1
  NGF (lambda): w=0.3, d=0.3, samp=0.1
  MSE (nu): w=0.0, d=0.0, samp=0.1
  ...
```

**Pros:** Clear user feedback, still flexible
**Cons:** None really

---

### **Strategy 3: Strict Mode (Forbid Conflicts)**
Raise error if both array and individual are provided:

```cpp
if (hasMetricsArray) {
    if (alpha >= 0 || lambda >= 0 || nu >= 0 || ...) {
        throw std::runtime_error("Cannot mix --metrics array with individual --alpha, --lambda, etc. "
                               "Use either array OR individual, not both.");
    }
}
```

**Pros:** Clear intent, prevents mistakes
**Cons:** Less flexible, forces choose one submission method

---

## **Recommended Approach: Strategy 2**

- User gets **warned** when they mix formats
- User can **still** override selectively if intentional
- **Clear output** shows final config
- **Backward compatible** with individual-only submissions

---

## **Usage Examples**

### **Example 1: Array-only (no conflicts)**
```bash
./3DRegAffine --fixed A.nii --moving B.nii --output T.txt \
  --metrics "1.0,0.5,0,0,0,0" \
  --metric-sampling "0.1,0.1,0.1,0.1,0.1,0.1"

# Output:
# Main Metrics:
#   MI (alpha): w=1.0, d=1.0, samp=0.1
#   NGF (lambda): w=0.5, d=0.5, samp=0.1
#   ... (rest zeros)
```

---

### **Example 2: Individual-only (no conflicts)**
```bash
./3DRegAffine --fixed A.nii --moving B.nii --output T.txt \
  --alpha 1.0 --alphaderivative 1.0 --mapercentage 0.15 \
  --lambda 0.5 --lambdaderivative 0.5 --ngfpercentage 0.1

# Output:
# Main Metrics:
#   MI (alpha): w=1.0, d=1.0, samp=0.15
#   NGF (lambda): w=0.5, d=0.5, samp=0.1
#   ... (rest default to 0)
```

---

### **Example 3: Conflict with warning (Strategy 2)**
```bash
./3DRegAffine --fixed A.nii --moving B.nii --output T.txt \
  --metrics "1.0,0.5,0,0,0,0" \
  --alpha 2.0 \
  --lambda 0.3

# Output:
# 
# WARNING: Both --metrics array and individual parameters provided: --alpha, --lambda. 
#          Individual parameters take precedence.
# 
# Main Metrics:
#   MI (alpha): w=2.0, d=2.0, samp=0.1        ← alpha overridden to 2.0
#   NGF (lambda): w=0.3, d=0.3, samp=0.1     ← lambda overridden to 0.3
#   MSE (nu): w=0.0, d=0.0, samp=0.1          ← MSE from array (not overridden)
#   ... (rest default)
```

---

## **Implementation in 3DRegAffine.cxx**

```cpp
#include "../MetricsConfig.h"

// ... in main() after parsing command-line ...

bool hasMetricsArray = (!vm["metrics"].as<std::string>().empty() ||
                        !vm["metric-derivatives"].as<std::string>().empty());

// Get base config (from array or defaults)
MetricsConfig::MainMetricsConfig metricsConfig;
if (!vm["metrics"].as<std::string>().empty()) {
    metricsConfig = MetricsConfig::ParseMainWeights(vm["metrics"].as<std::string>());
} else if (!vm["preset"].as<std::string>().empty()) {
    metricsConfig = MetricsConfig::GetMainPreset(vm["preset"].as<std::string>());
} else {
    metricsConfig = MetricsConfig::MainMetricsConfig();  // All defaults (zeros)
}

// **Detect and warn on conflicts**
MetricsConfig::DetectConflicts(
    hasMetricsArray,
    vm.count("alpha") ? vm["alpha"].as<double>() : -1,
    vm.count("lambda") ? vm["lambda"].as<double>() : -1,
    vm.count("nu") ? vm["nu"].as<double>() : -1,
    vm.count("rho") ? vm["rho"].as<double>() : -1,
    vm.count("yota") ? vm["yota"].as<double>() : -1,
    vm.count("sigma") ? vm["sigma"].as<double>() : -1,
    true  // verbose: print warnings
);

// **Apply individual overrides (if provided)**
metricsConfig = MetricsConfig::MergeIndividual(
    metricsConfig,
    vm.count("alpha") ? vm["alpha"].as<double>() : -1,
    vm.count("alphaderivative") ? vm["alphaderivative"].as<double>() : -1,
    vm.count("lambda") ? vm["lambda"].as<double>() : -1,
    vm.count("lambdaderivative") ? vm["lambdaderivative"].as<double>() : -1,
    vm.count("nu") ? vm["nu"].as<double>() : -1,
    vm.count("nuderivative") ? vm["nuderivative"].as<double>() : -1,
    vm.count("rho") ? vm["rho"].as<double>() : -1,
    vm.count("rhoderivative") ? vm["rhoderivative"].as<double>() : -1,
    vm.count("yota") ? vm["yota"].as<double>() : -1,
    vm.count("yotaderivative") ? vm["yotaderivative"].as<double>() : -1,
    vm.count("sigma") ? vm["sigma"].as<double>() : -1,
    vm.count("sigmaderivative") ? vm["sigmaderivative"].as<double>() : -1
);

// **Apply any metric-specific sampling overrides**
if (!vm["metric-sampling"].as<std::string>().empty()) {
    metricsConfig = MetricsConfig::ParseMainSampling(
        metricsConfig,
        vm["metric-sampling"].as<std::string>()
    );
}

// **Print final config to user**
metricsConfig.Print("  ");

// **Apply to metric**
metric->SetAlpha(metricsConfig.mi.weight);
metric->SetAlphaDerivative(metricsConfig.mi.derivative);
metric->SetLambda(metricsConfig.ngf.weight);
metric->SetLambdaDerivative(metricsConfig.ngf.derivative);
// ... etc for all metrics ...

const unsigned int numberOfPixels = fixedImage->GetLargestPossibleRegion().GetNumberOfPixels();
metric->SetMANumberOfSamples(static_cast<unsigned int>(numberOfPixels * metricsConfig.mi.samplingPercent));
metric->SetNGFNumberOfSamples(static_cast<unsigned int>(numberOfPixels * metricsConfig.ngf.samplingPercent));
// ... etc ...
```

---

## **Decision Tree**

```
User input: --metrics X && --alpha Y?
  ├─ NO  → Use whichever is provided (individual or array)
  │        No warning needed
  │
  └─ YES → DetectConflicts() warns user
            MergeIndividual() applies override
            Print final config
            Proceed (Strategy 2)
```
