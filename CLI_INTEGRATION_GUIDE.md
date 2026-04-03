# CLI Update: Complete Integration Guide for All Executables

This guide shows how to update all 4 registration executables with the new simplified CLI format.

## Files to Update

1. `3DRegAffine/src/3DRegAffine.cxx`
2. `3DRegAffineMultiLevel/src/3DRegAffineMultiLevel.cxx`
3. `3DRegBsplines/src/3DRegBsplines.cxx`
4. `3DRegSimilarity/src/3DRegSimilarity.cxx` 

---

## Step 1: Add Includes (Same for All 4 Files)

**Location:** After `#include "../../Version.h"`

```cpp
#include "../../MetricsConfig.h"
#include "../../LabelWeightsParser.h"
```

---

## Step 2: Add New CLI Options (Same Pattern for All 4)

**Location:** In `po::options_description desc()` section, after the manual weight options (--alpha, --lambda, --nu, etc.)

### Remove duplicates:
Delete these **duplicate** lines (they appear twice):
```cpp
("msepercentage",   po::value<double>()->default_value(0.1), "MSE percentage of pixels used (0.1 = 10%)")
("ngfpercentage",   po::value<double>()->default_value(0.1), "NGF percentage of pixels used (0.1 = 10%)")
```

### Add new array-based options:
```cpp
// ───── NEW SIMPLIFIED CLI (arrays + presets) ─────
("preset", po::value<std::string>()->default_value(""), 
 "Metric preset: 'multimodal' (MI+NGF), 'singlemodal' (MSE+NC), 'rigid', or empty for custom")

("metrics", po::value<std::string>()->default_value(""), 
 "Metric weights array: alpha,lambda,nu,rho,yota,sigma (e.g., '1.0,0.5,0,0,0,0')")

("metric-derivatives", po::value<std::string>()->default_value(""), 
 "Metric derivatives array: alpha_d,lambda_d,nu_d,rho_d,yota_d,sigma_d")

("metric-sampling", po::value<std::string>()->default_value(""), 
 "Metric sampling percentages: ma%,ngf%,mse%,gd%,nc%,nmi% (label → --labelsamples)")

("label-weights", po::value<std::string>()->default_value(""), 
 "Per-label weights (alternative to --labelkappa): comma-separated list (e.g., '0.5,0.3,0.2')")

("label-derivatives", po::value<std::string>()->default_value(""), 
 "Per-label derivatives (auto-derived from weights if not provided)")
```

---

## Step 3: Add Auto-Detection Logic (Same for All 4)

**Location:** After `po::store(po::parse_command_line(argc, argv, desc), vm);` and `po::notify(vm);`

**Before any metric setup:**

```cpp
// ─────────────────────────────────────────────────────────────────────────────
//  SIMPLIFIED CLI: Auto-detect format (array vs individual) and parse metrics
// ─────────────────────────────────────────────────────────────────────────────

bool hasMetricsArray = (!vm["metrics"].as<std::string>().empty() ||
                        !vm["metric-derivatives"].as<std::string>().empty() ||
                        !vm["preset"].as<std::string>().empty());

// Get base metrics config
MetricsConfig::MainMetricsConfig metricsConfig;

if (!vm["preset"].as<std::string>().empty())
{
    metricsConfig = MetricsConfig::GetMainPreset(vm["preset"].as<std::string>());
    std::cout << "\n[CLI] Using preset: " << vm["preset"].as<std::string>() << std::endl;
}
else if (!vm["metrics"].as<std::string>().empty())
{
    metricsConfig = MetricsConfig::ParseMainWeights(vm["metrics"].as<std::string>());
    std::cout << "\n[CLI] Parsed weights array" << std::endl;
}

// Apply explicit metric derivatives if provided
if (!vm["metric-derivatives"].as<std::string>().empty())
{
    metricsConfig = MetricsConfig::ParseMainDerivatives(
        metricsConfig,
        vm["metric-derivatives"].as<std::string>()
    );
}

// Warn on conflicts and merge individual parameter overrides
MetricsConfig::DetectConflicts(
    hasMetricsArray,
    (vm.count("alpha") && vm["alpha"].as<double>() != 1.0) ? vm["alpha"].as<double>() : -1,
    (vm.count("lambda") && vm["lambda"].as<double>() != 1.0) ? vm["lambda"].as<double>() : -1,
    (vm.count("nu") && vm["nu"].as<double>() != 1.0) ? vm["nu"].as<double>() : -1,
    vm.count("rho") ? vm["rho"].as<double>() : -1,
    vm.count("yota") ? vm["yota"].as<double>() : -1,
    vm.count("sigma") ? vm["sigma"].as<double>() : -1,
    true  // verbose
);

// Apply individual overrides
metricsConfig = MetricsConfig::MergeIndividual(
    metricsConfig,
    (vm.count("alpha") && vm["alpha"].as<double>() >= 0) ? vm["alpha"].as<double>() : -1,
    (vm.count("alphaderivative") && vm["alphaderivative"].as<double>() >= 0) ? vm["alphaderivative"].as<double>() : -1,
    (vm.count("lambda") && vm["lambda"].as<double>() >= 0) ? vm["lambda"].as<double>() : -1,
    (vm.count("lambdaderivative") && vm["lambdaderivative"].as<double>() >= 0) ? vm["lambdaderivative"].as<double>() : -1,
    (vm.count("nu") && vm["nu"].as<double>() >= 0) ? vm["nu"].as<double>() : -1,
    (vm.count("nuderivative") && vm["nuderivative"].as<double>() >= 0) ? vm["nuderivative"].as<double>() : -1,
    (vm.count("rho") && vm["rho"].as<double>() >= 0) ? vm["rho"].as<double>() : -1,
    (vm.count("rhoderivative") && vm["rhoderivative"].as<double>() >= 0) ? vm["rhoderivative"].as<double>() : -1,
    (vm.count("yota") && vm["yota"].as<double>() >= 0) ? vm["yota"].as<double>() : -1,
    (vm.count("yotaderivative") && vm["yotaderivative"].as<double>() >= 0) ? vm["yotaderivative"].as<double>() : -1,
    (vm.count("sigma") && vm["sigma"].as<double>() >= 0) ? vm["sigma"].as<double>() : -1,
    (vm.count("sigmaderivative") && vm["sigmaderivative"].as<double>() >= 0) ? vm["sigmaderivative"].as<double>() : -1
);

// Apply metric-specific sampling overrides
if (!vm["metric-sampling"].as<std::string>().empty())
{
    metricsConfig = MetricsConfig::ParseMainSampling(
        metricsConfig,
        vm["metric-sampling"].as<std::string>()
    );
}

// Print final configuration to user
std::cout << "\n[Metrics Configuration]" << std::endl;
metricsConfig.Print("  ");

// Handle label weights (separate from main metrics)
LabelWeightsParser::LabelWeights labelWeights;

std::string labelWeightsStr = vm["label-weights"].as<std::string>();
if (!labelWeightsStr.empty())
{
    labelWeights = LabelWeightsParser::ParseVector(labelWeightsStr);
    std::cout << "[Label Weights] Parsed vector from --label-weights" << std::endl;
}
else if (vm["labelkappa"].as<double>() > 1e-6 || !vm["fixedlabelmap"].as<std::string>().empty())
{
    labelWeights = LabelWeightsParser::ScalarLabelWeights(
        vm["labelkappa"].as<double>(),
        vm["labelkappaderiv"].as<double>()
    );
    
    // Auto-expand if labelmap is provided
    std::string fixedLabelMapPath = vm["fixedlabelmap"].as<std::string>();
    if (fixedLabelMapPath != "N" && !fixedLabelMapPath.empty())
    {
        unsigned int numLabels = LabelWeightsParser::DetectNumberOfLabels(fixedLabelMapPath);
        if (numLabels > 0)
        {
            labelWeights = LabelWeightsParser::ExpandToVector(labelWeights, numLabels);
        }
    }
}
else
{
    labelWeights = LabelWeightsParser::Disabled();
}

std::cout << "[Label Weights Configuration]" << std::endl;
labelWeights.Print("  ");
```

---

## Step 4: Apply to Metric Setup (Same Pattern for All 4)

**Location:** Replace the individual `metric->SetAlpha()`, `metric->SetLambda()`, etc. calls

**OLD (current):**
```cpp
metric->SetAlpha(ALPHA);
metric->SetAlphaDerivative(ALPHADERIVATIVE);
metric->SetLambda(LAMBDA);
metric->SetLambdaDerivative(LAMBDADERIVATIVE);
metric->SetNu(NU);
metric->SetNuDerivative(NUDERIVATIVE);
// ... etc
```

**NEW:**
```cpp
// Apply main metrics from parsed config
metric->SetAlpha(metricsConfig.mi.weight);
metric->SetAlphaDerivative(metricsConfig.mi.derivative);
metric->SetLambda(metricsConfig.ngf.weight);
metric->SetLambdaDerivative(metricsConfig.ngf.derivative);
metric->SetNu(metricsConfig.mse.weight);
metric->SetNuDerivative(metricsConfig.mse.derivative);
metric->SetRho(metricsConfig.gd.weight);
metric->SetRhoDerivative(metricsConfig.gd.derivative);
metric->SetYota(metricsConfig.nc.weight);
metric->SetYotaDerivative(metricsConfig.nc.derivative);
metric->SetSigma(metricsConfig.nmi.weight);
metric->SetSigmaDerivative(metricsConfig.nmi.derivative);

// Apply label weights
if (labelWeights.IsEnabled())
{
    metric->SetLabelKappa(labelWeights.GetScalarKappa());
    // If your Mplus supports vector label weights:
    // for (size_t i = 0; i < labelWeights.kappaValues.size(); ++i)
    //     metric->SetLabelKappa(i, labelWeights.kappaValues[i]);
}
```

---

## Step 5: Update Sampling Percentages

**Location:** Replace individual sampling percentage assignments

**OLD:**
```cpp
metric->SetMANumberOfSamples(static_cast<unsigned int>(numberOfPixels * MAPERCENTAGE));
metric->SetNGFNumberOfSamples(static_cast<unsigned int>(numberOfPixels * NGFPERCENTAGE));
metric->SetMSENumberOfSamples(static_cast<unsigned int>(numberOfPixels * MSEPERCENTAGE));
// ... etc
```

**NEW:**
```cpp
metric->SetMANumberOfSamples(static_cast<unsigned int>(numberOfPixels * metricsConfig.mi.samplingPercent));
metric->SetNGFNumberOfSamples(static_cast<unsigned int>(numberOfPixels * metricsConfig.ngf.samplingPercent));
metric->SetMSENumberOfSamples(static_cast<unsigned int>(numberOfPixels * metricsConfig.mse.samplingPercent));
metric->SetGDNumberOfSamples(static_cast<unsigned int>(numberOfPixels * metricsConfig.gd.samplingPercent));
metric->SetNCNumberOfSamples(static_cast<unsigned int>(numberOfPixels * metricsConfig.nc.samplingPercent));
metric->SetNMINumberOfSamples(static_cast<unsigned int>(numberOfPixels * metricsConfig.nmi.samplingPercent));
```

---

## Summary of Changes Per File

| File | Includes | CLI Options | Parsing Logic | Metric Setup | Sampling |
|------|----------|-----------|---------------|-------------|----------|
| 3DRegAffine.cxx | ✅ | ✅ | ✅ | ✅ | ✅ |
| 3DRegAffineMultiLevel.cxx | ✅ | ✅ | ✅ | ✅ | ✅ |
| 3DRegBsplines.cxx | ✅ | ✅ | ✅ | ✅ | ✅ |
| 3DRegSimilarity.cxx | ✅ | ✅ | ✅ | ✅ | ✅ |

---

## Testing

After updating, test each executable:

```bash
# Test new preset format
./3DRegAffine --fixed A.nii --moving B.nii --output T.txt --preset multimodal

# Test array format
./3DRegAffine --fixed A.nii --moving B.nii --output T.txt \
  --metrics "1.0,0.5,0,0,0,0,0" \
  --metric-sampling "0.1,0.1,0.1,0.1,0.1,0.1,0.1"

# Test backward compatibility (individual params)
./3DRegAffine --fixed A.nii --moving B.nii --output T.txt \
  --alpha 1.0 --lambda 0.5 --nu 0

# Test conflict warning
./3DRegAffine --fixed A.nii --moving B.nii --output T.txt \
  --metrics "1.0,0.5,0,0,0,0,0" \
  --alpha 2.0  # Should warn
```
