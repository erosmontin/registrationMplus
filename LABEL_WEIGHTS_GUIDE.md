# Label Weights CLI Integration Guide

## Overview

`LabelWeightsParser.h` handles both scalar and vector label weights for multi-label registration. It automatically detects which you're using and expands scalars to vectors as needed.

---

## Usage Patterns

### **Pattern 1: Scalar Label Weight (Simplest)**

```bash
# Single weight applied to ALL labels
./3DRegAffine --fixed fixed.nii --moving moving.nii --output transform.txt \
  --fixed-labelmap labels.nii \
  --labelkappa 0.5
```

**What happens:**
1. CLI parses `--labelkappa 0.5` → scalar weight
2. LabelMap is loaded (3 unique structures detected)
3. Scalar expands → `[0.5, 0.5, 0.5]` (one per structure)
4. All structures equally weighted

**Pros:**
- User doesn't need to know structure count
- Simplest for beginners
- Works if structures are equally important

---

### **Pattern 2: Per-Label Weights (Advanced)**

```bash
# Different weight per structure (e.g., hippocampus, ventricles, tumor)
./3DRegAffine --fixed fixed.nii --moving moving.nii --output transform.txt \
  --fixed-labelmap labels.nii \
  --label-weights "0.5,0.3,0.2"
```

**What happens:**
1. CLI parses `--label-weights "0.5,0.3,0.2"` → vector
2. Weights applied to labels 1, 2, 3 respectively
3. Label 1 (hippocampus): 0.5 weight
4. Label 2 (ventricles): 0.3 weight
5. Label 3 (other): 0.2 weight

**Pros:**
- Fine-grained control
- Can emphasize important structures

**Cons:**
- User must know label count and order
- More explicit

---

### **Pattern 3: Disable Labels (Default)**

```bash
# No label map penalty
./3DRegAffine --fixed fixed.nii --moving moving.nii --output transform.txt
```

**What happens:**
- `labelkappa` defaults to 0 → labels disabled
- No label penalty applied

---

## CLI Options

Add these to your `po::options_description`:

```cpp
("fixed-labelmap", po::value<std::string>()->default_value("N"), 
 "Fixed image label map (N or empty = none)")
("moving-labelmap", po::value<std::string>()->default_value("N"), 
 "Moving image label map (N or empty = none)")

("labelkappa", po::value<double>()->default_value(0.0), 
 "Scalar label weight (applied to ALL structures)")
("label-weights", po::value<std::string>()->default_value(""), 
 "Per-label weights: comma-separated list (e.g., '0.5,0.3,0.2')")

("labelkappadervative", po::value<double>()->default_value(-1), 
 "Label derivative (default: same as labelkappa, or 0 if disabled)")
("label-derivatives", po::value<std::string>()->default_value(""), 
 "Per-label derivatives (default: auto-derived from weights)")
```

---

## Integration Code

After parsing command-line arguments:

```cpp
#include "../LabelWeightsParser.h"

// Check if user provided label weights
std::string labelWeightsStr = vm["label-weights"].as<std::string>();
std::string labelDerivativesStr = vm["label-derivatives"].as<std::string>();
double scalarKappa = vm["labelkappa"].as<double>();
double scalarKappaDeriv = vm["labelkappadervative"].as<double>();

// Parse label weights
LabelWeightsParser::LabelWeights labelWeights;

if (!labelWeightsStr.empty())
{
    // User provided explicit per-label weights
    labelWeights = LabelWeightsParser::ParseVector(labelWeightsStr);
    std::cout << "Parsed label weights array" << std::endl;
}
else if (scalarKappa > 1e-6 || vm["fixed-labelmap"].as<std::string>() != "N")
{
    // User provided scalar weight, or labelmap without explicit weights
    labelWeights = LabelWeightsParser::ScalarLabelWeights(scalarKappa, scalarKappaDeriv);
    std::cout << "Using scalar label weight" << std::endl;

    // If labelmap is loaded, detect label count and expand scalar → vector
    std::string fixedLabelMapPath = vm["fixed-labelmap"].as<std::string>();
    if (fixedLabelMapPath != "N" && !fixedLabelMapPath.empty())
    {
        unsigned int numLabels = LabelWeightsParser::DetectNumberOfLabels(fixedLabelMapPath);
        if (numLabels > 0)
        {
            labelWeights = LabelWeightsParser::ExpandToVector(labelWeights, numLabels);
            std::cout << "Expanded scalar weight to " << numLabels << " labels" << std::endl;
        }
    }
}
else
{
    // No labels
    labelWeights = LabelWeightsParser::Disabled();
    std::cout << "Label weighting disabled" << std::endl;
}

// Print for user feedback
labelWeights.Print("  ");

// Apply to metric
if (labelWeights.IsEnabled())
{
    // Set per-label weights
    for (size_t i = 0; i < labelWeights.kappaValues.size(); ++i)
    {
        // Your metric may need a SetLabelKappa(label_index, weight) method
        // or you can store as vector and pass during Initialize()
        metric->SetLabelKappa(labelWeights.kappaValues[i]);
        metric->SetLabelKappaDerivative(labelWeights.kappaDerivatives[i]);
    }
}
```

---

## Real Example: Brain ROI Registration

```bash
# Register brain with 3 structures: hippocampus, ventricles, thalamus
# Prioritize hippocampus > ventricles > thalamus

./3DRegAffine \
  --fixed brain_template.nii \
  --moving patient_brain.nii \
  --fixed-labelmap template_rois.nii \
  --output patient_to_template.txt \
  --preset multimodal \
  --label-weights "0.8,0.5,0.3" \
  --threads 8
```

**Configuration:**
- Hippocampus (label 1): 0.8 weight — most important
- Ventricles (label 2): 0.5 weight — moderately important  
- Thalamus (label 3): 0.3 weight — less important
- Mi + NGF still active (`--preset multimodal`)
- Total objective: MI + NGF + weighted(ROI penalties)

---

## FAQ

**Q: How do I know my label map structure order?**
A: Use an image viewer (ITK-SNAP, 3D Slicer) or run:
```python
import nibabel as nib
import numpy as np
img = nib.load('labels.nii.gz')
unique = np.unique(img.get_fdata())
for label in unique:
    if label > 0:
        print(f"Label {int(label)}")
```

**Q: What if I have 10 labels but only care about 3?**
A: Set zero weights for the others:
```bash
--label-weights "0.5,0.3,0,0,0,0,0,0,0,0"
```

**Q: Can I mix scalar and vector weights?**
A: No. Use either:
- `--labelkappa 0.5` (scalar, gets expanded)
- `--label-weights "0.5,0.3,0.2"` (vector, explicit)

**Q: What if label count doesn't match my vector?**
A: Parser will error with clear message:
```
ERROR: Label weights vector (3 values) doesn't match label count (5 detected).
       Use --label-weights "0.5,0.3,0.2,0,0" to specify all 5 labels.
```

---

## Implementation Notes

### Detecting Label Count

`LabelWeightsParser::DetectNumberOfLabels()` is a stub. To implement properly:

```cpp
// In your imageUtils or similar
#include "itkImage.h"
#include "itkImageFileReader.h"

template <typename LabelImageType>
unsigned int CountUniqueLabels(const std::string& path)
{
    auto reader = itk::ImageFileReader<LabelImageType>::New();
    reader->SetFileName(path);
    reader->Update();

    auto img = reader->GetOutput();
    std::set<typename LabelImageType::PixelType> unique;

    itk::ImageRegionConstIterator<LabelImageType> it(img, img->GetLargestPossibleRegion());
    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
        if (it.Get() > 0)
            unique.insert(it.Get());

    return unique.size();
}
```

Then call from `LabelWeightsParser::DetectNumberOfLabels()`.

### Vector Storage in Metric

Your `Mplus` metric may need to support vector label weights. Options:

**Option A:** Store as vector member
```cpp
class Mplus {
    std::vector<double> m_LabelKappas;
    void SetLabelKappas(const std::vector<double>& weights);
};
```

**Option B:** Call SetLabelKappa multiple times
```cpp
for (double w : labelWeights.kappaValues)
    metric->SetLabelKappa(w);
```

**Option C:** Expand in label metric, not in registration driver
