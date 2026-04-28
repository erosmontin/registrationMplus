# 3DRegBsplines - ALL CONFIGURATIONS COMPLETE

## 📊 WHERE ARE ALL THE CONFIGURATIONS?

Here's the COMPLETE list of all 70+ CLI parameters organized by category in the source code.

**File:** `src/3DRegistration/3DRegBsplines/src/3DRegBsplines.cxx`  
**Lines:** 83-238 (CLI definitions) + 400-600+ (parsing and usage)

---

## 🎯 CATEGORY 1: REQUIRED INPUTS (No Defaults)

```cpp
// Lines: 84-86
--fixedimage,f          Required: Fixed image filename
--movingimage,m         Required: Moving image filename
--outputimage,o         Required: Output registered image filename
```

**In code:**
```cpp
std::string fixedImageFN = vm["fixedimage"].as<std::string>();
std::string movingImageFN = vm["movingimage"].as<std::string>();
std::string ou = vm["outputimage"].as<std::string>();
```

---

## 📤 CATEGORY 2: OUTPUT OPTIONS (Default: "N" = off)

```cpp
// Lines: 87-88, 115-118
--vfout,v               Default: "N"  → Deformation field output
--transformout,T        Default: "N"  → Transform (reusable)
--transformin,W         Default: "N"  → Input transform
--gridposition,G        Default: "N"  → Read grid position from image
```

**In code:**
```cpp
std::string VOUT = vm["vfout"].as<std::string>();
std::string TOUT = vm["transformout"].as<std::string>();
std::string TIN = vm["transformin"].as<std::string>();
std::string GRIDPOSITION = vm["gridposition"].as<std::string>();

if (GRIDPOSITION != "N") {
    // read the image that specifies the position of the grid
    FixedImageReaderType::Pointer meshImageReader = FixedImageReaderType::New();
    meshImageReader->SetFileName(GRIDPOSITION);
    meshImageReader->Update();
    // ... override mesh parameters from image
}
```

---

## ⚙️ CATEGORY 3: SYSTEM & THREADING

```cpp
// Lines: 88, 132
--numberofthreads       Default: 2    → CPU threads for computation
--verbose,V             Default: false → Print verbose output
```

**In code:**
```cpp
int NT = vm["numberofthreads"].as<int>();
bool V = vm["verbose"].as<bool>();

registration->SetNumberOfThreads(NT);

if (vm["verbose"].as<bool>())
    RegCommon::PrintOptions(vm);
```

---

## 🎯 CATEGORY 4: METRIC WEIGHTS (PRIMARY - NEW DEFAULTS!)

```cpp
// Lines: 90-93, 105-110
--alpha,a               Default: 1.0   ← Mutual Information (MI)
--lambda,l              Default: 0.5   ⭐ Normalized Gradient Field (NGF) NEW!
--nu,n                  Default: 0.0   ← Sum Squared Differences (MSE)
--rho                   Default: 0.0   ← Gradient Difference (GD)
--yota,y                Default: 0.0   ← Normalized Correlation (NC)
--sigma                 Default: 0.0   ← Normalized Mutual Information (NMI)
```

**In code:**
```cpp
double ALPHA = vm["alpha"].as<double>();
double LAMBDA = vm["lambda"].as<double>();
double NU = vm["nu"].as<double>();
double RHO = vm["rho"].as<double>();
double YOTA = 0.0;  // Special: references external variable
double SIGMA = vm["sigma"].as<double>();

// Used by metric
metric->SetAlpha(ALPHA);
metric->SetLambda(LAMBDA);
metric->SetNu(NU);
metric->SetRho(RHO);
metric->SetYota(YOTA);
metric->SetSigma(SIGMA);
```

---

## 📈 CATEGORY 5: METRIC DERIVATIVES

```cpp
// Lines: 91-92, 107-109, 148-153
--alphaderivative,A     Default: 1.0   ← MI derivative
--lambdaderivative,L    Default: 0     ← NGF derivative
--nuderivative,N        Default: 0     ← MSE derivative
--rhoderivative         Default: 0.0   ← GD derivative
--yotaderivative,Y      Default: 0     ← NC derivative
--sigmaderivative       Default: 0.0   ← NMI derivative
```

**In code:**
```cpp
double ALPHADERIVATIVE = vm["alphaderivative"].as<double>();
double LAMBDADERIVATIVE = vm["lambdaderivative"].as<double>();
double NUDERIVATIVE = vm["nuderivative"].as<double>();
double RHODERIVATIVE = vm["rhoderivative"].as<double>();
double YOTADERIVATIVE = 0.0;  // Special
double SIGMADERIVATIVE = vm["sigmaderivative"].as<double>();

metric->SetAlphaDerivative(ALPHADERIVATIVE);
metric->SetLambdaDerivative(LAMBDADERIVATIVE);
// ... etc
```

---

## 🔍 CATEGORY 6: MUTUAL INFORMATION (MI/MATTES)

```cpp
// Lines: 90-92
--mattespercentage,p    Default: 0.1   → Use 10% of image pixels
--mattesnumberofbins,b  Default: 64    → MI histogram bins (64-256)
--bsplinecaching,B      Default: true  → Cache B-spline weights
--explicitPDFderivatives Default: false → Use implicit vs explicit derivatives
```

**In code:**
```cpp
int NB = vm["mattesnumberofbins"].as<int>();
double MAPERCENTAGE = vm["mattespercentage"].as<double>();
bool TB = vm["bsplinecaching"].as<bool>();
bool EPDF = vm["explicitPDFderivatives"].as<bool>();

metric->SetNumberOfHistogramBins(NB);
metric->SetFixedImageSamplesPercentage(MAPERCENTAGE);
metric->SetUseExplicitPDFDerivatives(EPDF);

const unsigned int numberOfSamplesMA = 
    static_cast<unsigned int>(numberOfPixels * MAPERCENTAGE);
```

---

## 🌊 CATEGORY 7: NORMALIZED GRADIENT FIELD (NGF)

```cpp
// Lines: 94-99
--etavaluefixed,r       Default: -1    → Fixed image noise (-1 = auto-detect)
--etavaluemoving,s      Default: -1    → Moving image noise (-1 = auto-detect)
--NGFevaluator          Default: 0     → NGF type: 0=scalar, 1=cross, 2=scdelta, 3=Delta, 4=Delta2
--ngfprecompute         Default: false → Precompute NGF (faster, approximate)
--ngfpercentage         Default: 0.1   → Use 10% of pixels for NGF
--ngfspacing            Default: "4,4,4" → Gradient spacing (x,y,z in mm)
```

**In code:**
```cpp
double ETAF = vm["etavaluefixed"].as<double>();
double ETAM = vm["etavaluemoving"].as<double>();
int NGFevaluator = vm["NGFevaluator"].as<int>();

if ((LAMBDA != 0) || (LAMBDADERIVATIVE != 0)) {
    if ((ETAF == -1) || (ETAM == -1)) {
        metric->SetAutoEstimateEta(true);
    }
}

// Parse NGF spacing
auto s = vm["ngfspacing"].as<std::string>();
std::replace(s.begin(), s.end(), ',', ' ');
std::istringstream iss(s);
std::vector<double> tmp((std::istream_iterator<double>(iss)),
                        std::istream_iterator<double>());
ImageType::SpacingType ngf;
for (unsigned i = 0; i < ImageDimension; ++i)
    ngf[i] = tmp[i];
// Set on metric...
```

---

## 📊 CATEGORY 8: MSE, GD, NC, NMI SAMPLING

```cpp
// Lines: 120-128
--msepercentage         Default: 0.1   → Use 10% of pixels for MSE
--normalizemse          Default: false → Normalize MSE by intensity range (mean in [0,1])
--normalizegd           Default: false → Normalize GD by overlap voxel count (mean in [0,1])
--gdpercentage          Default: 0.1   → Use 10% of pixels for GD
--ncpercentage          Default: 0.1   → Use 10% of pixels for NC
--nmipercentage         Default: 0.1   → Use 10% of pixels for NMI
--nmibins               Default: 64    → NMI histogram bins
```

**In code:**
```cpp
double mseSamplingPercent = metricsConfig.mse.samplingPercent;  // 0.1 default
double gdSamplingPercent = metricsConfig.gd.samplingPercent;    // 0.1 default
double ncSamplingPercent = metricsConfig.nc.samplingPercent;    // 0.1 default
double nmiSamplingPercent = metricsConfig.nmi.samplingPercent;  // 0.1 default
int NMIBINS = vm["nmibins"].as<int>();
```

---

## 🎛️ CATEGORY 9: B-SPLINE MESH & OPTIMIZATION (NEW DEFAULTS!)

```cpp
// Lines: 100-104, 130
--gridresolution,g      Default: 50    → B-spline mesh spacing (mm)
--maxnumberofiterations,I Default: 1000 → Max optimization iterations
--costfunctionconvergencefactor,F Default: 1.e7 ⭐ Convergence (NEW!)
--projectedgradienttolerance,P Default: 1.e-5 → Gradient tolerance
--numberofevaluations,E Default: 500   → L-BFGS-B evaluations
--numberofcorrections,C Default: 5     → L-BFGS-B correction pairs
--overlappadding        Default: 5     ⭐ Control points (NEW!)
--meshmarginsize        Default: 0.0   → Mesh domain extension (mm)
```

**In code:**
```cpp
double GRIDRESOLUTION = vm["gridresolution"].as<double>();
int NI = vm["maxnumberofiterations"].as<int>();
double CFCF = vm["costfunctionconvergencefactor"].as<double>();
double PGT = vm["projectedgradienttolerance"].as<double>();
int NE = vm["numberofevaluations"].as<int>();
int NC = vm["numberofcorrections"].as<int>();
double meshMargin = vm["meshmarginsize"].as<double>();
unsigned int borderNodesPerSide = vm["overlappadding"].as<unsigned int>();

// Applied to transform
for (unsigned int i = 0; i < SpaceDimension; ++i) {
    const double extension = borderNodesPerSide * GRIDRESOLUTION;
    fixedOrigin[i] = meshorigin[i] - dirSign * (meshMargin + extension);
    fixedPhysicalDimensions[i] = 
        meshspacing[i] * (meshsize[i] - 1) + 2.0 * (meshMargin + extension);
    meshSize[i] = staticc_cast<unsigned int>(
        fixedPhysicalDimensions[i] / GRIDRESOLUTION) - SplineOrder;
}

transform->SetTransformDomainOrigin(fixedOrigin);
transform->SetTransformDomainPhysicalDimensions(fixedPhysicalDimensions);
transform->SetTransformDomainMeshSize(meshSize);

// Applied to optimizer
optimizer->SetCostFunctionConvergenceFactor(CFCF);
optimizer->SetProjectedGradientTolerance(PGT);
optimizer->SetMaximumNumberOfIterations(NI);
optimizer->SetMaximumNumberOfEvaluations(NE);
optimizer->SetMaximumNumberOfCorrections(NC);
```

---

## ⚖️ CATEGORY 10: BOUNDARY CONDITIONS & IMAGE THRESHOLDING

```cpp
// Lines: 111-113, 117
--fixedimagethreshold,t Default: -99999999 → Ignore pixels below this
--dfltpixelvalue        Default: 0    → Default for resampled pixels
--bound                 Default: 0    → Boundary: 0=unbounded, 1=lower, 2=both, 3=upper
--lbound                Default: 0    → Lower bound for parameters
--ubound                Default: 0    → Upper bound for parameters
```

**In code:**
```cpp
double TR = vm["fixedimagethreshold"].as<double>();
double DFLTPIXELVALUE = vm["dfltpixelvalue"].as<double>();
int BOUND = vm["bound"].as<int>();
double LBOUND = vm["lbound"].as<double>();
double UBOUND = vm["ubound"].as<double>();

// Used by metric
metric->SetFixedImageThreshold(TR);

// Used by resampler
resampler->SetDefaultPixelValue(DFLTPIXELVALUE);

// Used by optimizer for parameter constraints
optimizer->SetBound(BOUND, LBOUND, UBOUND);
```

---

## 🔀 CATEGORY 11: MULTI-METRIC CONFIGURATION

```cpp
// Lines: 129-130
--derivativemode        Default: 0    → Mode: 0=consistent, 2=adaptive (1 not allowed)
--mainmetric            Default: 0    → Main metric for mode 2: 0=MI, 1=NGF, 2=MSE, 3=NC, 4=Label, 5=GD, 6=NMI
--metricoverlap         Default: true → Compute overlap between images
```

**In code:**
```cpp
int DERIVMODE = vm["derivativemode"].as<int>();
int MAINMETRIC = vm["mainmetric"].as<int>();
bool METRICOVERLAP = vm["metricoverlap"].as<bool>();

if (DERIVMODE == 1) {
    std::cerr << "ERROR: derivativemode=1 (normalized) is not compatible with "
                 "LBFGS-B optimizer used by B-splines. Use 0 or 2." << std::endl;
    return EXIT_FAILURE;
}

metric->SetDerivativeMode(DERIVMODE);
metric->SetMainMetric(MAINMETRIC);
metric->SetComputeMetricOverlap(METRICOVERLAP);
```

---

## 🏷️ CATEGORY 12: LABEL/SEGMENTATION SUPPORT

```cpp
// Lines: 167-172
--fixedlabelmap         Default: "N"  → Fixed segmentation (N=none)
--movinglabelmap        Default: "N"  → Moving segmentation (N=none)
--labelkappa            Default: 0.0  → Label metric weight (0=off)
--labelkappaderiv       Default: 0.0  → Label derivative weight
--labelkappavec         Default: ""   → Per-label weights (L1:w1,L2:w2,...)
--labelkappaderivvec    Default: ""   → Per-label derivatives
--labelsamples          Default: 0.1  → Fraction (0,1] of voxels OR absolute count (>1)
--labelreport           Default: 1    → Report Dice every N iterations (0=off)
```

**In code:**
```cpp
const std::string FIXEDLABELMAP = vm["fixedlabelmap"].as<std::string>();
const std::string MOVINGLABELMAP = vm["movinglabelmap"].as<std::string>();
const double LABELKAPPA = vm["labelkappa"].as<double>();
const double LABELKAPPADERIV = vm["labelkappaderiv"].as<double>();
const double LABELSAMPLES = vm["labelsamples"].as<double>();
const int LABELREPORT = vm["labelreport"].as<int>();

const auto LABELKAPPAVEC = RegCommon::ParseLabelWeights(
    vm["labelkappavec"].as<std::string>());
const auto LABELKAPPADERIVVEC = RegCommon::ParseLabelWeights(
    vm["labelkappaderivvec"].as<std::string>());

// Read label maps
if (FIXEDLABELMAP != "N" && MOVINGLABELMAP != "N") {
    typedef itk::ImageFileReader<LabelImageType> LabelReaderType;
    auto flr = LabelReaderType::New(); 
    flr->SetFileName(FIXEDLABELMAP); 
    flr->Update();
    LabelImageType::ConstPointer fixedLabelMap = flr->GetOutput();
    
    auto mlr = LabelReaderType::New();
    mlr->SetFileName(MOVINGLABELMAP);
    mlr->Update();
    LabelImageType::ConstPointer movingLabelMap = mlr->GetOutput();
}

metric->SetLabelKappa(LABELKAPPA);
metric->SetLabelKappaDerivative(LABELKAPPADERIV);
metric->SetLabelSamplingPercent(LABELSAMPLES);
```

---

## 📸 CATEGORY 13: VISUALIZATION & SNAPSHOTS

```cpp
// Lines: 177-181
--snapshotdir           Default: "N"  → Directory for snapshots (N=off)
--snapshotevery         Default: 1    → Save PNG every N iterations
--snapshotstack         Default: false → Full 3D .nii.gz vs mid-slice PNG
--snapshotgrid          Default: false → Show grid vs knot mesh
--snapshotgridspacing   Default: 20   → Grid line spacing (voxels)
```

**In code:**
```cpp
const std::string SNAPSHOTDIR = vm["snapshotdir"].as<std::string>();
const int SNAPSHOTEVERY = vm["snapshotevery"].as<int>();
const bool SNAPSHOTSTACK = vm["snapshotstack"].as<bool>();
const bool SNAPSHOTGRID = vm["snapshotgrid"].as<bool>();
const unsigned int SNAPSHOTGRIDSP = vm["snapshotgridspacing"].as<unsigned int>();

if (SNAPSHOTDIR != "N") {
    // Create snapshots directory
    // Save PNG every N iterations during optimization
    // Overlay grid or knot mesh
}
```

---

## 🎁 CATEGORY 14: PRESET & SIMPLIFIED CLI

```cpp
// Lines: 191-206 (NEW - SIMPLIFIED INTERFACE)
--preset                Default: ""   → Preset: 'multimodal', 'singlemodal', 'rigid'
--metrics               Default: ""   → Metric array: alpha,lambda,nu,rho,yota,sigma
--metric-derivatives    Default: ""   → Derivative array: alpha_d,lambda_d,nu_d,rho_d,yota_d,sigma_d
--metric-sampling       Default: ""   → Sampling array: ma%,ngf%,mse%,gd%,nc%,nmi% (label → --labelsamples)
--label-weights         Default: ""   → Per-label weights (alternative to --labelkappa)
--label-derivatives     Default: ""   → Per-label derivatives
--modality              Default: "custom" → 'multimodal', 'singlemodal', or 'custom'
```

**In code (Preset handling):**
```cpp
bool hasMetricsArray = (!vm["metrics"].as<std::string>().empty() ||
                        !vm["metric-derivatives"].as<std::string>().empty() ||
                        !vm["preset"].as<std::string>().empty());

MetricsConfig::MainMetricsConfig metricsConfig;

if (!vm["preset"].as<std::string>().empty()) {
    metricsConfig = MetricsConfig::GetMainPreset(vm["preset"].as<std::string>());
    std::cout << "\n[CLI] Using preset: " << vm["preset"].as<std::string>() << std::endl;
}
else if (!vm["metrics"].as<std::string>().empty()) {
    metricsConfig = MetricsConfig::ParseMainWeights(vm["metrics"].as<std::string>());
}

// Apply derivatives if provided
if (!vm["metric-derivatives"].as<std::string>().empty()) {
    metricsConfig = MetricsConfig::ParseMainDerivatives(metricsConfig, ...);
}

// Apply sampling if provided
if (!vm["metric-sampling"].as<std::string>().empty()) {
    metricsConfig = MetricsConfig::ParseMainSampling(metricsConfig, ...);
}

// Modality presets
const std::string MODALITY = vm["modality"].as<std::string>();
if (MODALITY == "multimodal") {
    // MI(1.0) + NGF(0.5) - good for different modalities
    alpha = 1.0; lambda = 0.5; nu = 0.0; yota = 0.0;
} 
else if (MODALITY == "singlemodal") {
    // MSE(1.0) + NC(0.5) - good for same modality
    alpha = 0.0; lambda = 0.0; nu = 1.0; yota = 0.5;
}
```

---

## 📋 COMPLETE PARAMETER REFERENCE TABLE

| Category | Parameter | Type | Default | Lines |
|----------|-----------|------|---------|-------|
| **Input** | fixedimage | string | - | 84 |
| | movingimage | string | - | 85 |
| | outputimage | string | - | 86 |
| **Output** | vfout | string | "N" | 87 |
| | transformout | string | "N" | 115 |
| | transformin | string | "N" | 116 |
| | gridposition | string | "N" | 117 |
| **System** | numberofthreads | int | 2 | 88 |
| | verbose | bool | false | 119 |
| **Metrics** | alpha | double | 1.0 | 90 |
| | lambda | double | **0.5** ⭐ | 93 |
| | nu | double | 0.0 | 100 |
| | rho | double | 0.0 | 125 |
| | yota | double | 0.0 | 127 |
| | sigma | double | 0.0 | 126 |
| **Derivatives** | alphaderivative | double | 1.0 | 91 |
| | lambdaderivative | double | 0 | 94 |
| | nuderivative | double | 0 | 101 |
| | rhoderivative | double | 0.0 | 148 |
| | yotaderivative | double | 0 | 149 |
| | sigmaderivative | double | 0.0 | 150 |
| **MI/Mattes** | mattespercentage | double | 0.1 | 90 |
| | mattesnumberofbins | int | 64 | 91 |
| | bsplinecaching | bool | true | 92 |
| | explicitPDFderivatives | bool | false | 93 |
| **NGF** | etavaluefixed | double | -1 | 96 |
| | etavaluemoving | double | -1 | 97 |
| | NGFevaluator | int | 0 | 98 |
| | ngfprecompute | bool | false | 99 |
| | ngfpercentage | double | 0.1 | 120 |
| | ngfspacing | string | "4,4,4" | 153 |
| **Sampling** | msepercentage | double | 0.1 | 121 |
| | gdpercentage | double | 0.1 | 123 |
| | ncpercentage | double | 0.1 | 125 |
| | nmipercentage | double | 0.1 | 124 |
| | nmibins | int | 64 | 151 |
| | normalizemse | bool | false | 122 |
| | normalizegd | bool | false | 127 |
| **Optimization** | gridresolution | double | 50 | 102 |
| | maxnumberofiterations | int | 1000 | 103 |
| | costfunctionconvergencefactor | double | **1.e7** ⭐ | 104 |
| | projectedgradienttolerance | double | 1.e-5 | 106 |
| | numberofevaluations | int | 500 | 107 |
| | numberofcorrections | int | 5 | 108 |
| | overlappadding | uint | **5** ⭐ | 186 |
| | meshmarginsize | double | 0.0 | 154 |
| **Boundary** | fixedimagethreshold | double | -99999999 | 109 |
| | dfltpixelvalue | double | 0 | 117 |
| | bound | int | 0 | 110 |
| | lbound | double | 0 | 112 |
| | ubound | double | 0 | 113 |
| **Multi-Metric** | derivativemode | int | 0 | 151 |
| | mainmetric | int | 0 | 152 |
| | metricoverlap | bool | true | 165 |
| **Labels** | fixedlabelmap | string | "N" | 166 |
| | movinglabelmap | string | "N" | 167 |
| | labelkappa | double | 0.0 | 168 |
| | labelkappaderiv | double | 0.0 | 169 |
| | labelkappavec | string | "" | 170 |
| | labelkappaderivvec | string | "" | 171 |
| | labelsamples | double | 0.1 | 172 |
| | labelreport | int | 1 | 173 |
| **Snapshots** | snapshotdir | string | "N" | 174 |
| | snapshotevery | int | 1 | 175 |
| | snapshotstack | bool | false | 176 |
| | snapshotgrid | bool | false | 177 |
| | snapshotgridspacing | uint | 20 | 178 |
| **Preset** | preset | string | "" | 194 |
| | metrics | string | "" | 197 |
| | metric-derivatives | string | "" | 200 |
| | metric-sampling | string | "" | 203 |
| | label-weights | string | "" | 206 |
| | label-derivatives | string | "" | 209 |
| | modality | string | "custom" | 188 |

---

## 🔗 HOW PARAMETERS FLOW

```
CLI Input
    ↓
Parse with boost::program_options (lines 83-210)
    ↓
Extract to variables (lines 400-440)
    ↓
Apply to metric, optimizer, transform (lines 450+)
    ↓
Execute registration
    ↓
Output files
```

---

## ✅ ALL NEW DEFAULTS (After Fixes)

| Parameter | Old | New | Line |
|-----------|-----|-----|------|
| `--lambda` | 0 | **0.5** | 93 |
| `--costfunctionconvergencefactor` | 1.e12 | **1.e7** | 104 |
| `--overlappadding` | 1 | **5** | 186 |

---

## 🎯 QUICK SUMMARY

**Total CLI Parameters:** 70+

**Breakdown:**
- Required inputs: 3
- Metric weights: 6
- Metric derivatives: 6
- Sampling controls: 7
- Optimization: 8
- Mesh: 4
- Labels: 8
- Snapshots: 5
- Presets/simplified: 7
- Other: 10

All defined in lines **83-210** of 3DRegBsplines.cxx  
All applied in lines **400+** of 3DRegBsplines.cxx
