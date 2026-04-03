/**
 * INTEGRATION GUIDE: New Simplified CLI Format
 * 
 * This document shows how to integrate the new CliParser into 3DRegAffine.cxx
 * with automatic format detection. Changes are minimal and backward-compatible.
 */

// ============================================================================
// STEP 1: Add include at the top of 3DRegAffine.cxx
// ============================================================================

#include "../CliParser.h"  // Add this near other includes


// ============================================================================
// STEP 2: Add new simplified options to po::options_description (around line 85)
// ============================================================================

// === REQUIRED: IO and Transform ===
("fixed,F", po::value<std::string>()->required(), "Fixed image filename")
("moving,M", po::value<std::string>()->required(), "Moving image filename")
("output,O", po::value<std::string>()->required(), "Output transform filename")

// === NEW FORMAT: Simplified Metrics (optional, auto-detected) ===
("preset", po::value<std::string>()->default_value(""), 
 "Metric preset: 'multimodal' (alpha=1,lambda=0.5), 'singlemodal' (nu=1,yota=0.5), 'rigid', or 'custom' (no defaults)")
("weights", po::value<std::string>()->default_value(""), 
 "Metric weights as comma-separated array: alpha,lambda,nu,rho,yota,kappa,sigma (e.g., '1.0,0.5,0,0,0,0,0')")
("metric-percentages", po::value<std::string>()->default_value(""), 
 "Sampling percentages as comma-separated array: ma,ngf,mse,gd,nc,nmi,label (default: all 0.1)")

// === OLD FORMAT: Kept for backward compatibility ===
// (Keep all existing --alpha, --lambda, --nu, etc. options)
// ... existing code ...

// === OPTIMIZER ===
("iterations,I", po::value<int>()->default_value(1000), "Max iterations")
("step-length,S", po::value<double>()->default_value(0.1), "Minimum step length")
("tolerance,G", po::value<double>()->default_value(1e-4), "Gradient magnitude tolerance")

// === OPTIONAL ADVANCED ===
("normalizemse", po::value<bool>()->default_value(false), "Normalize MSE by intensity-range^2")
("derivativemode", po::value<int>()->default_value(0), "0=consistent, 1=normalized (RSGD), 2=adaptive")
("threads,T", po::value<unsigned int>()->default_value(1), "Number of threads")
("verbose,V", po::value<bool>()->default_value(false), "Verbose output")

// REMOVED (have sensible defaults now):
// - mapercentage, ngfpercentage, msepercentage, gdpercentage, ncpercentage, nmipercentage (use metric-percentages array or default 0.1)
// - alphaderivative, lambdaderivative, nuderivative, etc. (auto-derived from weights)
// - etavaluefixed, etavaluemoving (auto-estimated)
// - ngfevaluator (default 0)
// - nmibins (default 64)
// - ngfspacing (inferred from image)
// - fixedimagethreshold, metricoverlap, etc. (kept as advanced, commented defaults)


// ============================================================================
// STEP 3: Add auto-detection logic after parsing (around line 200)
// ============================================================================

po::store(po::parse_command_line(argc, argv, desc), vm);
po::notify(vm);

// === AUTO-DETECT FORMAT ===
std::string weightsStr = vm["weights"].as<std::string>();
std::string presetStr = vm["preset"].as<std::string>();
bool isNewFormat = CliParser::IsNewFormat(weightsStr, presetStr);

CliParser::MetricWeights weights;
CliParser::MetricDerivatives derivatives;
CliParser::SamplingPercentages sampling = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};

if (isNewFormat)
{
    // NEW FORMAT: parse preset or weights array
    if (!presetStr.empty())
    {
        weights = CliParser::GetPreset(presetStr);
        derivatives = CliParser::DeriveFromWeights(weights);
        std::cout << "Using preset: " << presetStr << std::endl;
    }
    else if (!weightsStr.empty())
    {
        weights = CliParser::ParseWeightsArray(weightsStr);
        derivatives = CliParser::DeriveFromWeights(weights);
        std::cout << "Using weights array: " << weightsStr << std::endl;
    }
    
    // Parse optional sampling percentages
    std::string samplingStr = vm["metric-percentages"].as<std::string>();
    if (!samplingStr.empty())
    {
        sampling = CliParser::ParseSamplingArray(samplingStr);
    }
}
else
{
    // OLD FORMAT: gather from individual parameters (backward compatible)
    weights.alpha = vm["alpha"].as<double>();
    weights.lambda = vm["lambda"].as<double>();
    weights.nu = vm["nu"].as<double>();
    weights.rho = vm["rho"].as<double>();
    weights.yota = vm["yota"].as<double>();
    weights.kappa = vm["labelkappa"].as<double>();
    weights.sigma = vm["sigma"].as<double>();
    
    derivatives.alphaDerivative = vm["alphaderivative"].as<double>();
    derivatives.lambdaDerivative = vm["lambdaderivative"].as<double>();
    derivatives.nuDerivative = vm["nuderivative"].as<double>();
    derivatives.rhoDerivative = vm["rhoderivative"].as<double>();
    derivatives.yotaDerivative = vm["yotaderivative"].as<double>();
    derivatives.kappaDerivative = vm["labelkappadervative"].as<double>();
    derivatives.sigmaDerivative = vm["sigmaderivative"].as<double>();
    
    sampling.ma = vm["mapercentage"].as<double>();
    sampling.ngf = vm["ngfpercentage"].as<double>();
    sampling.mse = vm["msepercentage"].as<double>();
    sampling.gd = vm["gdpercentage"].as<double>();
    sampling.nc = vm["ncpercentage"].as<double>();
    sampling.nmi = vm["nmipercentage"].as<double>();
    sampling.label = vm["labelsamples"].as<double>();
}


// ============================================================================
// STEP 4: Replace individual weight-setting calls with parsed values
// ============================================================================

// OLD (approx line 405):
// metric->SetAlpha(ALPHA);
// metric->SetAlphaDerivative(ALPHADERIVATIVE);
// metric->SetNu(NU);
// ... etc

// NEW:
metric->SetAlpha(weights.alpha);
metric->SetAlphaDerivative(derivatives.alphaDerivative);
metric->SetLambda(weights.lambda);
metric->SetLambdaDerivative(derivatives.lambdaDerivative);
metric->SetNu(weights.nu);
metric->SetNuDerivative(derivatives.nuDerivative);
metric->SetRho(weights.rho);
metric->SetRhoDerivative(derivatives.rhoDerivative);
metric->SetYota(weights.yota);
metric->SetYotaDerivative(derivatives.yotaDerivative);
metric->SetLabelKappa(weights.kappa);
metric->SetLabelKappaDerivative(derivatives.kappaDerivative);
metric->SetSigma(weights.sigma);
metric->SetSigmaDerivative(derivatives.sigmaDerivative);

// Sampling
const unsigned int numberOfPixels = fixedImage->GetLargestPossibleRegion().GetNumberOfPixels();
metric->SetMANumberOfSamples(static_cast<unsigned int>(numberOfPixels * sampling.ma));
metric->SetNGFNumberOfSamples(static_cast<unsigned int>(numberOfPixels * sampling.ngf));
metric->SetMSENumberOfSamples(static_cast<unsigned int>(numberOfPixels * sampling.mse));
metric->SetGDNumberOfSamples(static_cast<unsigned int>(numberOfPixels * sampling.gd));
metric->SetNCNumberOfSamples(static_cast<unsigned int>(numberOfPixels * sampling.nc));
metric->SetNMINumberOfSamples(static_cast<unsigned int>(numberOfPixels * sampling.nmi));


// ============================================================================
// STEP 5: Handle inferred parameters (NEW)
// ============================================================================

// Auto-estimate NGF noise (already exists as flag, just enable by default)
metric->SetAutoEstimateEta(true);   // Instead of requiring --etavaluefixed --etavaluemoving

// Infer NGF spacing from fixed image (NEW)
typename FixedImageType::SpacingType ngfSpacing = fixedImage->GetSpacing();
ngfSpacing *= 2.0;  // Scale to ~2x voxel spacing
metric->SetNGFSpacing(ngfSpacing);

// Use sensible defaults
metric->SetNGFPrecomputeGradient(false);
metric->SetDerivativeMode(vm["derivativemode"].as<int>());
metric->SetComputeOverlap(true);


// ============================================================================
// EXAMPLE USAGE
// ============================================================================

/*
// Minimalist (all defaults except images):
./3DRegAffine --fixed fixed.nii --moving moving.nii --output transform.txt

// With preset:
./3DRegAffine --fixed fixed.nii --moving moving.nii --output transform.txt \
  --preset multimodal

// Explicit weights (new format):
./3DRegAffine --fixed fixed.nii --moving moving.nii --output transform.txt \
  --weights "1.0,0.5,0,0,0,0,0"

// With custom sampling:
./3DRegAffine --fixed fixed.nii --moving moving.nii --output transform.txt \
  --preset singlemodal --metric-percentages "0.2,0.1,0.15,0.1,0.1,0.1,0.1"

// Advanced (old format, still works):
./3DRegAffine --fixed fixed.nii --moving moving.nii --output transform.txt \
  --alpha 1.0 --lambda 0.5 --nu 0 --iterations 500 --threads 4

// Tuning optimizer:
./3DRegAffine --fixed fixed.nii --moving moving.nii --output transform.txt \
  --preset rigid --iterations 2000 --step-length 0.05 --tolerance 1e-5
*/
