#define _USE_MATH_DEFINES
// 3DRegCalibrate.cxx
// Metric calibration tool for itkMplus multi-metric registration.
//
// Sweeps a Similarity3D (rigid + isotropic scale) transform over a
// configurable range of perturbations around the identity transform and
// records each sub-metric's weighted value and total gradient norm at every
// sample point.
//
// The resulting statistics (range, mean gradient) are used to suggest weight
// ratios that normalise the contribution of each active metric so that none
// dominates accidentally due to differing numerical scales.
//
// Output
//   --output-csv  : CSV file (axis, param_value, MI, NGF, MSE, NC, NMI, GD, Label, Total, GradNorm)
//   stdout        : human-readable summary table + suggested per-metric weights
//
// Eros Montin, 2024

#include "itkImageRegistrationMethod.h"
#include "../../../Metrics/Mplus/itkMplus.h"
#include "itkSimilarity3DTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImageFileReader.h"
#include "itkImageMaskSpatialObject.h"
#include "itkImage.h"
#include "itkVersor.h"

#include "../../MetricsConfig.h"
#include "../../LabelWeightsParser.h"
#include "../../../includes/RegistrationCommon.h"
#include "../../../includes/imageUtils.h"

#include <boost/program_options.hpp>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <iostream>
#include <functional>
#include <map>

namespace po = boost::program_options;

static const unsigned int ImageDimension = 3;
typedef float PixelType;
typedef itk::Image<PixelType, ImageDimension> ImageType;
typedef itk::Mplus<ImageType, ImageType> MetricType;
typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
typedef itk::Similarity3DTransform<double> TransformType;

// ─────────────────────────────────────────────────────────────────────────────
// One sample in a perturbation sweep
// ─────────────────────────────────────────────────────────────────────────────
struct SweepResult
{
    std::string axis;
    double paramValue;
    double total, mi, ngf, mse, nc, nmi, gd, label;
    double gradNorm;
};

// ─────────────────────────────────────────────────────────────────────────────
// Aggregate statistics over one set of SweepResults
// ─────────────────────────────────────────────────────────────────────────────
struct MetricStats
{
    double minVal, maxVal, range;
    double meanGradNorm;
};

static std::vector<double> Linspace(double a, double b, int n)
{
    std::vector<double> v(n);
    if (n == 1) { v[0] = (a + b) * 0.5; return v; }
    for (int i = 0; i < n; ++i)
        v[i] = a + i * (b - a) / (n - 1);
    return v;
}

// ─────────────────────────────────────────────────────────────────────────────
// Build Similarity3DTransform parameters from physical inputs.
//   rotDeg[3]: rotation in degrees around image X, Y, Z axes
//   trans[3]:  translation in mm
//   scale:     isotropic scale factor (1.0 = no scale)
//   center:    physical point used as centre of rotation
// ─────────────────────────────────────────────────────────────────────────────
static TransformType::ParametersType MakeParams(
    double rotXdeg, double rotYdeg, double rotZdeg,
    double tx, double ty, double tz,
    double scale,
    const TransformType::InputPointType& center)
{
    TransformType::Pointer T = TransformType::New();
    T->SetIdentity();
    T->SetCenter(center);

    const double deg2rad = M_PI / 180.0;

    itk::Vector<double, 3> axisX, axisY, axisZ;
    axisX.Fill(0); axisX[0] = 1;
    axisY.Fill(0); axisY[1] = 1;
    axisZ.Fill(0); axisZ[2] = 1;

    itk::Versor<double> vx, vy, vz;
    vx.Set(axisX, rotXdeg * deg2rad);
    vy.Set(axisY, rotYdeg * deg2rad);
    vz.Set(axisZ, rotZdeg * deg2rad);
    T->SetRotation(vz * vy * vx);

    TransformType::OutputVectorType t;
    t[0] = tx; t[1] = ty; t[2] = tz;
    T->SetTranslation(t);
    T->SetScale(scale);

    return T->GetParameters();
}

// ─────────────────────────────────────────────────────────────────────────────
// Evaluate metric at a given parameter set. Returns false on exception.
// ─────────────────────────────────────────────────────────────────────────────
static bool EvaluateMetric(
    MetricType* metric,
    const TransformType::ParametersType& params,
    bool computeDerivative,
    double& total,
    double& mi, double& ngf, double& mse,
    double& nc, double& nmi, double& gd, double& label,
    double& gradNorm)
{
    try
    {
        if (computeDerivative)
        {
            MetricType::MeasureType val;
            MetricType::DerivativeType deriv(params.Size());
            metric->GetValueAndDerivative(params, val, deriv);
            total = static_cast<double>(val);
            // Gradient norm (L2)
            double s = 0.0;
            for (unsigned i = 0; i < deriv.Size(); ++i)
                s += deriv[i] * deriv[i];
            gradNorm = std::sqrt(s);
        }
        else
        {
            total = static_cast<double>(metric->GetValue(params));
            gradNorm = 0.0;
        }

        mi    = metric->GetLastValMI();
        ngf   = metric->GetLastValNGF();
        mse   = metric->GetLastValMSE();
        nc    = metric->GetLastValNC();
        nmi   = metric->GetLastValNMI();
        gd    = metric->GetLastValGD();
        label = metric->GetLastValLabel();
        return true;
    }
    catch (const std::exception& ex)
    {
        std::cerr << "[Calibrate] Metric evaluation failed: " << ex.what() << "\n";
        return false;
    }
    catch (...)
    {
        std::cerr << "[Calibrate] Metric evaluation failed (unknown exception)\n";
        return false;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Run a 1D perturbation sweep along one axis.
//   paramsFn: maps a scalar perturbation value to a full parameter vector
// ─────────────────────────────────────────────────────────────────────────────
static std::vector<SweepResult> RunSweep(
    MetricType* metric,
    const std::string& axisName,
    const std::vector<double>& values,
    bool computeDerivative,
    std::function<TransformType::ParametersType(double)> paramsFn)
{
    std::vector<SweepResult> results;
    results.reserve(values.size());

    std::cout << "  Sweeping " << axisName << " (" << values.size() << " pts)..." << std::flush;

    for (double v : values)
    {
        SweepResult r;
        r.axis       = axisName;
        r.paramValue = v;

        bool ok = EvaluateMetric(metric, paramsFn(v), computeDerivative,
                                 r.total, r.mi, r.ngf, r.mse,
                                 r.nc, r.nmi, r.gd, r.label, r.gradNorm);
        if (ok)
            results.push_back(r);
    }

    std::cout << " done (" << results.size() << " valid)\n";
    return results;
}

// ─────────────────────────────────────────────────────────────────────────────
// Compute statistics for one metric field across all sweep results
// ─────────────────────────────────────────────────────────────────────────────
static MetricStats ComputeStats(
    const std::vector<SweepResult>& all,
    std::function<double(const SweepResult&)> valueGetter)
{
    MetricStats s;
    s.minVal = std::numeric_limits<double>::max();
    s.maxVal = -std::numeric_limits<double>::max();
    double gradSum = 0.0;

    for (const auto& r : all)
    {
        double v = valueGetter(r);
        if (std::isfinite(v))
        {
            s.minVal = std::min(s.minVal, v);
            s.maxVal = std::max(s.maxVal, v);
        }
        if (std::isfinite(r.gradNorm))
            gradSum += r.gradNorm;
    }

    s.range        = (all.empty()) ? 0.0 : s.maxVal - s.minVal;
    s.meanGradNorm = all.empty() ? 0.0 : gradSum / static_cast<double>(all.size());

    if (!std::isfinite(s.minVal)) s.minVal = 0.0;
    if (!std::isfinite(s.maxVal)) s.maxVal = 0.0;
    return s;
}

// ─────────────────────────────────────────────────────────────────────────────
// Write CSV
// ─────────────────────────────────────────────────────────────────────────────
static void WriteCSV(const std::string& path, const std::vector<SweepResult>& all)
{
    std::ofstream f(path);
    if (!f.is_open())
    {
        std::cerr << "[Calibrate] Cannot open CSV output file: " << path << "\n";
        return;
    }
    f << "axis,param_value,MI,NGF,MSE,NC,NMI,GD,Label,Total,GradNorm\n";
    f << std::fixed << std::setprecision(8);
    for (const auto& r : all)
        f << r.axis << "," << r.paramValue << ","
          << r.mi << "," << r.ngf << "," << r.mse << "," << r.nc << ","
          << r.nmi << "," << r.gd << "," << r.label << "," << r.total << ","
          << r.gradNorm << "\n";
    std::cout << "[Calibrate] CSV written to: " << path << "\n";
}

// ─────────────────────────────────────────────────────────────────────────────
// Print summary table and suggested weights
// ─────────────────────────────────────────────────────────────────────────────
static void PrintSummary(
    const std::vector<SweepResult>& all,
    bool computeDerivative)
{
    struct MetricEntry {
        std::string name;
        std::function<double(const SweepResult&)> getter;
    };

    const std::vector<MetricEntry> metrics = {
        {"MI",    [](const SweepResult& r){ return r.mi;    }},
        {"NGF",   [](const SweepResult& r){ return r.ngf;   }},
        {"MSE",   [](const SweepResult& r){ return r.mse;   }},
        {"NC",    [](const SweepResult& r){ return r.nc;    }},
        {"NMI",   [](const SweepResult& r){ return r.nmi;   }},
        {"GD",    [](const SweepResult& r){ return r.gd;    }},
        {"Label", [](const SweepResult& r){ return r.label; }},
        {"Total", [](const SweepResult& r){ return r.total; }},
    };

    std::cout << "\n";
    std::cout << std::string(80, '=') << "\n";
    std::cout << "  CALIBRATION SUMMARY\n";
    std::cout << std::string(80, '=') << "\n";
    std::cout << std::left
              << std::setw(8)  << "Metric"
              << std::setw(14) << "Min"
              << std::setw(14) << "Max"
              << std::setw(14) << "Range";
    if (computeDerivative)
        std::cout << std::setw(14) << "MeanGradNorm";
    std::cout << std::setw(14) << "SuggestedW"
              << "\n";
    std::cout << std::string(80, '-') << "\n";

    std::cout << std::fixed << std::setprecision(5);

    // Suggested weight = 1 / range (value range normalisation).
    // For derivative mode: 1 / meanGradNorm.
    // The weight for MI is used as the reference (=1.0); all others are scaled to MI.
    double miRange = 1.0;
    for (const auto& e : metrics)
    {
        if (e.name == "MI")
        {
            auto s = ComputeStats(all, e.getter);
            miRange = (s.range > 1e-10) ? s.range : 1.0;
            break;
        }
    }

    for (const auto& e : metrics)
    {
        auto s = ComputeStats(all, e.getter);
        const double refRange = (s.range > 1e-10) ? s.range : 1.0;
        const double suggestedW = miRange / refRange;  // relative to MI

        std::cout << std::setw(8)  << e.name
                  << std::setw(14) << s.minVal
                  << std::setw(14) << s.maxVal
                  << std::setw(14) << s.range;
        if (computeDerivative)
            std::cout << std::setw(14) << s.meanGradNorm;
        if (e.name == "Total")
            std::cout << std::setw(14) << "(combined)";
        else
            std::cout << std::setw(14) << suggestedW;
        std::cout << "\n";
    }

    std::cout << std::string(80, '=') << "\n";
    std::cout << "  Note: SuggestedW is relative to MI (MI=1.0), computed as\n"
              << "        range(MI)/range(metric).  Use these as a starting point\n"
              << "        for --alpha/--lambda/--nu/--rho/--yota/--sigma.\n";
    std::cout << std::string(80, '=') << "\n";
}

// ─────────────────────────────────────────────────────────────────────────────
// main
// ─────────────────────────────────────────────────────────────────────────────
int main(int argc, char* argv[])
{
    po::options_description desc(
        "3DRegCalibrate — Multi-metric calibration tool\n"
        "Sweeps a rigid/similarity perturbation and measures each sub-metric's\n"
        "value range and gradient norm to suggest balanced weight ratios.\n\n"
        "Allowed options");

    desc.add_options()
        ("help,h", "produce help message")
        ("fixedimage,f",  po::value<std::string>(), "Fixed image filename (required)")
        ("movingimage,m", po::value<std::string>(), "Moving image filename (required)")
        ("focusroi",      po::value<std::string>()->default_value("N"),
             "Focus ROI mask (N = none)")
        ("fixedlabelmap",  po::value<std::string>()->default_value("N"),
             "Fixed label map (N = none)")
        ("movinglabelmap", po::value<std::string>()->default_value("N"),
             "Moving label map (N = none)")

        // Metric weights (same names as registration tools)
        ("alpha",          po::value<double>()->default_value(1.0),  "MI weight")
        ("alphaderivative",po::value<double>()->default_value(1.0),  "MI derivative weight")
        ("mattesnumberofbins,b", po::value<int>()->default_value(64),"Mattes bins")
        ("mattespercentage,p", po::value<double>()->default_value(0.05), "Mattes sampling fraction")
        ("lambda",         po::value<double>()->default_value(0.0),  "NGF weight")
        ("lambdaderivative",po::value<double>()->default_value(0.0), "NGF derivative weight")
        ("ngfpercentage",  po::value<double>()->default_value(0.05), "NGF sampling fraction")
        ("NGFevaluator",   po::value<int>()->default_value(0),       "NGF evaluator (0=scalar)")
        ("etavaluefixed",  po::value<double>()->default_value(-1),   "NGF eta fixed (-1=auto)")
        ("etavaluemoving", po::value<double>()->default_value(-1),   "NGF eta moving (-1=auto)")
        ("nu",             po::value<double>()->default_value(0.0),  "MSE weight")
        ("nuderivative",   po::value<double>()->default_value(0.0),  "MSE derivative weight")
        ("msepercentage",  po::value<double>()->default_value(0.05), "MSE sampling fraction")
        ("rho",            po::value<double>()->default_value(0.0),  "GD weight")
        ("rhoderivative",  po::value<double>()->default_value(0.0),  "GD derivative weight")
        ("gdpercentage",   po::value<double>()->default_value(0.05), "GD sampling fraction")
        ("yota",           po::value<double>()->default_value(0.0),  "NC weight")
        ("yotaderivative", po::value<double>()->default_value(0.0),  "NC derivative weight")
        ("ncpercentage",   po::value<double>()->default_value(0.05), "NC sampling fraction")
        ("sigma",          po::value<double>()->default_value(0.0),  "NMI weight")
        ("sigmaderivative",po::value<double>()->default_value(0.0),  "NMI derivative weight")
        ("nmibins",        po::value<int>()->default_value(64),      "NMI bins")
        ("nmipercentage",  po::value<double>()->default_value(0.05), "NMI sampling fraction")
        ("normalizemse",   po::value<bool>()->default_value(false),  "Normalise MSE")
        ("normalizegd",    po::value<bool>()->default_value(false),  "Normalise GD")
        ("labelkappa",     po::value<double>()->default_value(0.0),  "Label kappa weight")
        ("labelkappaderiv",po::value<double>()->default_value(0.0),  "Label kappa deriv weight")
        ("labelkappavec",  po::value<std::string>()->default_value(""), "Per-label kappa weights")
        ("labelkappaderivvec", po::value<std::string>()->default_value(""), "Per-label kappa deriv")
        ("labelsamples",   po::value<double>()->default_value(0.05), "Label sampling fraction")
        ("labeldistmax",   po::value<double>()->default_value(20.0), "Label distance clamp (mm)")

        // Sweep configuration
        ("rot-max",    po::value<double>()->default_value(20.0),
             "Max rotation per axis in degrees (sweep from -rot-max to +rot-max)")
        ("rot-steps",  po::value<int>()->default_value(9),
             "Number of rotation sweep points per axis (odd → includes 0)")
        ("trans-max",  po::value<double>()->default_value(10.0),
             "Max translation per axis in mm (sweep from -trans-max to +trans-max)")
        ("trans-steps",po::value<int>()->default_value(9),
             "Number of translation sweep points per axis")
        ("scale-max",  po::value<double>()->default_value(0.1),
             "Max scale deviation (sweep from 1-scale-max to 1+scale-max)")
        ("scale-steps",po::value<int>()->default_value(5),
             "Number of scale sweep points (0 = skip scale sweep)")
        ("compute-derivative", po::value<bool>()->default_value(false),
             "Also compute gradient norms (slower; requires derivative weights > 0)")

        // Misc
        ("numberofthreads", po::value<int>()->default_value(2), "Number of threads")
        ("bsplinecaching",  po::value<bool>()->default_value(false),
             "B-spline weight caching (false recommended for calibration)")
        ("workingresolution", po::value<std::string>()->default_value("0,0,0"),
             "Resample inputs to this spacing in mm (0,0,0 = native)")
        ("output-csv", po::value<std::string>()->default_value(""),
             "Path for CSV output file (empty = no CSV)")
        ("metric-overlap", po::value<bool>()->default_value(true),
             "Restrict metric to overlapping image region")
        ("ngfspacing", po::value<std::string>()->default_value("4,4,4"),
             "NGF spacing (x,y,z mm)")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help") || !vm.count("fixedimage") || !vm.count("movingimage"))
    {
        std::cout << desc << "\n";
        return 1;
    }

    // ── Load images ───────────────────────────────────────────────────────────
    typedef itk::ImageFileReader<ImageType> ReaderType;

    ReaderType::Pointer fixedReader  = ReaderType::New();
    ReaderType::Pointer movingReader = ReaderType::New();
    fixedReader->SetFileName(vm["fixedimage"].as<std::string>());
    movingReader->SetFileName(vm["movingimage"].as<std::string>());
    fixedReader->Update();
    movingReader->Update();

    ImageType::ConstPointer fixedImage  = fixedReader->GetOutput();
    ImageType::ConstPointer movingImage = movingReader->GetOutput();

    // Optional working-resolution resampling
    ImageType::SpacingType workingSpacing;
    bool useWorkingRes = false;
    if (!RegCommon::ParseOptionalSpacing<ImageType::SpacingType, ImageDimension>(
            vm["workingresolution"].as<std::string>(), "workingresolution",
            workingSpacing, useWorkingRes))
        return EXIT_FAILURE;

    if (useWorkingRes)
    {
        std::cout << "[Calibrate] Resampling to working resolution " << workingSpacing << "\n";
        if (!RegCommon::SpacingEquals(fixedImage->GetSpacing(), workingSpacing))
            fixedImage = RegCommon::ResampleScalarImageToSpacing<ImageType>(fixedImage, workingSpacing);
        if (!RegCommon::SpacingEquals(movingImage->GetSpacing(), workingSpacing))
            movingImage = RegCommon::ResampleScalarImageToSpacing<ImageType>(
                movingImage, workingSpacing, 0.0);
    }

    // ── Compute image centre for rotation ────────────────────────────────────
    TransformType::InputPointType imageCenter;
    {
        const auto& sp = fixedImage->GetSpacing();
        const auto& sz = fixedImage->GetLargestPossibleRegion().GetSize();
        const auto& orig = fixedImage->GetOrigin();
        const auto& dir  = fixedImage->GetDirection();
        for (unsigned j = 0; j < ImageDimension; ++j)
        {
            imageCenter[j] = orig[j];
            for (unsigned i = 0; i < ImageDimension; ++i)
                imageCenter[j] += dir[j][i] * sp[i] * (sz[i] - 1) * 0.5;
        }
    }
    std::cout << "[Calibrate] Rotation centre: ("
              << imageCenter[0] << ", " << imageCenter[1] << ", " << imageCenter[2] << ")\n";

    // ── Set up Mplus metric ──────────────────────────────────────────────────
    MetricType::Pointer     metric       = MetricType::New();
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    TransformType::Pointer  transform    = TransformType::New();
    transform->SetIdentity();
    transform->SetCenter(imageCenter);

    metric->SetFixedImage(fixedImage);
    metric->SetMovingImage(const_cast<ImageType*>(movingImage.GetPointer()));
    metric->SetTransform(transform);
    metric->SetInterpolator(interpolator);
    metric->SetFixedImageRegion(fixedImage->GetBufferedRegion());
    metric->SetNumberOfThreads(vm["numberofthreads"].as<int>());
    metric->SetUseCachingOfBSplineWeights(vm["bsplinecaching"].as<bool>());
    metric->SetComputeOverlap(vm["metric-overlap"].as<bool>());

    // NGF spacing
    {
        auto s = vm["ngfspacing"].as<std::string>();
        std::replace(s.begin(), s.end(), ',', ' ');
        std::istringstream iss(s);
        std::vector<double> tmp(std::istream_iterator<double>{iss}, {});
        if (tmp.size() == ImageDimension)
        {
            ImageType::SpacingType ngfSp;
            for (unsigned i = 0; i < ImageDimension; ++i) ngfSp[i] = tmp[i];
            metric->SetNGFSpacing(ngfSp);
        }
    }

    // Sub-metric weights
    metric->SetAlpha(vm["alpha"].as<double>());
    metric->SetAlphaDerivative(vm["alphaderivative"].as<double>());
    metric->SetLambda(vm["lambda"].as<double>());
    metric->SetLambdaDerivative(vm["lambdaderivative"].as<double>());
    metric->SetNu(vm["nu"].as<double>());
    metric->SetNuDerivative(vm["nuderivative"].as<double>());
    metric->SetRho(vm["rho"].as<double>());
    metric->SetRhoDerivative(vm["rhoderivative"].as<double>());
    metric->SetYota(vm["yota"].as<double>());
    metric->SetYotaDerivative(vm["yotaderivative"].as<double>());
    metric->SetSigma(vm["sigma"].as<double>());
    metric->SetSigmaDerivative(vm["sigmaderivative"].as<double>());
    metric->SetNormalizeMSE(vm["normalizemse"].as<bool>());
    metric->SetNormalizeGD(vm["normalizegd"].as<bool>());

    if (vm["etavaluefixed"].as<double>() == -1 ||
        vm["etavaluemoving"].as<double>() == -1)
        metric->SetAutoEstimateEta(true);
    else
    {
        metric->SetFixedEta(vm["etavaluefixed"].as<double>());
        metric->SetMovingEta(vm["etavaluemoving"].as<double>());
    }
    metric->SetEvaluator(vm["NGFevaluator"].as<int>());

    // Sampling counts
    const unsigned int nPix = fixedImage->GetLargestPossibleRegion().GetNumberOfPixels();
    metric->SetMANumberOfSamples(static_cast<unsigned int>(nPix * vm["mattespercentage"].as<double>()));
    metric->SetBinNumbers(vm["mattesnumberofbins"].as<int>());
    metric->SetNGFNumberOfSamples(static_cast<unsigned int>(nPix * vm["ngfpercentage"].as<double>()));
    metric->SetMSENumberOfSamples(static_cast<unsigned int>(nPix * vm["msepercentage"].as<double>()));
    metric->SetGDNumberOfSamples(static_cast<unsigned int>(nPix * vm["gdpercentage"].as<double>()));
    metric->SetNCNumberOfSamples(static_cast<unsigned int>(nPix * vm["ncpercentage"].as<double>()));
    metric->SetNMIBinNumbers(vm["nmibins"].as<int>());
    metric->SetNMINumberOfSamples(static_cast<unsigned int>(nPix * vm["nmipercentage"].as<double>()));
    metric->SetLabelNumberOfSamples(RegCommon::ResolveLabelSampleCount(
        vm["labelsamples"].as<double>(), nPix));
    metric->SetLabelDistanceMax(vm["labeldistmax"].as<double>());

    // Optional focus ROI mask
    const std::string focusROI = vm["focusroi"].as<std::string>();
    if (focusROI != "N" && !focusROI.empty())
    {
        typedef itk::ImageFileReader<itk::Image<unsigned char, ImageDimension>> MaskReaderType;
        MaskReaderType::Pointer mr = MaskReaderType::New();
        mr->SetFileName(focusROI);
        mr->Update();
        typedef itk::ImageMaskSpatialObject<ImageDimension> MaskSOType;
        MaskSOType::Pointer maskSO = MaskSOType::New();
        maskSO->SetImage(mr->GetOutput());
        maskSO->Update();
        metric->SetFixedImageMask(maskSO);
        std::cout << "[Calibrate] Focus ROI mask: " << focusROI << "\n";
    }

    // Optional label maps
    typedef itk::Image<short, ImageDimension> LabelImageType;
    const std::string fixedLabelFN  = vm["fixedlabelmap"].as<std::string>();
    const std::string movingLabelFN = vm["movinglabelmap"].as<std::string>();
    if (fixedLabelFN != "N" && movingLabelFN != "N")
    {
        typedef itk::ImageFileReader<LabelImageType> LabelReaderType;
        auto flr = LabelReaderType::New(); flr->SetFileName(fixedLabelFN);  flr->Update();
        auto mlr = LabelReaderType::New(); mlr->SetFileName(movingLabelFN); mlr->Update();
        LabelImageType::Pointer fl = flr->GetOutput(); fl->DisconnectPipeline();
        LabelImageType::Pointer ml = mlr->GetOutput(); ml->DisconnectPipeline();
        metric->SetFixedLabelMap(fl.GetPointer());
        metric->SetMovingLabelMap(ml.GetPointer());
        metric->SetLabelKappa(vm["labelkappa"].as<double>());
        metric->SetLabelKappaDerivative(vm["labelkappaderiv"].as<double>());
        const auto kv  = RegCommon::ParseLabelWeights(vm["labelkappavec"].as<std::string>());
        const auto kdv = RegCommon::ParseLabelWeights(vm["labelkappaderivvec"].as<std::string>());
        if (!kv.empty())  metric->SetLabelKappaWeights(kv);
        if (!kdv.empty()) metric->SetLabelKappaDerivativeWeights(kdv);
        std::cout << "[Calibrate] Label maps loaded\n";
    }

    // ── Initialize metric ─────────────────────────────────────────────────────
    std::cout << "[Calibrate] Initialising metric...\n";
    try { metric->Initialize(); }
    catch (const itk::ExceptionObject& e)
    {
        std::cerr << "Metric initialisation failed:\n" << e.GetDescription() << "\n";
        return EXIT_FAILURE;
    }
    std::cout << "[Calibrate] Metric initialised. " << nPix << " voxels.\n\n";

    // ── Sweep configuration ───────────────────────────────────────────────────
    const double rotMax    = vm["rot-max"].as<double>();
    const int    rotSteps  = vm["rot-steps"].as<int>();
    const double transMax  = vm["trans-max"].as<double>();
    const int    transSteps = vm["trans-steps"].as<int>();
    const double scaleMax  = vm["scale-max"].as<double>();
    const int    scaleSteps = vm["scale-steps"].as<int>();
    const bool   doDerivative = vm["compute-derivative"].as<bool>();

    const auto rotVals   = Linspace(-rotMax,   +rotMax,   rotSteps);
    const auto transVals = Linspace(-transMax, +transMax, transSteps);
    const auto scaleVals = (scaleSteps > 0)
        ? Linspace(1.0 - scaleMax, 1.0 + scaleMax, scaleSteps)
        : std::vector<double>{};

    // ── Run sweeps ────────────────────────────────────────────────────────────
    std::cout << "[Calibrate] Running perturbation sweeps...\n";

    std::vector<SweepResult> allResults;

    auto append = [&allResults](std::vector<SweepResult>&& v)
    {
        allResults.insert(allResults.end(),
                          std::make_move_iterator(v.begin()),
                          std::make_move_iterator(v.end()));
    };

    // Rotation sweeps
    append(RunSweep(metric, "rotX", rotVals, doDerivative,
        [&](double deg){ return MakeParams(deg, 0, 0, 0, 0, 0, 1.0, imageCenter); }));
    append(RunSweep(metric, "rotY", rotVals, doDerivative,
        [&](double deg){ return MakeParams(0, deg, 0, 0, 0, 0, 1.0, imageCenter); }));
    append(RunSweep(metric, "rotZ", rotVals, doDerivative,
        [&](double deg){ return MakeParams(0, 0, deg, 0, 0, 0, 1.0, imageCenter); }));

    // Translation sweeps
    append(RunSweep(metric, "transX", transVals, doDerivative,
        [&](double mm){ return MakeParams(0, 0, 0, mm, 0,  0,  1.0, imageCenter); }));
    append(RunSweep(metric, "transY", transVals, doDerivative,
        [&](double mm){ return MakeParams(0, 0, 0, 0,  mm, 0,  1.0, imageCenter); }));
    append(RunSweep(metric, "transZ", transVals, doDerivative,
        [&](double mm){ return MakeParams(0, 0, 0, 0,  0,  mm, 1.0, imageCenter); }));

    // Scale sweep (optional)
    if (!scaleVals.empty())
        append(RunSweep(metric, "scale", scaleVals, doDerivative,
            [&](double s){ return MakeParams(0, 0, 0, 0, 0, 0, s, imageCenter); }));

    std::cout << "[Calibrate] Total sweep results: " << allResults.size() << "\n";

    // ── CSV output ───────────────────────────────────────────────────────────
    const std::string csvPath = vm["output-csv"].as<std::string>();
    if (!csvPath.empty())
        WriteCSV(csvPath, allResults);

    // ── Print summary ─────────────────────────────────────────────────────────
    PrintSummary(allResults, doDerivative);

    return EXIT_SUCCESS;
}
