#define _USE_MATH_DEFINES
// 3DRegCalibrate.cxx
// Metric calibration tool for itkMplus multi-metric registration.
//
// Sweeps a Similarity3D (rigid + isotropic scale) transform over a
// configurable range of perturbations around the identity transform and
// records each sub-metric's unit-weight value and derivative scale at every
// sample point.
//
// The resulting value and derivative statistics are used to suggest value and
// derivative weight ratios so that no metric dominates accidentally due to
// differing numerical scales.
//
// Output
//   --output-csv  : CSV file with values, absolute values, ranges, and derivative norms
//   stdout        : human-readable summary table + relative value/derivative weights
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
#include <array>
#include <string>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <iostream>
#include <functional>
#include <map>
#include <limits>
#include <stdexcept>

namespace po = boost::program_options;

static const unsigned int ImageDimension = 3;
typedef float PixelType;
typedef itk::Image<PixelType, ImageDimension> ImageType;
typedef itk::Mplus<ImageType, ImageType> MetricType;
typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
typedef itk::Similarity3DTransform<double> TransformType;

static const unsigned int MetricCount = 7;
enum MetricId
{
    MetricMI = 0,
    MetricNGF = 1,
    MetricMSE = 2,
    MetricGD = 3,
    MetricNC = 4,
    MetricNMI = 5,
    MetricLabel = 6
};

static const std::array<std::string, MetricCount> MetricNames = {{
    "MI", "NGF", "MSE", "GD", "NC", "NMI", "Label"
}};

static bool IsNonZero(double v)
{
    return std::isfinite(v) && std::fabs(v) > 1.0e-12;
}

struct CalibrationSelection
{
    std::array<bool, MetricCount> value;
    std::array<bool, MetricCount> derivative;

    CalibrationSelection()
    {
        value.fill(false);
        derivative.fill(false);
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// One sample in a perturbation sweep
// ─────────────────────────────────────────────────────────────────────────────
struct SweepResult
{
    std::string axis;
    double paramValue;
    double total;
    double absTotal;
    double gradNorm;

    std::array<double, MetricCount> value;
    std::array<double, MetricCount> absValue;
    std::array<double, MetricCount> derivNorm;
    std::array<double, MetricCount> derivMeanAbs;
    std::array<double, MetricCount> derivAbsRange;

    SweepResult() : paramValue(0.0), total(0.0), absTotal(0.0), gradNorm(0.0)
    {
        value.fill(0.0);
        absValue.fill(0.0);
        derivNorm.fill(0.0);
        derivMeanAbs.fill(0.0);
        derivAbsRange.fill(0.0);
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// Aggregate statistics over one set of SweepResults
// ─────────────────────────────────────────────────────────────────────────────
struct MetricStats
{
    double minVal, maxVal, range;
    double minAbs, maxAbs, absRange, meanAbs;
    double meanDerivNorm, maxDerivNorm;
    double meanDerivMeanAbs, meanDerivAbsRange;
    unsigned int count;

    MetricStats()
        : minVal(std::numeric_limits<double>::max()),
          maxVal(-std::numeric_limits<double>::max()),
          range(0.0),
          minAbs(std::numeric_limits<double>::max()),
          maxAbs(0.0),
          absRange(0.0),
          meanAbs(0.0),
          meanDerivNorm(0.0),
          maxDerivNorm(0.0),
          meanDerivMeanAbs(0.0),
          meanDerivAbsRange(0.0),
          count(0)
    {}
};

struct SummaryEntry
{
    unsigned int index;
    std::string name;
    bool valueEnabled;
    bool derivativeEnabled;
    MetricStats stats;
    double valueScale;
    double derivativeScale;
    double relValueWeight;
    double relDerivativeWeight;

    SummaryEntry()
        : index(0), valueEnabled(false), derivativeEnabled(false),
          valueScale(0.0), derivativeScale(0.0),
          relValueWeight(0.0), relDerivativeWeight(0.0)
    {}
};

static std::vector<double> Linspace(double a, double b, int n)
{
    std::vector<double> v(n);
    if (n == 1) { v[0] = (a + b) * 0.5; return v; }
    for (int i = 0; i < n; ++i)
        v[i] = a + i * (b - a) / (n - 1);
    return v;
}

static std::string Trim(const std::string& s)
{
    const auto begin = s.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos) return "";
    const auto end = s.find_last_not_of(" \t\r\n");
    return s.substr(begin, end - begin + 1);
}

static std::vector<double> ParseDoubleList(const std::string& raw,
                                           unsigned int minCount,
                                           unsigned int maxCount,
                                           const std::string& name)
{
    std::vector<double> out;
    std::istringstream ss(raw);
    std::string token;
    while (std::getline(ss, token, ','))
    {
        token = Trim(token);
        if (token.empty())
            throw std::runtime_error("Empty value in " + name);
        out.push_back(std::stod(token));
    }

    if (out.size() < minCount || out.size() > maxCount)
    {
        std::ostringstream msg;
        msg << name << " expects ";
        if (minCount == maxCount)
            msg << minCount;
        else
            msg << minCount << "-" << maxCount;
        msg << " values, got " << out.size();
        throw std::runtime_error(msg.str());
    }
    return out;
}

static void ApplySelectorArray(const std::vector<double>& values,
                               std::array<bool, MetricCount>& target)
{
    target.fill(false);
    for (unsigned int i = 0; i < values.size() && i < MetricCount; ++i)
        target[i] = IsNonZero(values[i]);
}

static void PrintSelection(const CalibrationSelection& selection)
{
    std::cout << "[Calibrate] Metric value selection:";
    bool any = false;
    for (unsigned int i = 0; i < MetricCount; ++i)
    {
        if (!selection.value[i]) continue;
        std::cout << (any ? "," : " ") << MetricNames[i];
        any = true;
    }
    if (!any) std::cout << " none";
    std::cout << "\n";

    std::cout << "[Calibrate] Metric derivative selection:";
    any = false;
    for (unsigned int i = 0; i < MetricCount; ++i)
    {
        if (!selection.derivative[i]) continue;
        std::cout << (any ? "," : " ") << MetricNames[i];
        any = true;
    }
    if (!any) std::cout << " none";
    std::cout << "\n";
}

static CalibrationSelection BuildSelection(const po::variables_map& vm)
{
    CalibrationSelection selection;

    selection.value[MetricMI]    = IsNonZero(vm["alpha"].as<double>());
    selection.value[MetricNGF]   = IsNonZero(vm["lambda"].as<double>());
    selection.value[MetricMSE]   = IsNonZero(vm["nu"].as<double>());
    selection.value[MetricGD]    = IsNonZero(vm["rho"].as<double>());
    selection.value[MetricNC]    = IsNonZero(vm["yota"].as<double>());
    selection.value[MetricNMI]   = IsNonZero(vm["sigma"].as<double>());
    selection.value[MetricLabel] = IsNonZero(vm["labelkappa"].as<double>()) ||
                                   !vm["labelkappavec"].as<std::string>().empty();

    selection.derivative[MetricMI]    = IsNonZero(vm["alphaderivative"].as<double>());
    selection.derivative[MetricNGF]   = IsNonZero(vm["lambdaderivative"].as<double>());
    selection.derivative[MetricMSE]   = IsNonZero(vm["nuderivative"].as<double>());
    selection.derivative[MetricGD]    = IsNonZero(vm["rhoderivative"].as<double>());
    selection.derivative[MetricNC]    = IsNonZero(vm["yotaderivative"].as<double>());
    selection.derivative[MetricNMI]   = IsNonZero(vm["sigmaderivative"].as<double>());
    selection.derivative[MetricLabel] = IsNonZero(vm["labelkappaderiv"].as<double>()) ||
                                        !vm["labelkappaderivvec"].as<std::string>().empty();

    const std::string metricsArray = vm["metrics"].as<std::string>();
    if (!metricsArray.empty())
    {
        const auto values = ParseDoubleList(metricsArray, 6, 7, "--metrics");
        ApplySelectorArray(values, selection.value);
        selection.derivative = selection.value;
        std::cout << "[Calibrate] Parsed --metrics as calibration selector mask\n";
    }

    const std::string derivativeArray = vm["metric-derivatives"].as<std::string>();
    if (!derivativeArray.empty())
    {
        const auto values = ParseDoubleList(derivativeArray, 6, 7, "--metric-derivatives");
        ApplySelectorArray(values, selection.derivative);
        std::cout << "[Calibrate] Parsed --metric-derivatives as derivative selector mask\n";
    }

    return selection;
}

static std::array<double, MetricCount> BuildSampling(const po::variables_map& vm)
{
    std::array<double, MetricCount> sampling = {{
        vm["mattespercentage"].as<double>(),
        vm["ngfpercentage"].as<double>(),
        vm["msepercentage"].as<double>(),
        vm["gdpercentage"].as<double>(),
        vm["ncpercentage"].as<double>(),
        vm["nmipercentage"].as<double>(),
        vm["labelsamples"].as<double>()
    }};

    const std::string samplingArray = vm["metric-sampling"].as<std::string>();
    if (!samplingArray.empty())
    {
        const auto values = ParseDoubleList(samplingArray, 6, 7, "--metric-sampling");
        for (unsigned int i = 0; i < values.size() && i < MetricCount; ++i)
            sampling[i] = values[i];
        std::cout << "[Calibrate] Parsed --metric-sampling array\n";
    }

    return sampling;
}

static double UnitIf(bool enabled)
{
    return enabled ? 1.0 : 0.0;
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
    SweepResult& result)
{
    try
    {
        if (computeDerivative)
        {
            MetricType::MeasureType val;
            MetricType::DerivativeType deriv(params.Size());
            metric->GetValueAndDerivative(params, val, deriv);
            result.total = static_cast<double>(val);
            // Gradient norm (L2)
            double s = 0.0;
            for (unsigned i = 0; i < deriv.Size(); ++i)
                s += deriv[i] * deriv[i];
            result.gradNorm = std::sqrt(s);

            const auto& stats = metric->GetLastDerivativeStats();
            for (unsigned int i = 0; i < MetricCount; ++i)
            {
                const auto it = stats.find(MetricNames[i]);
                if (it == stats.end()) continue;
                result.derivNorm[i] = it->second.norm;
                result.derivMeanAbs[i] = it->second.meanAbs;
                result.derivAbsRange[i] = it->second.absRange;
            }
        }
        else
        {
            result.total = static_cast<double>(metric->GetValue(params));
            result.gradNorm = 0.0;
        }

        result.value[MetricMI]    = metric->GetLastValMI();
        result.value[MetricNGF]   = metric->GetLastValNGF();
        result.value[MetricMSE]   = metric->GetLastValMSE();
        result.value[MetricGD]    = metric->GetLastValGD();
        result.value[MetricNC]    = metric->GetLastValNC();
        result.value[MetricNMI]   = metric->GetLastValNMI();
        result.value[MetricLabel] = metric->GetLastValLabel();

        for (unsigned int i = 0; i < MetricCount; ++i)
            result.absValue[i] = std::fabs(result.value[i]);
        result.absTotal = std::fabs(result.total);
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

        bool ok = EvaluateMetric(metric, paramsFn(v), computeDerivative, r);
        if (ok)
            results.push_back(r);
    }

    std::cout << " done (" << results.size() << " valid)\n";
    return results;
}

// ─────────────────────────────────────────────────────────────────────────────
// Compute statistics for one metric field across all sweep results
// ─────────────────────────────────────────────────────────────────────────────
static MetricStats ComputeMetricStats(
    const std::vector<SweepResult>& all,
    unsigned int metricIndex)
{
    MetricStats s;
    double absSum = 0.0;
    double derivNormSum = 0.0;
    double derivMeanAbsSum = 0.0;
    double derivAbsRangeSum = 0.0;

    for (const auto& r : all)
    {
        const double v = r.value[metricIndex];
        const double av = r.absValue[metricIndex];
        if (std::isfinite(v))
        {
            s.minVal = std::min(s.minVal, v);
            s.maxVal = std::max(s.maxVal, v);
            s.minAbs = std::min(s.minAbs, av);
            s.maxAbs = std::max(s.maxAbs, av);
            absSum += av;
            ++s.count;
        }

        if (std::isfinite(r.derivNorm[metricIndex]))
        {
            derivNormSum += r.derivNorm[metricIndex];
            s.maxDerivNorm = std::max(s.maxDerivNorm, r.derivNorm[metricIndex]);
        }
        if (std::isfinite(r.derivMeanAbs[metricIndex]))
            derivMeanAbsSum += r.derivMeanAbs[metricIndex];
        if (std::isfinite(r.derivAbsRange[metricIndex]))
            derivAbsRangeSum += r.derivAbsRange[metricIndex];
    }

    if (s.count == 0)
    {
        s.minVal = s.maxVal = 0.0;
        s.minAbs = s.maxAbs = 0.0;
        return s;
    }

    s.range = s.maxVal - s.minVal;
    s.absRange = s.maxAbs - s.minAbs;
    s.meanAbs = absSum / static_cast<double>(s.count);
    s.meanDerivNorm = derivNormSum / static_cast<double>(s.count);
    s.meanDerivMeanAbs = derivMeanAbsSum / static_cast<double>(s.count);
    s.meanDerivAbsRange = derivAbsRangeSum / static_cast<double>(s.count);
    return s;
}

static double ChooseValueScale(const MetricStats& s)
{
    if (s.absRange > 1.0e-12) return s.absRange;
    if (s.range > 1.0e-12) return s.range;
    if (s.maxAbs > 1.0e-12) return s.maxAbs;
    if (s.meanAbs > 1.0e-12) return s.meanAbs;
    return 0.0;
}

static double ChooseDerivativeScale(const MetricStats& s)
{
    if (s.meanDerivNorm > 1.0e-12) return s.meanDerivNorm;
    if (s.meanDerivAbsRange > 1.0e-12) return s.meanDerivAbsRange;
    if (s.maxDerivNorm > 1.0e-12) return s.maxDerivNorm;
    return 0.0;
}

static std::vector<SummaryEntry> BuildSummary(
    const std::vector<SweepResult>& all,
    const CalibrationSelection& selection)
{
    std::vector<SummaryEntry> summary;
    summary.reserve(MetricCount);

    for (unsigned int i = 0; i < MetricCount; ++i)
    {
        SummaryEntry e;
        e.index = i;
        e.name = MetricNames[i];
        e.valueEnabled = selection.value[i];
        e.derivativeEnabled = selection.derivative[i];
        e.stats = ComputeMetricStats(all, i);
        e.valueScale = ChooseValueScale(e.stats);
        e.derivativeScale = ChooseDerivativeScale(e.stats);
        summary.push_back(e);
    }

    auto findReference = [&summary](bool derivative) -> double {
        const unsigned int preferred = MetricMI;
        const SummaryEntry& mi = summary[preferred];
        const bool miEnabled = derivative ? mi.derivativeEnabled : mi.valueEnabled;
        const double miScale = derivative ? mi.derivativeScale : mi.valueScale;
        if (miEnabled && miScale > 1.0e-12)
            return miScale;

        for (const auto& e : summary)
        {
            const bool enabled = derivative ? e.derivativeEnabled : e.valueEnabled;
            const double scale = derivative ? e.derivativeScale : e.valueScale;
            if (enabled && scale > 1.0e-12)
                return scale;
        }
        return 0.0;
    };

    const double refValue = findReference(false);
    const double refDerivative = findReference(true);

    for (auto& e : summary)
    {
        if (e.valueEnabled && refValue > 1.0e-12 && e.valueScale > 1.0e-12)
            e.relValueWeight = refValue / e.valueScale;
        if (e.derivativeEnabled && refDerivative > 1.0e-12 && e.derivativeScale > 1.0e-12)
            e.relDerivativeWeight = refDerivative / e.derivativeScale;
    }

    return summary;
}

static const SummaryEntry* FindSummaryEntry(
    const std::vector<SummaryEntry>& summary,
    unsigned int metricIndex)
{
    for (const auto& e : summary)
        if (e.index == metricIndex)
            return &e;
    return nullptr;
}

// ─────────────────────────────────────────────────────────────────────────────
// Write CSV
// ─────────────────────────────────────────────────────────────────────────────
static void WriteCSV(const std::string& path,
                     const std::vector<SweepResult>& all,
                     const std::vector<SummaryEntry>& summary)
{
    std::ofstream f(path);
    if (!f.is_open())
    {
        std::cerr << "[Calibrate] Cannot open CSV output file: " << path << "\n";
        return;
    }
    f << "axis,param_value";
    for (unsigned int i = 0; i < MetricCount; ++i)
    {
        const std::string& n = MetricNames[i];
        f << "," << n
          << "," << n << "_abs"
          << "," << n << "_range"
          << "," << n << "_abs_range"
          << "," << n << "_rel_value_weight"
          << "," << n << "_deriv_norm"
          << "," << n << "_deriv_mean_abs"
          << "," << n << "_deriv_abs_range"
          << "," << n << "_rel_deriv_weight";
    }
    f << ",Total,Total_abs,TotalGradNorm\n";

    f << std::fixed << std::setprecision(8);
    for (const auto& r : all)
    {
        f << r.axis << "," << r.paramValue << ",";
        for (unsigned int i = 0; i < MetricCount; ++i)
        {
            const SummaryEntry* e = FindSummaryEntry(summary, i);
            const double range = e ? e->stats.range : 0.0;
            const double absRange = e ? e->stats.absRange : 0.0;
            const double relValue = e ? e->relValueWeight : 0.0;
            const double relDerivative = e ? e->relDerivativeWeight : 0.0;
            f << r.value[i] << ","
              << r.absValue[i] << ","
              << range << ","
              << absRange << ","
              << relValue << ","
              << r.derivNorm[i] << ","
              << r.derivMeanAbs[i] << ","
              << r.derivAbsRange[i] << ","
              << relDerivative;
            f << ",";
        }
        f << r.total << "," << r.absTotal << "," << r.gradNorm << "\n";
    }
    std::cout << "[Calibrate] CSV written to: " << path << "\n";
}

// ─────────────────────────────────────────────────────────────────────────────
// Print summary table and relative value/derivative weights
// ─────────────────────────────────────────────────────────────────────────────
static void PrintSummary(
    const std::vector<SweepResult>& all,
    const CalibrationSelection& selection,
    bool computeDerivative)
{
    const auto summary = BuildSummary(all, selection);

    std::cout << "\n";
    std::cout << std::string(160, '=') << "\n";
    std::cout << "  CALIBRATION SUMMARY (unit-weight intrinsic metric scales)\n";
    std::cout << std::string(160, '=') << "\n";
    std::cout << std::left
              << std::setw(8)  << "Metric"
              << std::setw(8)  << "Value?"
              << std::setw(14) << "Min"
              << std::setw(14) << "Max"
              << std::setw(14) << "Range"
              << std::setw(14) << "AbsRange"
              << std::setw(14) << "MeanAbs"
              << std::setw(14) << "RelValueW";
    if (computeDerivative)
        std::cout << std::setw(8)  << "Deriv?"
                  << std::setw(14) << "MeanDerNorm"
                  << std::setw(14) << "MeanDerAbs"
                  << std::setw(14) << "DerAbsRange"
                  << std::setw(14) << "RelDerivW";
    std::cout << "\n";
    std::cout << std::string(160, '-') << "\n";

    std::cout << std::fixed << std::setprecision(5);

    for (const auto& e : summary)
    {
        std::cout << std::setw(8)  << e.name
                  << std::setw(8)  << (e.valueEnabled ? "yes" : "no")
                  << std::setw(14) << e.stats.minVal
                  << std::setw(14) << e.stats.maxVal
                  << std::setw(14) << e.stats.range
                  << std::setw(14) << e.stats.absRange
                  << std::setw(14) << e.stats.meanAbs
                  << std::setw(14) << (e.valueEnabled ? e.relValueWeight : 0.0);
        if (computeDerivative)
            std::cout << std::setw(8)  << (e.derivativeEnabled ? "yes" : "no")
                      << std::setw(14) << e.stats.meanDerivNorm
                      << std::setw(14) << e.stats.meanDerivMeanAbs
                      << std::setw(14) << e.stats.meanDerivAbsRange
                      << std::setw(14) << (e.derivativeEnabled ? e.relDerivativeWeight : 0.0);
        std::cout << "\n";
    }

    std::cout << std::string(160, '=') << "\n";
    std::cout << "  RelValueW uses each metric's absolute value scale, preferably AbsRange.\n"
              << "  RelDerivW uses each metric derivative's mean L2 norm. Reference is MI\n"
              << "  when MI is active; otherwise the first active non-zero metric is used.\n"
              << "  Use RelValueW for --alpha/--lambda/--nu/--rho/--yota/--sigma and\n"
              << "  RelDerivW for the corresponding derivative weights.\n";
    std::cout << std::string(160, '=') << "\n";
}

// ─────────────────────────────────────────────────────────────────────────────
// main
// ─────────────────────────────────────────────────────────────────────────────
int main(int argc, char* argv[])
{
    po::options_description desc(
        "3DRegCalibrate — Multi-metric calibration tool\n"
        "Sweeps a rigid/similarity perturbation and measures each sub-metric's\n"
        "unit-weight value and derivative scale to suggest balanced ratios.\n\n"
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

        // Calibration selector arrays. Values are masks, not final weights:
        // any non-zero value enables that metric at unit calibration scale.
        ("metrics", po::value<std::string>()->default_value(""),
             "Metric selector array: MI,NGF,MSE,GD,NC,NMI[,Label]; non-zero = calibrate")
        ("metric-derivatives", po::value<std::string>()->default_value(""),
             "Derivative selector array: MI,NGF,MSE,GD,NC,NMI[,Label]; defaults to --metrics")
        ("metric-sampling", po::value<std::string>()->default_value(""),
             "Metric sampling fractions: MI,NGF,MSE,GD,NC,NMI[,Label]")

        // Legacy individual selectors (same flag names as registration tools).
        // Numeric values are not calibration weights; non-zero means selected.
        ("alpha",          po::value<double>()->default_value(1.0),  "MI value selector")
        ("alphaderivative",po::value<double>()->default_value(1.0),  "MI derivative selector")
        ("mattesnumberofbins,b", po::value<int>()->default_value(64),"Mattes bins")
        ("mattespercentage,p", po::value<double>()->default_value(0.05), "Mattes sampling fraction")
        ("lambda",         po::value<double>()->default_value(0.0),  "NGF value selector")
        ("lambdaderivative",po::value<double>()->default_value(0.0), "NGF derivative selector")
        ("ngfpercentage",  po::value<double>()->default_value(0.05), "NGF sampling fraction")
        ("NGFevaluator",   po::value<int>()->default_value(0),       "NGF evaluator (0=scalar)")
        ("ngfprecompute",  po::value<bool>()->default_value(false),
             "Precompute moving-image NGF once and resample vector field each iteration")
        ("etavaluefixed",  po::value<double>()->default_value(-1),   "NGF eta fixed (-1=auto)")
        ("etavaluemoving", po::value<double>()->default_value(-1),   "NGF eta moving (-1=auto)")
        ("nu",             po::value<double>()->default_value(0.0),  "MSE value selector")
        ("nuderivative",   po::value<double>()->default_value(0.0),  "MSE derivative selector")
        ("msepercentage",  po::value<double>()->default_value(0.05), "MSE sampling fraction")
        ("rho",            po::value<double>()->default_value(0.0),  "GD value selector")
        ("rhoderivative",  po::value<double>()->default_value(0.0),  "GD derivative selector")
        ("gdpercentage",   po::value<double>()->default_value(0.05), "GD sampling fraction")
        ("yota",           po::value<double>()->default_value(0.0),  "NC value selector")
        ("yotaderivative", po::value<double>()->default_value(0.0),  "NC derivative selector")
        ("ncpercentage",   po::value<double>()->default_value(0.05), "NC sampling fraction")
        ("sigma",          po::value<double>()->default_value(0.0),  "NMI value selector")
        ("sigmaderivative",po::value<double>()->default_value(0.0),  "NMI derivative selector")
        ("nmibins",        po::value<int>()->default_value(64),      "NMI bins")
        ("nmipercentage",  po::value<double>()->default_value(0.05), "NMI sampling fraction")
        ("normalizemse",   po::value<bool>()->default_value(false),  "Normalise MSE")
        ("normalizegd",    po::value<bool>()->default_value(false),  "Normalise GD")
        ("labelkappa",     po::value<double>()->default_value(0.0),  "Label value selector")
        ("labelkappaderiv",po::value<double>()->default_value(0.0),  "Label derivative selector")
        ("labelkappavec",  po::value<std::string>()->default_value(""), "Per-label kappa weights")
        ("labelkappaderivvec", po::value<std::string>()->default_value(""), "Per-label kappa deriv")
        ("labelsamples",   po::value<double>()->default_value(0.05), "Label sampling fraction")
        ("labeldistmax",   po::value<double>()->default_value(20.0), "Label distance clamp (mm)")
        ("labelnarrowband",po::value<bool>()->default_value(false),
             "Use narrow-band label loss near boundaries only")
        ("labelbandwidth", po::value<double>()->default_value(5.0), "Narrow-band half-width in mm")
        ("labelhuber",     po::value<bool>()->default_value(false),
             "Use Huber loss for normalized label residuals")
        ("labelhuberdelta",po::value<double>()->default_value(0.25), "Huber delta")

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
             "Also compute per-metric derivative norms (slower; requires selected derivatives)")

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
    try
    {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
    }
    catch (const std::exception& e)
    {
        std::cerr << "[Calibrate] CLI parse error: " << e.what() << "\n";
        std::cerr << "Use --help to list supported options.\n";
        return EXIT_FAILURE;
    }

    if (vm.count("help") || !vm.count("fixedimage") || !vm.count("movingimage"))
    {
        std::cout << desc << "\n";
        return 1;
    }

    const bool doDerivative = vm["compute-derivative"].as<bool>();
    CalibrationSelection selection;
    std::array<double, MetricCount> sampling;
    try
    {
        selection = BuildSelection(vm);
        sampling = BuildSampling(vm);
    }
    catch (const std::exception& e)
    {
        std::cerr << "[Calibrate] CLI parse error: " << e.what() << "\n";
        return EXIT_FAILURE;
    }
    PrintSelection(selection);

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

    // Calibrate intrinsic scales: selected metrics are evaluated with unit
    // weights. The output summary suggests the real registration weights.
    metric->SetAlpha(UnitIf(selection.value[MetricMI]));
    metric->SetAlphaDerivative(UnitIf(doDerivative && selection.derivative[MetricMI]));
    metric->SetLambda(UnitIf(selection.value[MetricNGF]));
    metric->SetLambdaDerivative(UnitIf(doDerivative && selection.derivative[MetricNGF]));
    metric->SetNu(UnitIf(selection.value[MetricMSE]));
    metric->SetNuDerivative(UnitIf(doDerivative && selection.derivative[MetricMSE]));
    metric->SetRho(UnitIf(selection.value[MetricGD]));
    metric->SetRhoDerivative(UnitIf(doDerivative && selection.derivative[MetricGD]));
    metric->SetYota(UnitIf(selection.value[MetricNC]));
    metric->SetYotaDerivative(UnitIf(doDerivative && selection.derivative[MetricNC]));
    metric->SetSigma(UnitIf(selection.value[MetricNMI]));
    metric->SetSigmaDerivative(UnitIf(doDerivative && selection.derivative[MetricNMI]));
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
    metric->SetNGFPrecomputeGradient(vm["ngfprecompute"].as<bool>());

    // Sampling counts
    const unsigned int nPix = fixedImage->GetLargestPossibleRegion().GetNumberOfPixels();
    metric->SetMANumberOfSamples(static_cast<unsigned int>(nPix * sampling[MetricMI]));
    metric->SetBinNumbers(vm["mattesnumberofbins"].as<int>());
    metric->SetNGFNumberOfSamples(static_cast<unsigned int>(nPix * sampling[MetricNGF]));
    metric->SetMSENumberOfSamples(static_cast<unsigned int>(nPix * sampling[MetricMSE]));
    metric->SetGDNumberOfSamples(static_cast<unsigned int>(nPix * sampling[MetricGD]));
    metric->SetNCNumberOfSamples(static_cast<unsigned int>(nPix * sampling[MetricNC]));
    metric->SetNMIBinNumbers(vm["nmibins"].as<int>());
    metric->SetNMINumberOfSamples(static_cast<unsigned int>(nPix * sampling[MetricNMI]));
    metric->SetLabelNumberOfSamples(RegCommon::ResolveLabelSampleCount(
        sampling[MetricLabel], nPix));
    metric->SetLabelDistanceMax(vm["labeldistmax"].as<double>());
    metric->SetLabelUseNarrowBand(vm["labelnarrowband"].as<bool>());
    metric->SetLabelNarrowBandWidth(vm["labelbandwidth"].as<double>());
    metric->SetLabelUseHuber(vm["labelhuber"].as<bool>());
    metric->SetLabelHuberDelta(vm["labelhuberdelta"].as<double>());

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
        metric->SetLabelKappa(UnitIf(selection.value[MetricLabel]));
        metric->SetLabelKappaDerivative(UnitIf(doDerivative && selection.derivative[MetricLabel]));
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
    const auto summary = BuildSummary(allResults, selection);

    // ── CSV output ───────────────────────────────────────────────────────────
    const std::string csvPath = vm["output-csv"].as<std::string>();
    if (!csvPath.empty())
        WriteCSV(csvPath, allResults, summary);

    // ── Print summary ─────────────────────────────────────────────────────────
    PrintSummary(allResults, selection, doDerivative);

    return EXIT_SUCCESS;
}
