#pragma once
// RegistrationCommon.h — Shared utility functions for registration programs
// Reduces code duplication across 3DRegAffine, 3DRegSimilarity, 3DRegAffineMultiLevel, 3DRegBsplines

#include <boost/program_options.hpp>
#include <itkIdentityTransform.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkResampleImageFilter.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <map>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

namespace RegCommon {

namespace po = boost::program_options;

// ── Parse "L1:w1,L2:w2,..." into a map ─────────────────────────────────────
inline std::map<short, double> ParseLabelWeights(const std::string & s)
{
    std::map<short, double> result;
    if (s.empty()) return result;
    std::istringstream ss(s);
    std::string token;
    while (std::getline(ss, token, ','))
    {
        const auto pos = token.find(':');
        if (pos != std::string::npos)
            result[static_cast<short>(std::stoi(token.substr(0, pos)))] =
                std::stod(token.substr(pos + 1));
    }
    return result;
}

// ── Convert label sample input to an absolute sample count ───────────────────
// Accepts either:
//   - fraction in (0, 1]  -> interpreted as percentage of image voxels
//   - absolute count > 1  -> interpreted as direct number of samples
inline unsigned int ResolveLabelSampleCount(double rawLabelSamples,
                                            std::size_t numberOfPixels)
{
    const std::size_t safePixelCount = std::max<std::size_t>(1, numberOfPixels);
    const std::size_t defaultCount =
        std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(0.1 * safePixelCount)));
    const auto clampToUInt = [](std::size_t v) -> unsigned int {
        constexpr std::size_t kMaxU = std::numeric_limits<unsigned int>::max();
        return static_cast<unsigned int>(std::min(v, kMaxU));
    };

    if (!std::isfinite(rawLabelSamples) || rawLabelSamples <= 0.0)
        return clampToUInt(defaultCount);

    if (rawLabelSamples <= 1.0)
    {
        const double asCount = std::ceil(rawLabelSamples * static_cast<double>(safePixelCount));
        const std::size_t v = std::max<std::size_t>(1, static_cast<std::size_t>(asCount));
        return clampToUInt(v);
    }

    const double capped = std::min(rawLabelSamples, static_cast<double>(safePixelCount));
    const std::size_t v = std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(capped)));
    return clampToUInt(v);
}

// ── Print all boost::program_options ────────────────────────────────────────
inline void PrintOptions(const po::variables_map & vm)
{
    for (const auto & it : vm)
    {
        std::cout << "Option: " << it.first.c_str() << "\nValue: ";
        auto & value = it.second.value();

        if (value.type() == typeid(int))
            std::cout << *boost::any_cast<int>(&value);
        else if (value.type() == typeid(double))
            std::cout << *boost::any_cast<double>(&value);
        else if (value.type() == typeid(std::string))
            std::cout << *boost::any_cast<std::string>(&value);
        else if (value.type() == typeid(bool))
            std::cout << std::boolalpha << *boost::any_cast<bool>(&value);
        else if (value.type() == typeid(unsigned int))
            std::cout << *boost::any_cast<unsigned int>(&value);
        else
            std::cout << "Unknown type";
        std::cout << "\n-------------------\n";
    }
}

// ── Validate a loaded image is non-null and non-empty ───────────────────────
template <typename TImage>
inline bool ValidateImage(const typename TImage::Pointer & img, const std::string & name)
{
    if (!img || img->GetLargestPossibleRegion().GetNumberOfPixels() == 0)
    {
        std::cerr << "Error: Failed to load " << name << " or image is empty." << std::endl;
        return false;
    }
    return true;
}

template <typename TImage>
inline bool ValidateImage(const typename TImage::ConstPointer & img, const std::string & name)
{
    if (!img || img->GetLargestPossibleRegion().GetNumberOfPixels() == 0)
    {
        std::cerr << "Error: Failed to load " << name << " or image is empty." << std::endl;
        return false;
    }
    return true;
}

// ── Common CLI options shared by all registration programs ──────────────────
inline void AddCommonOptions(po::options_description & desc,
                             std::string & /*unused method ref kept for compat*/)
{
    desc.add_options()
        ("version", "Print version and exit")
        ("overlappadding", po::value<unsigned int>()->default_value(20),
            "Overlap padding in voxels (default 20)")
    ;
}

// ── Handle --help and --version early, before any heavy I/O ─────────────────
inline bool HandleHelpAndVersion(const po::variables_map & vm,
                                 const po::options_description & desc,
                                 int versionMajor, int versionMinor)
{
    if (vm.count("version"))
    {
        std::cout << "Registration Suite v" << versionMajor << "." << versionMinor << std::endl;
        return true; // caller should exit
    }

    if (vm.count("help") || !vm.count("fixedimage") ||
        !vm.count("movingimage") || !vm.count("outputimage"))
    {
        std::cout << desc << "\n";
        return true; // caller should exit
    }

    return false; // continue
}

// ── Wire up the Mplus metric from CLI variables ─────────────────────────────
template <typename MetricPointer, typename VariablesMap, typename SpacingType>
inline void WireMetric(MetricPointer & metric, const VariablesMap & vm,
                       const SpacingType & ngfSpacing,
                       unsigned int numberOfPixels)
{
    const double MAPERCENTAGE  = vm["mattespercentage"].template as<double>();
    const double MSEPERCENTAGE = vm["msepercentage"].template as<double>();
    const double NGFPERCENTAGE = vm["ngfpercentage"].template as<double>();
    const double NCPERCENTAGE  = vm["ncpercentage"].template as<double>();
    const double GDPERCENTAGE  = vm["gdpercentage"].template as<double>();
    const double NMIPERCENTAGE = vm["nmipercentage"].template as<double>();

    metric->SetComputeOverlap(vm["metricoverlap"].template as<bool>());
    metric->SetOverlapPadding(vm["overlappadding"].template as<unsigned int>());
    metric->SetUseExplicitPDFDerivatives(vm["explicitPDFderivatives"].template as<bool>());
    metric->SetNumberOfThreads(vm["numberofthreads"].template as<int>());

    const int derivMode = vm["derivativemode"].template as<int>();
    const int mainMetric = vm["mainmetric"].template as<int>();
    metric->SetDerivativeMode(derivMode);
    metric->SetMainMetricIndex(mainMetric);

    metric->SetAlpha(vm["alpha"].template as<double>());
    metric->SetAlphaDerivative(vm["alphaderivative"].template as<double>());
    metric->SetMANumberOfSamples(static_cast<unsigned int>(numberOfPixels * MAPERCENTAGE));
    metric->SetBinNumbers(vm["mattesnumberofbins"].template as<int>());

    metric->SetFixedEta(vm["etavaluefixed"].template as<double>());
    metric->SetMovingEta(vm["etavaluemoving"].template as<double>());
    metric->SetEvaluator(vm["NGFevaluator"].template as<int>());

    metric->SetLambda(vm["lambda"].template as<double>());
    metric->SetLambdaDerivative(vm["lambdaderivative"].template as<double>());
    metric->SetNGFNumberOfSamples(static_cast<unsigned int>(numberOfPixels * NGFPERCENTAGE));
    metric->SetNGFSpacing(ngfSpacing);

    metric->SetMSENumberOfSamples(static_cast<unsigned int>(numberOfPixels * MSEPERCENTAGE));
    metric->SetNu(vm["nu"].template as<double>());
    metric->SetNuDerivative(vm["nuderivative"].template as<double>());

    metric->SetYota(vm["yota"].template as<double>());
    metric->SetYotaDerivative(vm["yotaderivative"].template as<double>());
    metric->SetNCNumberOfSamples(static_cast<unsigned int>(numberOfPixels * NCPERCENTAGE));

    metric->SetRho(vm["rho"].template as<double>());
    metric->SetRhoDerivative(vm["rhoderivative"].template as<double>());
    metric->SetGDNumberOfSamples(static_cast<unsigned int>(numberOfPixels * GDPERCENTAGE));

    metric->SetSigma(vm["sigma"].template as<double>());
    metric->SetSigmaDerivative(vm["sigmaderivative"].template as<double>());
    metric->SetNMIBinNumbers(vm["nmibins"].template as<int>());
    metric->SetNMINumberOfSamples(static_cast<unsigned int>(numberOfPixels * NMIPERCENTAGE));

    // Threshold
    const double TR = vm["fixedimagethreshold"].template as<double>();
    if (TR != -99999999)
        metric->SetFixedImageThreshold(TR);
}

// ── Read label maps and wire them to the metric ─────────────────────────────
template <typename MetricPointer, typename LabelImageType, unsigned int Dim>
inline void WireLabelMaps(MetricPointer & metric, const po::variables_map & vm)
{
    const std::string fixedLabelFN  = vm["fixedlabelmap"].template as<std::string>();
    const std::string movingLabelFN = vm["movinglabelmap"].template as<std::string>();

    if (fixedLabelFN == "N" || movingLabelFN == "N") return;

    typedef itk::ImageFileReader<LabelImageType> LabelReaderType;
    auto flr = LabelReaderType::New(); flr->SetFileName(fixedLabelFN); flr->Update();
    typename LabelImageType::Pointer fixTmp = flr->GetOutput(); fixTmp->DisconnectPipeline();
    auto mlr = LabelReaderType::New(); mlr->SetFileName(movingLabelFN); mlr->Update();
    typename LabelImageType::Pointer movTmp = mlr->GetOutput(); movTmp->DisconnectPipeline();

    std::cout << "[Label] Fixed: " << fixedLabelFN << " Moving: " << movingLabelFN << std::endl;

    metric->SetFixedLabelMap(fixTmp);
    metric->SetMovingLabelMap(movTmp);
    metric->SetLabelKappa(vm["labelkappa"].template as<double>());
    metric->SetLabelKappaDerivative(vm["labelkappaderiv"].template as<double>());
    // Use a 64-bit pixel count so huge label maps don't truncate.
    const std::size_t numberOfPixels =
        static_cast<std::size_t>(fixTmp->GetLargestPossibleRegion().GetNumberOfPixels());
    metric->SetLabelNumberOfSamples(
        ResolveLabelSampleCount(vm["labelsamples"].template as<double>(), numberOfPixels));
    if (vm.count("labeldistmax"))
        metric->SetLabelDistanceMax(vm["labeldistmax"].template as<double>());
    if (vm.count("labelnarrowband"))
        metric->SetLabelUseNarrowBand(vm["labelnarrowband"].template as<bool>());
    if (vm.count("labelbandwidth"))
        metric->SetLabelNarrowBandWidth(vm["labelbandwidth"].template as<double>());
    if (vm.count("labelhuber"))
        metric->SetLabelUseHuber(vm["labelhuber"].template as<bool>());
    if (vm.count("labelhuberdelta"))
        metric->SetLabelHuberDelta(vm["labelhuberdelta"].template as<double>());

    const auto kappaVec      = ParseLabelWeights(vm["labelkappavec"].template as<std::string>());
    const auto kappaDerivVec = ParseLabelWeights(vm["labelkappaderivvec"].template as<std::string>());
    if (!kappaVec.empty())      metric->SetLabelKappaWeights(kappaVec);
    if (!kappaDerivVec.empty()) metric->SetLabelKappaDerivativeWeights(kappaDerivVec);
}

// ── Parse NGF spacing string → SpacingType ──────────────────────────────────
template <typename SpacingType, unsigned int Dim>
inline bool ParseNGFSpacing(const std::string & raw, SpacingType & out)
{
    std::string s = raw;
    std::replace(s.begin(), s.end(), ',', ' ');
    std::istringstream iss(s);
    std::vector<double> v{std::istream_iterator<double>(iss),
                          std::istream_iterator<double>()};
    if (v.size() != Dim)
    {
        std::cerr << "Error: ngfspacing must have " << Dim
                  << " comma-separated values, got " << v.size() << std::endl;
        return false;
    }
    for (unsigned i = 0; i < Dim; ++i) out[i] = v[i];
    return true;
}

// ── Parse optional working resolution: "sx,sy,sz" or disabled via "0,0,0" ──
template <typename SpacingType, unsigned int Dim>
inline bool ParseOptionalSpacing(const std::string & raw,
                                 const std::string & optionName,
                                 SpacingType & out,
                                 bool & enabled)
{
    std::string s = raw;
    std::replace(s.begin(), s.end(), ',', ' ');
    std::istringstream iss(s);
    std::vector<double> v{std::istream_iterator<double>(iss),
                          std::istream_iterator<double>()};

    if (v.size() != Dim)
    {
        std::cerr << "Error: " << optionName << " must have " << Dim
                  << " comma-separated values, got " << v.size() << std::endl;
        return false;
    }

    const bool allZero = std::all_of(v.begin(), v.end(),
        [](const double value) { return value == 0.0; });

    if (allZero)
    {
        enabled = false;
        return true;
    }

    for (const double value : v)
    {
        if (value <= 0.0)
        {
            std::cerr << "Error: " << optionName
                      << " values must all be > 0, or all be 0 to disable."
                      << std::endl;
            return false;
        }
    }

    for (unsigned int i = 0; i < Dim; ++i)
    {
        out[i] = v[i];
    }
    enabled = true;
    return true;
}

template <typename SpacingType>
inline bool SpacingEquals(const SpacingType & lhs, const SpacingType & rhs,
                          const double tolerance = 1e-6)
{
    for (unsigned int i = 0; i < lhs.Size(); ++i)
    {
        if (std::abs(lhs[i] - rhs[i]) > tolerance)
        {
            return false;
        }
    }
    return true;
}

template <typename TImage>
inline typename TImage::SizeType ComputeSizeForSpacing(
    const typename TImage::ConstPointer & image,
    const typename TImage::SpacingType & outputSpacing)
{
    typename TImage::SizeType outputSize;
    const typename TImage::SizeType inputSize =
        image->GetLargestPossibleRegion().GetSize();
    const typename TImage::SpacingType inputSpacing = image->GetSpacing();

    for (unsigned int i = 0; i < TImage::ImageDimension; ++i)
    {
        if (inputSize[i] <= 1)
        {
            outputSize[i] = 1;
            continue;
        }

        const double physicalExtent =
            (static_cast<double>(inputSize[i]) - 1.0) * inputSpacing[i];
        const double scaledSize = std::llround(physicalExtent / outputSpacing[i]) + 1.0;
        outputSize[i] = static_cast<typename TImage::SizeType::SizeValueType>(
            std::max(1.0, scaledSize));
    }

    return outputSize;
}

template <typename TImage, typename TInterpolator>
inline typename TImage::Pointer ResampleImageToSpacing(
    const typename TImage::ConstPointer & image,
    const typename TImage::SpacingType & outputSpacing,
    const typename TInterpolator::Pointer & interpolator,
    const double defaultPixelValue = 0.0)
{
    typedef itk::IdentityTransform<double, TImage::ImageDimension> IdentityTransformType;
    typedef itk::ResampleImageFilter<TImage, TImage> ResampleFilterType;

    typename IdentityTransformType::Pointer identity = IdentityTransformType::New();
    identity->SetIdentity();

    typename ResampleFilterType::Pointer resample = ResampleFilterType::New();
    resample->SetInput(image);
    resample->SetTransform(identity);
    resample->SetInterpolator(interpolator);
    resample->SetSize(ComputeSizeForSpacing<TImage>(image, outputSpacing));
    resample->SetOutputOrigin(image->GetOrigin());
    resample->SetOutputSpacing(outputSpacing);
    resample->SetOutputDirection(image->GetDirection());
    resample->SetDefaultPixelValue(static_cast<typename TImage::PixelType>(defaultPixelValue));
    resample->Update();

    typename TImage::Pointer output = resample->GetOutput();
    output->DisconnectPipeline();
    return output;
}

template <typename TImage>
inline typename TImage::Pointer ResampleScalarImageToSpacing(
    const typename TImage::ConstPointer & image,
    const typename TImage::SpacingType & outputSpacing,
    const double defaultPixelValue = 0.0)
{
    typedef itk::LinearInterpolateImageFunction<TImage, double> InterpolatorType;
    typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
    return ResampleImageToSpacing<TImage, InterpolatorType>(
        image, outputSpacing, interpolator, defaultPixelValue);
}

template <typename TImage>
inline typename TImage::Pointer ResampleNearestNeighborImageToSpacing(
    const typename TImage::ConstPointer & image,
    const typename TImage::SpacingType & outputSpacing,
    const double defaultPixelValue = 0.0)
{
    typedef itk::NearestNeighborInterpolateImageFunction<TImage, double> InterpolatorType;
    typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
    return ResampleImageToSpacing<TImage, InterpolatorType>(
        image, outputSpacing, interpolator, defaultPixelValue);
}

} // namespace RegCommon
