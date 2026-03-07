#pragma once
// RegistrationCommon.h — Shared utility functions for registration programs
// Reduces code duplication across 3DRegAffine, 3DRegSimilarity, 3DRegAffineMultiLevel, 3DRegBsplines

#include <boost/program_options.hpp>
#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

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

    metric->SetSigma(vm["sigma"].template as<double>());
    metric->SetSigmaDerivative(vm["sigmaderivative"].template as<double>());
    metric->SetNMIBinNumbers(vm["nmibins"].template as<int>());

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
    metric->SetLabelNumberOfSamples(vm["labelsamples"].template as<unsigned int>());

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

} // namespace RegCommon
