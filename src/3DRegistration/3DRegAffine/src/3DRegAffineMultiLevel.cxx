#include "itkMultiResolutionImageRegistrationMethod.h"
#include "../../../Metrics/Mplus/itkMplus.h"
#include "itkBSplineTransform.h"
#include "itkLBFGSBOptimizer.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkAffineTransform.h"
#include "itkIdentityTransform.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkImageMaskSpatialObject.h"

#include "itkTransformToDeformationFieldSource.h"
#include "itkTransformFileWriter.h"
#include "itkTransformFileReader.h"
#include "../../Version.h"

#include "../../../Metrics/NGF/NGFImageMetric/NGFImageToImageMetric/Code/itkGetImageNoiseFunction.h"
#include "../../../includes/imageUtils.h"
#include "../../../includes/registrationUtils.h"
#include "../../../includes/RegistrationCommon.h"

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"

#include "itkCenteredTransformInitializer.h"
#include "itkRegularStepGradientDescentOptimizer.h"

#include <boost/program_options.hpp>
#include <sstream>
#include <map>
namespace po = boost::program_options;

#include <algorithm>      // for std::replace
#include <iterator>       // for std::istream_iterator
#include <sstream>        // for std::istringstream

const unsigned int ImageDimension = 3;

	typedef  float          PixelType;
	typedef itk::Image< PixelType, ImageDimension >  FixedImageType;

	const unsigned int SpaceDimension = ImageDimension;
	const unsigned int SplineOrder = 3;
	typedef double CoordinateRepType;
	
	typedef itk::AffineTransform<double, 3> TransformType;


	typedef itk::Mplus<
			FixedImageType,
			FixedImageType >    MetricType;

	typedef itk:: LinearInterpolateImageFunction<
			FixedImageType,
			double          >    InterpolatorType;

	typedef itk::MultiResolutionImageRegistrationMethod<
			FixedImageType,
			FixedImageType >    RegistrationType;

typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
int main( int argc, char *argv[] )
{
    po::options_description desc("B-spline Registration\n"
	"Dr. Eros Montin Ph.D., 2014\n"
	"eros.montin@gmail.com\n\n"
	"cite us:\n\nMontin, E., Belfatto, A., Bologna, M., Meroni, S., Cavatorta, C., Pecori, E., Diletto, B., Massimino, M., Oprandi, M. C., Poggi, G., Arrigoni, F., Peruzzo, D., Pignoli, E., Gandola, L., Cerveri, P., & Mainardi, L. (2020). A multi-metric registration strategy for the alignment of longitudinal brain images in pediatric oncology. Medical & biological engineering & computing, 58(4), 843–855. https://doi.org/10.1007/s11517-019-02109-4\n\n"
	"Allowed options for alpha MI + lambda NGF +  nu MSE + yota NMI ");
    std::string method;
    int NL=2;
	double YOTA=0.1;
	double YOTADERIVATIVE=0;
	desc.add_options()
	    ("help,h", "produce help message")
        ("fixedimage,f", po::value<std::string>(), "Fixed image filename")
        ("movingimage,m", po::value<std::string>(), "Moving image filename")
        ("outputimage,o", po::value<std::string>(), "Output registered imagefilename")
        ("vfout,v", po::value<std::string>()->default_value("N"), "VF output filename")
        ("numberofthreads", po::value<int>()->default_value(2), "Number of threads 2")
	    ("alpha,a", po::value<double>()->default_value(1.0), "alpha value MI 1.0")
		("alphaderivative,A", po::value<double>()->default_value(1.0), "alpha derivative MI 1.0")
		("subtype", po::value<std::string>(&method)->default_value("affine"), "Subtype (translation, rotation, scaling, affine)")
		("mattespercentage,p", po::value<double>()->default_value(0.1), "Mattes percentage 0.1")
        ("mattesnumberofbins,b", po::value<int>()->default_value(64), "Mattes number of bins 64")
		("explicitPDFderivatives", po::value<bool>()->default_value(false), "Explicit PDF derivatives, 0 for false")
		("lambda,l", po::value<double>()->default_value(1.0), "lambda value NGF 1.0")
        ("lambdaderivative,L", po::value<double>()->default_value(0), "Lambda derivative NGF 0 no derivatives")
        ("yota,y", po::value<double>(&YOTA)->default_value(0), "Yota value NC (Normalized Correlation) weight")
        ("yotaderivative,Y", po::value<double>(&YOTADERIVATIVE)->default_value(0), "Yota derivative NC, 0 = no derivatives") 
        ("etavaluefixed,r", po::value<double>()->default_value(-1), "Eta value fixed image(NGF noise) -1 (autodetermine)")
        ("etavaluemoving,s", po::value<double>()->default_value(-1), "Eta value moving image (NGF noise) -1 (autodetermine)")
        ("NGFevaluator", po::value<int>()->default_value(0), "NGF Evaluator (0 scalar,1cross,2scdelta,3Delta,4Delta2)")
	    ("nu,n", po::value<double>()->default_value(1.0), "nu value MSE 1.0")
		("nuderivative,N", po::value<double>()->default_value(1.0), "nu MSE derivative 1.0")
        ("maxnumberofiterations,I", po::value<int>()->default_value(1000), "Max number of Iterations 1000")
		("minimumsteplength,S", po::value<double>()->default_value(0.1), "Minimum step length")
		("maximumsteplength,X", po::value<double>()->default_value(1.0), "Maximum step length")
		("relaxationfactor,R", po::value<double>()->default_value(0.5), "Relaxation factor")
		("gradientmagnitudetolerance,G", po::value<double>()->default_value(1e-4), "Gradient magnitude tolerance")
		("fixedimagethreshold,t", po::value<double>()->default_value(-99999999), "Fixed image threshold")
        ("transformout,T", po::value<std::string>()->default_value("N"), "Output for transform")
        ("transformin,W", po::value<std::string>()->default_value("N"), "Input no rigid transform for transform")
		("dfltpixelvalue,P", po::value<double>()->default_value(0), "Default pixel value")
		("verbose,V", po::value<bool>()->default_value(false), "verbose")
		("derivativemode", po::value<int>()->default_value(0), "Derivative merge mode: 0=consistent, 1=normalized (RSGD), 2=main-metric adaptive")
		("mainmetric", po::value<int>()->default_value(0), "Main metric index for mode 2: 0=MI, 1=NGF, 2=MSE, 3=NC, 4=Label, 5=GD, 6=NMI")
        ("numberoflevels,U", po::value<int>(&NL)->default_value(2), "Number of levels")
		("msepercentage",   po::value<double>()->default_value(0.1), "MSE percentage of pixels used (0.1 = 10%)")
		("ngfpercentage",   po::value<double>()->default_value(0.1), "NGF percentage of pixels used (0.1 = 10%)")
		("gdpercentage",    po::value<double>()->default_value(0.1), "GD percentage of pixels used (0.1 = 10%)")
		("nmipercentage",   po::value<double>()->default_value(0.1), "NMI percentage of pixels used (0.1 = 10%)")
		("ncpercentage",    po::value<double>()->default_value(0.1), "NC percentage of pixels used (0.1 = 10%)")
        ("nmipercentage",   po::value<double>()->default_value(0.1), "NMI percentage of pixels used (0.1 = 10%)")
        ("mipercentage",    po::value<double>()->default_value(0.1), "Histogram‐MI percentage of pixels used (0.1 = 10%)")
        ("rho",             po::value<double>()->default_value(0.0), "Rho weight for Gradient Difference (GD)")
        ("rhoderivative",   po::value<double>()->default_value(0.0), "Rho derivative for GD")
        ("sigma",           po::value<double>()->default_value(0.0), "Sigma weight for Normalized Mutual Information (NMI)")
        ("sigmaderivative", po::value<double>()->default_value(0.0), "Sigma derivative for NMI")
        ("nmibins",         po::value<int>()->default_value(64),     "Number of histogram bins for NMI")
        ("ngfspacing",      po::value<std::string>()->default_value("4,4,4"), "NGF spacing per dimension (x,y,z)")
        ("workingresolution", po::value<std::string>()->default_value("0,0,0"),
            "Internal registration spacing in mm (x,y,z). Use 0,0,0 to keep the input spacing.")
        ("ngfprecompute", po::value<bool>()->default_value(false), "Precompute moving-image NGF once and resample vector field each iteration (faster, approximate)")
	("metricoverlap", po::value<bool>()->default_value(true), "Compute overlap between fixed and moving image (default true)")
	("fixedlabelmap",  po::value<std::string>()->default_value("N"), "Fixed label map filename (N = none)")
	("movinglabelmap", po::value<std::string>()->default_value("N"), "Moving label map filename (N = none)")
	("labelkappa",     po::value<double>()->default_value(0.0),       "Global kappa weight for label-map distance metric (0 = off)")
	("labelkappaderiv",po::value<double>()->default_value(0.0),       "Global kappa weight for label-map derivative")
	("labelkappavec",  po::value<std::string>()->default_value(""),   "Per-label kappa weights: 'L1:w1,L2:w2,...'")
	("labelkappaderivvec", po::value<std::string>()->default_value(""),"Per-label kappa derivative weights")
	("labelsamples",   po::value<double>()->default_value(0.1), "Label metric percentage of pixels used (0.1 = 10%)")
	("labelreport",    po::value<int>()->default_value(1),            "Report Dice every N iterations (0 = off)")
	("snapshotdir",    po::value<std::string>()->default_value("N"), "Directory for iteration snapshots (N = off)")
	("snapshotevery",  po::value<int>()->default_value(1),            "Save snapshot every N iterations")
	("snapshotstack",  po::value<bool>()->default_value(false),       "Save full 3D .nii.gz instead of mid-slice PNG")
	("snapshotgrid",   po::value<bool>()->default_value(true),        "Overlay warped grid on snapshot panels (default on)")
	("snapshotgridspacing", po::value<unsigned int>()->default_value(20), "Grid line spacing in voxels")
	("snapshotlinewidth", po::value<double>()->default_value(0.0),    "Overlay line width in pixels (0 = auto)")
	("version", "Print version and exit")
	("overlappadding", po::value<unsigned int>()->default_value(20), "Overlap padding in voxels")
	("modality", po::value<std::string>()->default_value("custom"),
		"Preset modality: 'multimodal' (MI+NGF), 'singlemodal' (MSE+NC), or 'custom' (manual weights)")
    ;
	


    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

	// ── Handle --help and --version BEFORE any file I/O ────────────────────────
	if (vm.count("version")) {
		std::cout << "3DRegAffineMultiLevel v" << VERSION_MAJOR << "." << VERSION_MINOR << std::endl;
		return EXIT_SUCCESS;
	}
	if (vm.count("help") || !vm.count("fixedimage") || !vm.count("movingimage") || !vm.count("outputimage")) {
		std::cout << desc << "\n";
		return 1;
	}

	// ── Modality presets ───────────────────────────────────────────────────────
	// Apply sensible defaults based on --modality BEFORE reading individual weights.
	// Any weight explicitly supplied on the command line will override the preset.
	const std::string MODALITY = vm["modality"].as<std::string>();
	if (MODALITY == "multimodal") {
		// Multimodal: MI dominates, NGF adds structural guidance, MSE/NC off
		if (vm["alpha"].defaulted())            const_cast<po::variable_value&>(vm["alpha"]).value()            = 1.0;
		if (vm["alphaderivative"].defaulted())  const_cast<po::variable_value&>(vm["alphaderivative"]).value()  = 1.0;
		if (vm["lambda"].defaulted())           const_cast<po::variable_value&>(vm["lambda"]).value()           = 0.5;
		if (vm["lambdaderivative"].defaulted()) const_cast<po::variable_value&>(vm["lambdaderivative"]).value() = 0.5;
		if (vm["nu"].defaulted())               const_cast<po::variable_value&>(vm["nu"]).value()               = 0.0;
		if (vm["nuderivative"].defaulted())     const_cast<po::variable_value&>(vm["nuderivative"]).value()     = 0.0;
		if (vm["yota"].defaulted())             const_cast<po::variable_value&>(vm["yota"]).value()             = 0.0;
		if (vm["yotaderivative"].defaulted())   const_cast<po::variable_value&>(vm["yotaderivative"]).value()   = 0.0;
		if (vm["rho"].defaulted())              const_cast<po::variable_value&>(vm["rho"]).value()              = 0.0;
		if (vm["rhoderivative"].defaulted())    const_cast<po::variable_value&>(vm["rhoderivative"]).value()    = 0.0;
		if (vm["sigma"].defaulted())            const_cast<po::variable_value&>(vm["sigma"]).value()            = 0.0;
		if (vm["sigmaderivative"].defaulted())  const_cast<po::variable_value&>(vm["sigmaderivative"]).value()  = 0.0;
		std::cout << "[Modality] multimodal preset: MI(alpha=" << vm["alpha"].as<double>()
		          << ") + NGF(lambda=" << vm["lambda"].as<double>() << ")" << std::endl;
	} else if (MODALITY == "singlemodal") {
		// Singlemodal: MSE + NC complement each other, MI/NGF off
		if (vm["alpha"].defaulted())            const_cast<po::variable_value&>(vm["alpha"]).value()            = 0.0;
		if (vm["alphaderivative"].defaulted())  const_cast<po::variable_value&>(vm["alphaderivative"]).value()  = 0.0;
		if (vm["lambda"].defaulted())           const_cast<po::variable_value&>(vm["lambda"]).value()           = 0.0;
		if (vm["lambdaderivative"].defaulted()) const_cast<po::variable_value&>(vm["lambdaderivative"]).value() = 0.0;
		if (vm["nu"].defaulted())               const_cast<po::variable_value&>(vm["nu"]).value()               = 1.0;
		if (vm["nuderivative"].defaulted())     const_cast<po::variable_value&>(vm["nuderivative"]).value()     = 1.0;
		if (vm["yota"].defaulted())             const_cast<po::variable_value&>(vm["yota"]).value()             = 0.5;
		if (vm["yotaderivative"].defaulted())   const_cast<po::variable_value&>(vm["yotaderivative"]).value()   = 0.5;
		if (vm["rho"].defaulted())              const_cast<po::variable_value&>(vm["rho"]).value()              = 0.0;
		if (vm["rhoderivative"].defaulted())    const_cast<po::variable_value&>(vm["rhoderivative"]).value()    = 0.0;
		if (vm["sigma"].defaulted())            const_cast<po::variable_value&>(vm["sigma"]).value()            = 0.0;
		if (vm["sigmaderivative"].defaulted())  const_cast<po::variable_value&>(vm["sigmaderivative"]).value()  = 0.0;
		std::cout << "[Modality] singlemodal preset: MSE(nu=" << vm["nu"].as<double>()
		          << ") + NC(yota=" << vm["yota"].as<double>() << ")" << std::endl;
	} else if (MODALITY != "custom") {
		std::cerr << "Error: --modality must be 'multimodal', 'singlemodal', or 'custom', got '" << MODALITY << "'" << std::endl;
		return EXIT_FAILURE;
	}

	MetricType::Pointer         metric        = MetricType::New();
	InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
	RegistrationType::Pointer   registration  = RegistrationType::New();

	int  DERIVMODE = vm["derivativemode"].as<int>();
	int  MAINMETRIC = vm["mainmetric"].as<int>();


	// ── label map options
	const std::string FIXEDLABELMAP   = vm["fixedlabelmap"].as<std::string>();
	const std::string MOVINGLABELMAP  = vm["movinglabelmap"].as<std::string>();
	const double      LABELKAPPA      = vm["labelkappa"].as<double>();
	const double      LABELKAPPADERIV = vm["labelkappaderiv"].as<double>();
	const double LABELSAMPLES   = vm["labelsamples"].as<double>();
	const int         LABELREPORT     = vm["labelreport"].as<int>();
	const std::string SNAPSHOTDIR     = vm["snapshotdir"].as<std::string>();
	const int         SNAPSHOTEVERY   = vm["snapshotevery"].as<int>();
	const bool        SNAPSHOTSTACK   = vm["snapshotstack"].as<bool>();
	const bool        SNAPSHOTGRID    = vm["snapshotgrid"].as<bool>();
	const unsigned int SNAPSHOTGRIDSP  = vm["snapshotgridspacing"].as<unsigned int>();
	const double      SNAPSHOTLINEW   = vm["snapshotlinewidth"].as<double>();
	const auto LABELKAPPAVEC      = RegCommon::ParseLabelWeights(vm["labelkappavec"].as<std::string>());
	const auto LABELKAPPADERIVVEC = RegCommon::ParseLabelWeights(vm["labelkappaderivvec"].as<std::string>());
	typedef itk::Image<short, ImageDimension> LabelImageType;
	LabelImageType::ConstPointer fixedLabelMap, movingLabelMap;
	if (FIXEDLABELMAP != "N" && MOVINGLABELMAP != "N") {
		typedef itk::ImageFileReader<LabelImageType> LR;
		auto flr=LR::New(); flr->SetFileName(FIXEDLABELMAP); flr->Update();
					{ LabelImageType::Pointer _tmp = flr->GetOutput(); _tmp->DisconnectPipeline(); fixedLabelMap = _tmp; }
					auto mlr=LR::New(); mlr->SetFileName(MOVINGLABELMAP); mlr->Update();
					{ LabelImageType::Pointer _tmp = mlr->GetOutput(); _tmp->DisconnectPipeline(); movingLabelMap = _tmp; }
		std::cout << "[Label] Fixed: " << FIXEDLABELMAP << " Moving: " << MOVINGLABELMAP << std::endl;
	}

	if (vm["verbose"].as<bool>())
		RegCommon::PrintOptions(vm);

    std::string fixedImageFN = vm["fixedimage"].as<std::string>();
    std::string movingImageFN = vm["movingimage"].as<std::string>();

	std::string ou = vm["outputimage"].as<std::string>();
	std::string VOUT = vm["vfout"].as<std::string>();

	int NT=vm["numberofthreads"].as<int>();
	int NB=vm["mattesnumberofbins"].as<int>();
    double MAPERCENTAGE = vm["mattespercentage"].as<double>();
	double ALPHA = vm["alpha"].as<double>();
	double ALPHADERIVATIVE = vm["alphaderivative"].as<double>();
	double NU = vm["nu"].as<double>();
	double NUDERIVATIVE = vm["nuderivative"].as<double>();
    double LAMBDA = vm["lambda"].as<double>();
    double LAMBDADERIVATIVE = vm["lambdaderivative"].as<double>();
    double ETAF = vm["etavaluefixed"].as<double>();
    double ETAM = vm["etavaluemoving"].as<double>();
	int NGFevaluator = vm["NGFevaluator"].as<int>();
	if (NGFevaluator < 0 || NGFevaluator > 4) {
    std::cerr << "Error: NGFevaluator must be between 0 and 4" << std::endl;
    return EXIT_FAILURE;
}
	int NI=vm["maxnumberofiterations"].as<int>();
	double TR=vm["fixedimagethreshold"].as<double>();

	bool EPDF=vm["explicitPDFderivatives"].as<bool>();

	std::string TOUT=vm["transformout"].as<std::string>();
	std::string TIN=vm["transformin"].as<std::string>();
	double DFLTPIXELVALUE=vm["dfltpixelvalue"].as<double>();
	double MINSTEP=vm["minimumsteplength"].as<double>();
	double MAXSTEP=vm["maximumsteplength"].as<double>();
	double RF=vm["relaxationfactor"].as<double>();
	double GMT=vm["gradientmagnitudetolerance"].as<double>();

    double MSEPERCENTAGE = vm["msepercentage"].as<double>();
    double NGFPERCENTAGE = vm["ngfpercentage"].as<double>();
    // double NMIPERCENTAGE = vm["nmipercentage"].as<double>();  // NMI not yet in Mplus
    double NCPERCENTAGE  = vm["ncpercentage"].as<double>();
    double RHO             = vm["rho"].as<double>();
    double RHODERIVATIVE   = vm["rhoderivative"].as<double>();
    double SIGMA           = vm["sigma"].as<double>();
    double SIGMADERIVATIVE = vm["sigmaderivative"].as<double>();
    int    NMIBINS         = vm["nmibins"].as<int>();

    // parse ngfspacing string → SpacingType
    {
        auto s = vm["ngfspacing"].as<std::string>();
        std::replace(s.begin(), s.end(), ',', ' ');
        std::istringstream iss(s);
        std::vector<double> tmp{
            std::istream_iterator<double>(iss),
            std::istream_iterator<double>()};
        if (tmp.size() != ImageDimension)
        {
            std::cerr << "Error: ngfspacing must have "
                      << ImageDimension << " comma-separated values\n";
            return EXIT_FAILURE;
        }
        FixedImageType::SpacingType ngf;
        for (unsigned i = 0; i < ImageDimension; ++i)
            ngf[i] = tmp[i];
        vm.insert({ "parsed_ngfspacing",
            po::variable_value(boost::any(ngf), false) });
    }

    FixedImageType::SpacingType workingSpacing;
    bool useWorkingResolution = false;
    if (!RegCommon::ParseOptionalSpacing<FixedImageType::SpacingType, ImageDimension>(
            vm["workingresolution"].as<std::string>(),
            "workingresolution",
            workingSpacing,
            useWorkingResolution))
    {
        return EXIT_FAILURE;
    }

	auto optimizer = OptimizerType::New();


	registration->SetMetric(        metric        );
	registration->SetOptimizer(     optimizer     );
	registration->SetInterpolator(  interpolator  );
	registration->SetNumberOfThreads(NT);
	registration->SetFixedImageRegion( FixedImageType::RegionType() );

	TransformType::Pointer  transform = TransformType::New();
	registration->SetTransform( transform );

	
	
	typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
	typedef itk::ImageFileReader< FixedImageType > MovingImageReaderType;

	FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
	MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

	fixedImageReader->SetFileName(  fixedImageFN );
	movingImageReader->SetFileName( movingImageFN );

	fixedImageReader->Update();
	movingImageReader->Update();

	FixedImageType::ConstPointer originalFixedImage = fixedImageReader->GetOutput();
	FixedImageType::ConstPointer originalMovingImage = movingImageReader->GetOutput();
	FixedImageType::ConstPointer fixedImage = originalFixedImage;
	FixedImageType::ConstPointer movingImage = originalMovingImage;

	// ── Input validation ──────────────────────────────────────────────────────
	if (!originalFixedImage || originalFixedImage->GetLargestPossibleRegion().GetNumberOfPixels() == 0) {
		std::cerr << "Error: Failed to load fixed image or image is empty." << std::endl;
		return EXIT_FAILURE;
	}
	if (!originalMovingImage || originalMovingImage->GetLargestPossibleRegion().GetNumberOfPixels() == 0) {
		std::cerr << "Error: Failed to load moving image or image is empty." << std::endl;
		return EXIT_FAILURE;
	}

	if (useWorkingResolution) {
		std::cout << "[WorkingResolution] Resampling registration inputs to spacing "
		          << workingSpacing << std::endl;
		if (!RegCommon::SpacingEquals(originalFixedImage->GetSpacing(), workingSpacing)) {
			fixedImage = RegCommon::ResampleScalarImageToSpacing<FixedImageType>(
			    originalFixedImage, workingSpacing);
		}
		if (!RegCommon::SpacingEquals(originalMovingImage->GetSpacing(), workingSpacing)) {
			movingImage = RegCommon::ResampleScalarImageToSpacing<FixedImageType>(
			    originalMovingImage, workingSpacing, DFLTPIXELVALUE);
		}
		if (fixedLabelMap && !RegCommon::SpacingEquals(fixedLabelMap->GetSpacing(), workingSpacing)) {
			fixedLabelMap = RegCommon::ResampleNearestNeighborImageToSpacing<LabelImageType>(
			    fixedLabelMap, workingSpacing);
		}
		if (movingLabelMap && !RegCommon::SpacingEquals(movingLabelMap->GetSpacing(), workingSpacing)) {
			movingLabelMap = RegCommon::ResampleNearestNeighborImageToSpacing<LabelImageType>(
			    movingLabelMap, workingSpacing);
		}
	}

  using FixedImagePyramidType =
    itk::MultiResolutionPyramidImageFilter<FixedImageType,
                                           FixedImageType>;
  using MovingImagePyramidType =
    itk::MultiResolutionPyramidImageFilter<FixedImageType,
                                           FixedImageType>;
 
  auto fixedImagePyramid = FixedImagePyramidType::New();
  auto movingImagePyramid = MovingImagePyramidType::New();

	registration->SetFixedImage(  fixedImage   );
    registration->SetFixedImagePyramid( fixedImagePyramid);
	registration->SetMovingImage(   movingImage);
    registration->SetMovingImagePyramid( movingImagePyramid);

	FixedImageType::RegionType fixedRegion = fixedImage->GetBufferedRegion();
	registration->SetFixedImageRegion( fixedRegion );

	
	
	
// #let's fix a few things
if ((LAMBDA!=0) || (LAMBDADERIVATIVE!=0))
{

	if ((ETAF==-1) || (ETAM==-1))
	{
		metric->SetAutoEstimateEta(true);
	}
}

	
// Duplicate CenteredTransformInitializer removed — kept in TIN/else block below

	  using OptimizerScalesType = OptimizerType::ScalesType;
  OptimizerScalesType optimizerScales(
    transform->GetNumberOfParameters());
	optimizerScales.Fill(1.0);
  // Translation parameters (indices 9-11) need a small scale so the optimizer
  // takes proportionally larger steps for them (gradient /= scale).
  const double translationScale = 1.0 / 1000.0;
  optimizerScales[9] = translationScale;
  optimizerScales[10] = translationScale;
  optimizerScales[11] = translationScale;
  // Large value used only to genuinely freeze parameters in subtype modes
  constexpr double kFreezeScale = 1e10;
  optimizer->SetNumberOfIterations(NI);
  optimizer->SetMinimumStepLength(MINSTEP);
  optimizer->SetRelaxationFactor(RF);
  optimizer->SetGradientMagnitudeTolerance(GMT);
  optimizer->SetMaximumStepLength(MAXSTEP);



if (method == "translation") {
    // Only optimize translation parameters — freeze matrix entries
    optimizerScales.Fill(1.0);
    for (int i = 0; i < 9; ++i) {
        optimizerScales[i] = kFreezeScale;
    }
} else if (method == "rotation") {
    // Only optimize rotation parameters
    optimizerScales.Fill(1.0);
    for (int i = 3; i < 12; ++i) {
        if (i < 9) {
            optimizerScales[i] = 1.0;
        } else {
            optimizerScales[i] = kFreezeScale;
        }
    }
} else if (method == "scaling") {
    // Only optimize scaling parameters
    optimizerScales.Fill(1.0);
    for (int i = 0; i < 12; ++i) {
        if (i == 0 || i == 4 || i == 8) {
            optimizerScales[i] = 1.0;
        } else {
            optimizerScales[i] = kFreezeScale;
        }
    }
}

  optimizer->SetScales(optimizerScales);

	typedef TransformType::ParametersType     ParametersType;

	const unsigned int numberOfParameters =
			transform->GetNumberOfParameters();

	ParametersType parametersLow( numberOfParameters );

	parametersLow.Fill( 0.0 );

	transform->SetParameters( parametersLow );

	registration->SetInitialTransformParameters( transform->GetParameters() );

	const unsigned int numberOfPixels = fixedImage->GetLargestPossibleRegion().GetNumberOfPixels();
	const unsigned int numberOfSamplesMA =static_cast< unsigned int >( numberOfPixels * MAPERCENTAGE );
    const unsigned int numberOfSamplesMSE = static_cast<unsigned int>(numberOfPixels * MSEPERCENTAGE);
    const unsigned int numberOfSamplesNGF = static_cast<unsigned int>(numberOfPixels * NGFPERCENTAGE);
    // const unsigned int numberOfSamplesNMI = static_cast<unsigned int>(numberOfPixels * NMIPERCENTAGE);
    // const unsigned int numberOfSamplesHMI = static_cast<unsigned int>(numberOfPixels * MIPERCENTAGE);
    const unsigned int numberOfSamplesNC  = static_cast<unsigned int>(numberOfPixels * NCPERCENTAGE);

	metric->SetAlpha(ALPHA);
	metric->SetAlphaDerivative(ALPHADERIVATIVE);
	metric->SetMSENumberOfSamples(numberOfSamplesMSE);
	metric->SetBinNumbers(NB);
	metric->SetMANumberOfSamples(numberOfSamplesMA);
	metric->SetUseExplicitPDFDerivatives(EPDF);
	metric->SetNu(NU);
	metric->SetNuDerivative(NUDERIVATIVE);
	metric->SetLambda(LAMBDA);
	metric->SetLambdaDerivative(LAMBDADERIVATIVE);
	metric->SetNGFNumberOfSamples(numberOfSamplesNGF);
	metric->SetDerivativeMode(DERIVMODE);
	metric->SetMainMetricIndex(MAINMETRIC);
	metric->SetComputeOverlap(vm["metricoverlap"].as<bool>());
	metric->SetOverlapPadding(vm["overlappadding"].as<unsigned int>());
	metric->SetNumberOfThreads(NT);
	metric->SetYota(YOTA);
	metric->SetYotaDerivative(YOTADERIVATIVE);
    // metric->SetNMINumberOfSamples(numberOfSamplesNMI);  // NMI not yet implemented in Mplus
    metric->SetNCNumberOfSamples(numberOfSamplesNC);
    metric->SetRho(RHO);
    metric->SetRhoDerivative(RHODERIVATIVE);
    metric->SetSigma(SIGMA);
    metric->SetSigmaDerivative(SIGMADERIVATIVE);
    metric->SetNMIBinNumbers(NMIBINS);
	
	metric->SetFixedEta(ETAF);
	metric->SetMovingEta(ETAM);
	metric->SetEvaluator(NGFevaluator);
	
	if (TR!=-99999999)
	{
		metric->SetFixedImageThreshold(TR);
	}

    // apply the parsed spacing
    {
        auto ngf = boost::any_cast<FixedImageType::SpacingType>(vm["parsed_ngfspacing"].value());
        metric->SetNGFSpacing(ngf);
    }
    metric->SetNGFPrecomputeGradient(vm["ngfprecompute"].as<bool>());

	if (fixedLabelMap && movingLabelMap) {
		metric->SetFixedLabelMap(fixedLabelMap);
		metric->SetMovingLabelMap(movingLabelMap);
		metric->SetLabelKappa(LABELKAPPA);
		metric->SetLabelKappaDerivative(LABELKAPPADERIV);
		metric->SetLabelNumberOfSamples(LABELSAMPLES);
		if (!LABELKAPPAVEC.empty())      metric->SetLabelKappaWeights(LABELKAPPAVEC);
		if (!LABELKAPPADERIVVEC.empty()) metric->SetLabelKappaDerivativeWeights(LABELKAPPADERIVVEC);
	}

	if (TIN!="N")
	{
		typedef itk::TransformFileReader TransformReaderType;
		TransformReaderType::Pointer transformReader = TransformReaderType::New();
		transformReader->SetFileName( TIN );
		transformReader->Update();

		auto *rawPtr = transformReader->GetTransformList()->front().GetPointer();
		TransformType *directCast = dynamic_cast<TransformType*>(rawPtr);

		if (directCast)
		{
			transform = directCast;
		}
		else
		{
			// The file may contain a different linear transform (e.g. Similarity3DTransform).
			// Convert it to AffineTransform via the common MatrixOffsetTransformBase interface.
			typedef itk::MatrixOffsetTransformBase<double, 3, 3> MatrixTransformBaseType;
			MatrixTransformBaseType *matBase = dynamic_cast<MatrixTransformBaseType*>(rawPtr);
			if (matBase)
			{
				std::cout << "[Info] Input transform is '" << rawPtr->GetTransformTypeAsString()
				          << "' — converting to AffineTransform." << std::endl;
				transform->SetMatrix(matBase->GetMatrix());
				transform->SetOffset(matBase->GetOffset());
				transform->SetCenter(matBase->GetCenter());
			}
			else
			{
				std::cerr << "Error: Transform read from '" << TIN
				          << "' is of type '" << rawPtr->GetTransformTypeAsString()
				          << "' and cannot be converted to AffineTransform." << std::endl;
				return EXIT_FAILURE;
			}
		}

	}else
	{
		transform->SetIdentity();
		// allign the center of the images
		typedef itk::CenteredTransformInitializer<
				TransformType,
				FixedImageType,
				FixedImageType >  TransformInitializerType;

		TransformInitializerType::Pointer initializer = TransformInitializerType::New();

		initializer->SetTransform(   transform );
		initializer->SetFixedImage(  fixedImage );
		initializer->SetMovingImage( movingImage );
		initializer->MomentsOn();

		initializer->InitializeTransform();


	}
		registration->SetInitialTransformParameters( transform->GetParameters() );
	std::cout << "Starting Registration "
			<< std::endl;

	// SaveImage<movingImageType>(ApplyTransform<FixedImageType, TransformType>(movingImageReader->GetOutput(), transform));
	
	std::cout<< "\n\n\n\n Affine Transform using itkMplus	\n\tThread: "<< metric->GetNumberOfThreads() <<
			"\nVariables: " <<transform->GetNumberOfParameters() <<

			std::endl;

			TransformType::ParametersType init_ = transform->GetParameters();
			std::cout << "Initial transform parameters: " << init_ << std::endl;

			init_ = registration->GetInitialTransformParameters();
			std::cout << "Initial transform parameters rec " << init_ << std::endl;

      // Add a time probe
    itk::TimeProbesCollectorBase   chronometer;
    itk::MemoryProbesCollectorBase memorymeter;




  using CommandType = RegistrationInterfaceCommand<RegistrationType>;
  auto command = CommandType::New();
  registration->AddObserver(itk::IterationEvent(), command);
  registration->SetNumberOfLevels(NL);
  // ── label Dice monitoring (attaches to the optimizer, not the multi-res registration)
  if (fixedLabelMap && movingLabelMap && LABELREPORT > 0) {
	  auto optimizer = static_cast<OptimizerType*>(registration->GetModifiableOptimizer());
	  LabelMapDiceObserver<TransformType, LabelImageType>::Pointer lo =
	      LabelMapDiceObserver<TransformType, LabelImageType>::New();
	  lo->SetFixedLabelMap(fixedLabelMap); lo->SetMovingLabelMap(movingLabelMap);
	  lo->SetTransform(transform);
	  lo->SetEvaluateEveryNIterations(static_cast<unsigned int>(LABELREPORT));
	  optimizer->AddObserver(itk::IterationEvent(), lo);
  }
  using SnapObsType = IterationSnapshotObserver<TransformType, FixedImageType>;
  SnapObsType::Pointer snapObs;
  if (SNAPSHOTDIR != "N") {
	  snapObs = SnapObsType::New();
	  snapObs->SetFixedImage(fixedImage);
	  snapObs->SetMovingImage(movingImage);
	  snapObs->SetTransform(transform);
	  snapObs->SetOutputDirectory(SNAPSHOTDIR);
	  snapObs->SetSaveEveryNIterations(static_cast<unsigned int>(SNAPSHOTEVERY));
	  snapObs->SetSaveStack(SNAPSHOTSTACK);
	  snapObs->SetShowDeformationGrid(SNAPSHOTGRID);
	  snapObs->SetGridSpacingPixels(SNAPSHOTGRIDSP);
	  snapObs->SetOverlayLineWidthPixels(SNAPSHOTLINEW);
	  snapObs->SetMetricValuesGetter([metric]() -> std::map<std::string,double> {
		  return {
			  {"Total", metric->GetLastValTotal()},
			  {"MI",    metric->GetLastValMI()},
			  {"NGF",   metric->GetLastValNGF()},
			  {"MSE",   metric->GetLastValMSE()},
			  {"NC",    metric->GetLastValNC()},
			  {"GD",    metric->GetLastValGD()},
			  {"NMI",   metric->GetLastValNMI()},
			  {"Label", metric->GetLastValLabel()},
		  };
	  });
	  optimizer->AddObserver(itk::IterationEvent(), snapObs);
  }


	try
	{
    memorymeter.Start("Registration");
    chronometer.Start("Registration");
 
    registration->Update();
 
    chronometer.Stop("Registration");
    memorymeter.Stop("Registration");
 
    std::cout << "Optimizer stop condition = "
              << registration->GetOptimizer()->GetStopConditionDescription()
              << std::endl;
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}

	if (snapObs) snapObs->FinalizeConvergencePlot();

  // Report the time and memory taken by the registration
  chronometer.Report(std::cout);
  memorymeter.Report(std::cout);
  
	transform->SetParameters( registration->GetLastTransformParameters() );

	typedef itk::ResampleImageFilter<
			FixedImageType,
			FixedImageType >    ResampleFilterType;

	ResampleFilterType::Pointer resample = ResampleFilterType::New();

	resample->SetTransform( transform );
	resample->SetInput( movingImageReader->GetOutput() );

	resample->SetSize(    originalFixedImage->GetLargestPossibleRegion().GetSize() );
	resample->SetOutputOrigin(  originalFixedImage->GetOrigin() );
	resample->SetOutputSpacing( originalFixedImage->GetSpacing() );
	resample->SetOutputDirection( originalFixedImage->GetDirection() );
	resample->SetDefaultPixelValue( DFLTPIXELVALUE );

	typedef itk::ImageFileWriter< FixedImageType >  WriterType;
	WriterType::Pointer      writer =  WriterType::New();


	writer->SetFileName( ou );
	writer->SetInput( resample->GetOutput()   );


	try
	{
		writer->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}


	if (VOUT!="N")
	{
		typedef itk::Vector< float,  ImageDimension >  VectorType2;
		typedef itk::Image< VectorType2,  ImageDimension >   OutputTransformationImageType2;
		typedef itk::TransformToDeformationFieldSource< OutputTransformationImageType2, double >TransformToDeformationFieldSourceType2;
		TransformToDeformationFieldSourceType2::Pointer td = TransformToDeformationFieldSourceType2::New();
		td->SetOutputParametersFromImage(movingImageReader->GetOutput());
		td->SetTransform( registration->GetOutput()->Get() );
		//std::cout<<registration->GetOutput()->Get()<<std::endl;
		typedef itk::ImageFileWriter< OutputTransformationImageType2>TransformToDeformationFieldSourceWriterType;
		TransformToDeformationFieldSourceWriterType::Pointer rtd = TransformToDeformationFieldSourceWriterType::New();
		rtd->SetInput(td->GetOutput());
		rtd->SetFileName(VOUT);
		rtd->Update();
	};


	if (TOUT!="N")
			{
			#if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 5) || ITK_VERSION_MAJOR > 4
				itk::TransformFileWriterTemplate<double>::Pointer writer =
						itk::TransformFileWriterTemplate<double>::New();
	#else
				itk::TransformFileWriter::Pointer writer = itk::TransformFileWriter::New();
	#endif
				writer->SetInput(registration->GetOutput()->Get());
				writer->SetFileName(TOUT);
				writer->Update();

};


	return EXIT_SUCCESS;
}
