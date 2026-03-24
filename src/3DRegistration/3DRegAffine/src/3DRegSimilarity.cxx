#include "itkImageRegistrationMethod.h"
#include "../../../Metrics/Mplus/itkMplus.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkAffineTransform.h"
#include "itkIdentityTransform.h"
#include "itkImageMaskSpatialObject.h"

#include "itkTransformToDeformationFieldSource.h"
#include "itkTransformFileWriter.h"
#include "itkTransformFileReader.h"
#include "../../Version.h"
#include "itkSimilarity3DTransform.h"
#include "../../../Metrics/NGF/NGFImageMetric/NGFImageToImageMetric/Code/itkGetImageNoiseFunction.h"
#include "../../../includes/imageUtils.h"
#include "../../../includes/registrationUtils.h"
#include "../../../includes/RegistrationCommon.h"

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"

#include "itkCenteredTransformInitializer.h"
#include "itkRegularStepGradientDescentOptimizer.h"

#include <boost/program_options.hpp>
#include <algorithm>      // for std::replace
#include <iterator>       // for std::istream_iterator
#include <sstream>        // for std::istringstream
#include <map>
namespace po = boost::program_options;

const unsigned int ImageDimension = 3;

	typedef  float          PixelType;
	typedef itk::Image< PixelType, ImageDimension >  FixedImageType;
	typedef itk::Image< PixelType, ImageDimension >  MovingImageType;

	const unsigned int SpaceDimension = ImageDimension;
	const unsigned int SplineOrder = 3;
	typedef double CoordinateRepType;
	
	typedef itk::Similarity3DTransform<double> TransformType;


	typedef itk::Mplus<
			FixedImageType,
			MovingImageType >    MetricType;

	typedef itk:: LinearInterpolateImageFunction<
			MovingImageType,
			double          >    InterpolatorType;

	typedef itk::ImageRegistrationMethod<
			FixedImageType,
			MovingImageType >    RegistrationType;

typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
int main( int argc, char *argv[] )
{
    po::options_description desc("B-spline Registration\n"
	"Dr. Eros Montin Ph.D., 2014\n"
	"eros.montin@gmail.com\n\n"
	"cite us:\n\nMontin, E., Belfatto, A., Bologna, M., Meroni, S., Cavatorta, C., Pecori, E., Diletto, B., Massimino, M., Oprandi, M. C., Poggi, G., Arrigoni, F., Peruzzo, D., Pignoli, E., Gandola, L., Cerveri, P., & Mainardi, L. (2020). A multi-metric registration strategy for the alignment of longitudinal brain images in pediatric oncology. Medical & biological engineering & computing, 58(4), 843–855. https://doi.org/10.1007/s11517-019-02109-4\n\n"
	"Allowed options for alpha MI + lambda NGF +  nu MSE +yota NMI\n\n");

	double YOTA=0;
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
		("mattespercentage,p", po::value<double>()->default_value(0.1), "Mattes percentage 0.1")
        ("mattesnumberofbins,b", po::value<int>()->default_value(64), "Mattes number of bins 64")
		("explicitPDFderivatives", po::value<bool>()->default_value(false), "Explicit PDF derivatives, 0 for false")
		("lambda,l", po::value<double>()->default_value(1.0), "lambda value NGF 1.0")
        ("lambdaderivative,L", po::value<double>()->default_value(0), "Lambda derivative NGF 0 no derivatives")
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
		// ("cir,c", po::value<boost::optional<std::vector<double>>>()->default_value(boost::none, ""), "CR 3D point index (x y z)")
		// ("CIR,C", po::value<std::vector<double>>()->multitoken()->default_value(boost::none, ""), "CR 3D index index (i j k)")
		("yota,y", po::value<double>(&YOTA)->default_value(0), "Yota value NC (Normalized Correlation) weight")
		("yotaderivative,Y", po::value<double>(&YOTADERIVATIVE)->default_value(0), "Yota derivative NC, 0 = no derivatives")        
		("msepercentage",   po::value<double>()->default_value(0.1), "MSE percentage of pixels used (0.1 = 10%)")
        ("ngfpercentage",   po::value<double>()->default_value(0.1), "NGF percentage of pixels used (0.1 = 10%)")
			("gdpercentage",    po::value<double>()->default_value(0.1), "GD percentage of pixels used (0.1 = 10%)")
			("nmipercentage",   po::value<double>()->default_value(0.1), "NMI percentage of pixels used (0.1 = 10%)")
		("ncpercentage", po::value<double>()->default_value(0.1), "NC percentage of pixels used (0.1 = 10%)")
		("rho", po::value<double>()->default_value(0.0), "Rho weight for Gradient Difference (GD)")
		("rhoderivative", po::value<double>()->default_value(0.0), "Rho derivative for GD")
		("sigma", po::value<double>()->default_value(0.0), "Sigma weight for Normalized Mutual Information (NMI)")
		("sigmaderivative", po::value<double>()->default_value(0.0), "Sigma derivative for NMI")
		("nmibins", po::value<int>()->default_value(64), "Number of histogram bins for NMI")
        ("ngfspacing",      po::value<std::string>()->default_value("4,4,4"),
                             "NGF spacing per dimension (x,y,z)")
        ("workingresolution", po::value<std::string>()->default_value("0,0,0"),
                             "Internal registration spacing in mm (x,y,z). Use 0,0,0 to keep the input spacing.")
        ("ngfprecompute", po::value<bool>()->default_value(false), "Precompute moving-image NGF once and resample vector field each iteration (faster, approximate)")
							 ("metricoverlap", po::value<bool>()->default_value(true), "Compute overlap between fixed and moving image (default true)")
	("fixedlabelmap",  po::value<std::string>()->default_value("N"), "Fixed label map filename (N = none)")
	("movinglabelmap", po::value<std::string>()->default_value("N"), "Moving label map filename (N = none)")
	("labelkappa",     po::value<double>()->default_value(0.0),       "Global kappa weight for label-map distance metric (0 = off)")
	("labelkappaderiv",po::value<double>()->default_value(0.0),       "Global kappa weight for label-map derivative")
	("labelkappavec",  po::value<std::string>()->default_value(""),   "Per-label kappa (value) weights: 'L1:w1,L2:w2,...'")
	("labelkappaderivvec", po::value<std::string>()->default_value(""),"Per-label kappa (derivative) weights: 'L1:w1,L2:w2,...'")
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
		std::cout << "3DRegSimilarity v" << VERSION_MAJOR << "." << VERSION_MINOR << std::endl;
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

	bool METRICOVERLAP = vm["metricoverlap"].as<bool>();
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

	// Parse NGF spacing
	auto ngfStr = vm["ngfspacing"].as<std::string>();
	std::replace(ngfStr.begin(), ngfStr.end(), ',', ' ');
	std::istringstream iss(ngfStr);
	std::vector<double> ngfVec((std::istream_iterator<double>(iss)),
							std::istream_iterator<double>());
	if (ngfVec.size() != ImageDimension) {
		std::cerr << "Error: --ngfspacing must have exactly " << ImageDimension
		          << " comma-separated values, got " << ngfVec.size() << std::endl;
		return EXIT_FAILURE;
	}
	FixedImageType::SpacingType ngf;
	for (unsigned i = 0; i < ImageDimension; ++i) ngf[i] = ngfVec[i];

	FixedImageType::SpacingType workingSpacing;
	bool useWorkingResolution = false;
	if (!RegCommon::ParseOptionalSpacing<FixedImageType::SpacingType, ImageDimension>(
	        vm["workingresolution"].as<std::string>(),
	        "workingresolution",
	        workingSpacing,
	        useWorkingResolution)) {
		return EXIT_FAILURE;
	}

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

	// bool iscir=false;
	// boost::optional<std::vector<double>> cir;
	// if (vm.count("cir")) {
	// 	iscir = true;
	// 	cir = vm["cir"].as<boost::optional<std::vector<double>>>();
	// 	if (cir && cir->size() != 3) {
	// 		std::cerr << "Error: cir must have exactly 3 coordinates" << std::endl;
	// 		return EXIT_FAILURE;
	// 	}
	// 	else if (vm.count("CIR"))
	// 	{

	// 		std::vector<double> CIR = vm["CIR"].as<std::vector<double>>();

	// 		if (CIR.size() != 3) {
	// 			std::cerr << "Error: CIR must have exactly 3 coordinates" << std::endl;
	// 			return EXIT_FAILURE;
	// 		}
			
	// 	}
		
	// }
   
	double NGFPERCENTAGE = vm["ngfpercentage"].as<double>();
	double MSEPERCENTAGE = vm["msepercentage"].as<double>();
	double NCPERCENTAGE = vm["ncpercentage"].as<double>();

	double RHO           = vm["rho"].as<double>();
	double RHODERIVATIVE = vm["rhoderivative"].as<double>();
	double SIGMA         = vm["sigma"].as<double>();
	double SIGMADERIVATIVE = vm["sigmaderivative"].as<double>();
	int    NMIBINS       = vm["nmibins"].as<int>();


	auto optimizer = OptimizerType::New();


	registration->SetMetric(        metric        );
	registration->SetOptimizer(     optimizer     );
	registration->SetInterpolator(  interpolator  );
	registration->SetNumberOfThreads(NT);
	registration->SetFixedImageRegion( FixedImageType::RegionType() );

	TransformType::Pointer  transform = TransformType::New();
	registration->SetTransform( transform );

	
	
	typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
	typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;

	FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
	MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

	fixedImageReader->SetFileName(  fixedImageFN );
	movingImageReader->SetFileName( movingImageFN );

	fixedImageReader->Update();
	movingImageReader->Update();

	FixedImageType::ConstPointer originalFixedImage = fixedImageReader->GetOutput();
	MovingImageType::ConstPointer originalMovingImage = movingImageReader->GetOutput();
	FixedImageType::ConstPointer fixedImage = originalFixedImage;
	MovingImageType::ConstPointer movingImage = originalMovingImage;

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
			movingImage = RegCommon::ResampleScalarImageToSpacing<MovingImageType>(
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

	registration->SetFixedImage(  fixedImage   );
	registration->SetMovingImage(   movingImage);

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

	
	transform->SetIdentity();
	// allign the center of the images
	typedef itk::CenteredTransformInitializer<
			TransformType,
			FixedImageType,
			MovingImageType >  TransformInitializerType;

	TransformInitializerType::Pointer initializer = TransformInitializerType::New();

	initializer->SetTransform(   transform );
	initializer->SetFixedImage(  fixedImage );
	initializer->SetMovingImage( movingImage );
	initializer->MomentsOn();
	initializer->InitializeTransform();



	  using OptimizerScalesType = OptimizerType::ScalesType;
  OptimizerScalesType optimizerScales(
    transform->GetNumberOfParameters());

// The serialization of the optimizable parameters is an array of 7 elements. 
// The first 3 elements are the components of the versor representation of 3D rotation. 
// The next 3 parameters defines the translation in each dimension. 
// The last parameter defines the isotropic scaling.

// The serialization of the fixed parameters is an array of 3 elements defining the center of rotation.

  optimizerScales.Fill(1.0);
  // Similarity3DTransform: [versor(3), translation(3), scale(1)]
  // Translation parameters (indices 3-5) need a small scale so the optimizer
  // takes proportionally larger steps for them (gradient /= scale).
  optimizerScales[3] = GMT;
  optimizerScales[4] = GMT;
  optimizerScales[5] = GMT;
  optimizer->SetNumberOfIterations(NI);
  optimizer->SetMinimumStepLength(MINSTEP);
  optimizer->SetRelaxationFactor(RF);
  optimizer->SetGradientMagnitudeTolerance(GMT);
  optimizer->SetMaximumStepLength(MAXSTEP);


  optimizer->SetScales(optimizerScales);


	typedef TransformType::ParametersType     ParametersType;

	const unsigned int numberOfParameters =
			transform->GetNumberOfParameters();

	ParametersType parametersLow( numberOfParameters );

	parametersLow.Fill( 0.0 );

	transform->SetParameters( parametersLow );


	// itk::Point<double, 3> center;
	// itk::Index<3> centerIndex;

	// for (int i = 0; i < 3; ++i) {
			// 	centerIndex[i] = sourceImage->GetLargestPossibleRegion().GetSize()[i] / 2;
			// }
			// sourceImage->TransformIndexToPhysicalPoint(centerIndex, center);
			// transform->SetCenter(center);



	registration->SetInitialTransformParameters( transform->GetParameters() );
	const unsigned int numberOfPixels = fixedImage->GetLargestPossibleRegion().GetNumberOfPixels();
	const unsigned int numberOfSamplesMA = static_cast<unsigned int>(numberOfPixels * MAPERCENTAGE);


	// const unsigned int numberOfSamplesGD =
	// 	static_cast<unsigned int>(numberOfPixels * GDPERCENTAGE);

	const unsigned int numberOfSamplesNGF = static_cast<unsigned int>(numberOfPixels * NGFPERCENTAGE);
	const unsigned int numberOfSamplesMSE = static_cast<unsigned int>(numberOfPixels * MSEPERCENTAGE);
	const unsigned int numberOfSamplesNC = static_cast<unsigned int>(numberOfPixels * NCPERCENTAGE);

	metric->SetComputeOverlap(METRICOVERLAP);
	metric->SetOverlapPadding(vm["overlappadding"].as<unsigned int>());
	metric->SetUseExplicitPDFDerivatives(EPDF);
	metric->SetNumberOfThreads(NT);
	metric->SetDerivativeMode(DERIVMODE);
	metric->SetMainMetricIndex(MAINMETRIC);

	metric->SetAlpha(ALPHA);
	metric->SetAlphaDerivative(ALPHADERIVATIVE);
	metric->SetMANumberOfSamples(numberOfSamplesMA);
	metric->SetBinNumbers(NB);
	metric->SetFixedEta(ETAF);
	metric->SetMovingEta(ETAM);
	metric->SetEvaluator(NGFevaluator);


	metric->SetLambda(LAMBDA);
	metric->SetLambdaDerivative(LAMBDADERIVATIVE);
	metric->SetNGFNumberOfSamples(numberOfSamplesNGF);
	metric->SetNGFSpacing(ngf);
	metric->SetNGFPrecomputeGradient(vm["ngfprecompute"].as<bool>());

	metric->SetMSENumberOfSamples(numberOfSamplesMSE);
	metric->SetNu(NU);
	metric->SetNuDerivative(NUDERIVATIVE);


	metric->SetYota(YOTA);
	metric->SetYotaDerivative(YOTADERIVATIVE);
	metric->SetNCNumberOfSamples(numberOfSamplesNC);

	metric->SetRho(RHO);
	metric->SetRhoDerivative(RHODERIVATIVE);

	metric->SetSigma(SIGMA);
	metric->SetSigmaDerivative(SIGMADERIVATIVE);
	metric->SetNMIBinNumbers(NMIBINS);

	if (TR!=-99999999)
	{
		metric->SetFixedImageThreshold(TR);
	}

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
		transform=dynamic_cast<TransformType*>(transformReader->GetTransformList()->front().GetPointer());

		
	}else
	{
		transform->SetIdentity();
		// allign the center of the images
		typedef itk::CenteredTransformInitializer<
				TransformType,
				FixedImageType,
				MovingImageType >  TransformInitializerType;

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

	// SaveImage<movingImageType>(ApplyTransform<MovingImageType, TransformType>(movingImageReader->GetOutput(), transform));
	
	std::cout<< "\n\n\n\n Affine Transform using itkMplus	\n\tThread: "<< metric->GetNumberOfThreads() <<
			"\nVariables: " <<transform->GetNumberOfParameters() <<

			std::endl;

			TransformType::ParametersType init_ = transform->GetParameters();
			std::cout << "Initial transform parameters: " << init_ << std::endl;

			init_ = registration->GetInitialTransformParameters();
			std::cout << "Initial transform parameters rec " << init_ << std::endl;

	RegularStepGradientDescentOptimizerCommandIterationUpdate::Pointer observer = RegularStepGradientDescentOptimizerCommandIterationUpdate::New();
	optimizer->AddObserver( itk::IterationEvent(), observer );
	if (fixedLabelMap && movingLabelMap && LABELREPORT > 0) {
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
      // Add a time probe
    itk::TimeProbesCollectorBase   chronometer;
    itk::MemoryProbesCollectorBase memorymeter;

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
			MovingImageType,
			FixedImageType >    ResampleFilterType;

	ResampleFilterType::Pointer resample = ResampleFilterType::New();

	resample->SetTransform( transform );
	resample->SetInput( movingImageReader->GetOutput() );

	resample->SetSize(    originalFixedImage->GetLargestPossibleRegion().GetSize() );
	resample->SetOutputOrigin(  originalFixedImage->GetOrigin() );
	resample->SetOutputSpacing( originalFixedImage->GetSpacing() );
	resample->SetOutputDirection( originalFixedImage->GetDirection() );
	resample->SetDefaultPixelValue( DFLTPIXELVALUE );

	typedef itk::ImageFileWriter< MovingImageType >  WriterType;
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
