// DeformableRegistration6
#include "itkImageRegistrationMethod.h"
#include "../../../Metrics/Mplus/itkMplus.h"
#include "itkBSplineTransform.h"
#include "itkLBFGSBOptimizer.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkCastImageFilter.h"
#include "itkBSplineResampleImageFunction.h"
#include "itkIdentityTransform.h"
#include "itkBSplineDecompositionImageFilter.h"
#include "itkImageMaskSpatialObject.h"

#include "itkTransformToDeformationFieldSource.h"
#include "itkTransformFileWriter.h"
#include "itkTransformFileReader.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkCompositeTransform.h"
#include "../../Version.h"

#include "../../../Metrics/NGF/NGFImageMetric/NGFImageToImageMetric/Code/itkGetImageNoiseFunction.h"
#include "../../../includes/imageUtils.h"
#include "../../../includes/registrationUtils.h"
#include "../../../includes/RegistrationCommon.h"

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"

#include <boost/program_options.hpp>
#include <vector> // add this include at top
namespace po = boost::program_options;

const unsigned int ImageDimension = 3;

typedef float PixelType;
typedef itk::Image<PixelType, ImageDimension> ImageType;

const unsigned int SpaceDimension = ImageDimension;
const unsigned int SplineOrder = 3;
typedef double CoordinateRepType;

typedef itk::BSplineTransform<
	CoordinateRepType,
	SpaceDimension,
	SplineOrder>
	TransformType;

typedef itk::LBFGSBOptimizer OptimizerType;

typedef itk::Mplus<
	ImageType,
	ImageType>
	MetricType;

typedef itk::LinearInterpolateImageFunction<
	ImageType,
	double>
	InterpolatorType;

typedef itk::ImageRegistrationMethod<
	ImageType,
	ImageType>
	RegistrationType;

typedef itk::Vector<double, ImageDimension> VectorType;
typedef itk::Image<VectorType, ImageDimension> DeformationTransformImageType;

int main(int argc, char *argv[])
{
	po::options_description desc("B-spline Registration\n"
								 "Dr. Eros Montin Ph.D., 2014\n"
								 "eros.montin@gmail.com\n\n"
								 "cite us:\n\nMontin, E., Belfatto, A., Bologna, M., Meroni, S., Cavatorta, C., Pecori, E., Diletto, B., Massimino, M., Oprandi, M. C., Poggi, G., Arrigoni, F., Peruzzo, D., Pignoli, E., Gandola, L., Cerveri, P., & Mainardi, L. (2020). A multi-metric registration strategy for the alignment of longitudinal brain images in pediatric oncology. Medical & biological engineering & computing, 58(4), 843–855. https://doi.org/10.1007/s11517-019-02109-4\n\n"
								 "Allowed options for alpha MI + lambda NGF +  nu MSE + yota NMI");
	double YOTA = 0.0;
	double YOTADERIVATIVE = 0.0;
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
    ("bsplinecaching,B", po::value<bool>()->default_value(true), "B-spline caching, 1 for true")
    ("explicitPDFderivatives", po::value<bool>()->default_value(false), "Explicit PDF derivatives, 0 for false")
    ("lambda,l", po::value<double>()->default_value(0), "lambda value NGF 1.0")
    ("lambdaderivative,L", po::value<double>()->default_value(0), "Lambda derivative NGF 0 no derivatives")
    ("etavaluefixed,r", po::value<double>()->default_value(-1), "Eta value fixed image(NGF noise) -1 (autodetermine)")
    ("etavaluemoving,s", po::value<double>()->default_value(-1), "Eta value moving image (NGF noise) -1 (autodetermine)")
    ("NGFevaluator", po::value<int>()->default_value(0), "NGF Evaluator (0 scalar,1cross,2scdelta,3Delta,4Delta2)")
    ("ngfprecompute", po::value<bool>()->default_value(false), "Precompute moving-image NGF once and resample vector field each iteration (faster, approximate)")
    ("nu,n", po::value<double>()->default_value(0), "nu value MSE 1.0")
    ("nuderivative,N", po::value<double>()->default_value(0), "nu MSE derivative 1.0")
    ("gridresolution,g", po::value<double>()->default_value(50), "Mesh resolution (mm)")
    ("maxnumberofiterations,I", po::value<int>()->default_value(1000), "Max number of Iterations 1000")
    ("costfunctionconvergencefactor,F", po::value<double>()->default_value(1.e12), 
       "CostFunctionConvergenceFactor 1e+12 for low accuracy; 1e+7 for moderate accuracy and 1e+1 for extremely high accuracy.")
    ("projectedgradienttolerance,P", po::value<double>()->default_value(1.e-5), 
       "ProjectedGradientTolerance. Algorithm terminates when the project gradient is below the tolerance. Default value is 1e-5.")
    ("numberofevaluations,E", po::value<int>()->default_value(500), "Number of Evaluations")
    ("numberofcorrections,C", po::value<int>()->default_value(5), "Number of Corrections")
    ("fixedimagethreshold,t", po::value<double>()->default_value(-99999999), "Fixed image threshold")
    ("bound", po::value<int>()->default_value(0), 
       "Set the boundary condition for each variable, where = 0 if x[i] is unbounded, 1 if x[i] has only a lower bound,"
       " 2 if x[i] has both lower and upper bounds and 3 if x[i] has only an upper bound")
    ("lbound", po::value<double>()->default_value(0), "Lower bound")
    ("ubound", po::value<double>()->default_value(0), "Upper bound")
    ("transformout,T", po::value<std::string>()->default_value("N"), "Output for transform")
    ("transformin,W", po::value<std::string>()->default_value("N"), "Input no rigid transform for transform")
    ("gridposition,G", po::value<std::string>()->default_value("N"), "Read the position of the grid from a file; if it contains non-zero voxels, use their bounding box")
    ("dfltpixelvalue", po::value<double>()->default_value(0), "Default pixel value")
    ("verbose,V", po::value<bool>()->default_value(false), "verbose")
    ("yota,y", po::value<double>(&YOTA)->default_value(0), "Yota value NC (Normalized Correlation) weight")
    ("yotaderivative,Y", po::value<double>(&YOTADERIVATIVE)->default_value(0), "Yota derivative NC, 0 = no derivatives")
	("ngfpercentage", po::value<double>()->default_value(0.1), "NGF percentage of pixels used (0.1 = 10%)")
	("msepercentage", po::value<double>()->default_value(0.1), "MSE percentage of pixels used (0.1 = 10%)")
	("gdpercentage",   po::value<double>()->default_value(0.1), "GD percentage of pixels used (0.1 = 10%)")
	("nmipercentage",  po::value<double>()->default_value(0.1), "NMI percentage of pixels used (0.1 = 10%)")
	("ncpercentage", po::value<double>()->default_value(0.1), "NC percentage of pixels used (0.1 = 10%)")
    ("rho", po::value<double>()->default_value(0.0), "Rho weight for Gradient Difference (GD)")
    ("rhoderivative", po::value<double>()->default_value(0.0), "Rho derivative for GD")
    ("sigma", po::value<double>()->default_value(0.0), "Sigma weight for Normalized Mutual Information (NMI)")
    ("sigmaderivative", po::value<double>()->default_value(0.0), "Sigma derivative for NMI")
    ("nmibins", po::value<int>()->default_value(64), "Number of histogram bins for NMI")
    ("derivativemode", po::value<int>()->default_value(0), "Derivative merge mode: 0=consistent, 2=main-metric adaptive (mode 1 not allowed for bsplines)")
    ("mainmetric", po::value<int>()->default_value(0), "Main metric index for mode 2: 0=MI, 1=NGF, 2=MSE, 3=NC, 4=Label, 5=GD, 6=NMI")
    ("ngfspacing", po::value<std::string>()->default_value("4,4,4"), "NGF spacing per dimension (x,y,z)")
    ("workingresolution", po::value<std::string>()->default_value("0,0,0"),
       "Internal registration spacing in mm (x,y,z). Use 0,0,0 to keep the input spacing.")
    ("meshmarginsize", po::value<double>()->default_value(0.0), "Margin (mm) to extend mesh domain")
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
	("snapshotgrid",   po::value<bool>()->default_value(false),       "Overlay regular warped pixel-grid on snapshot panels (default off; when off, the real B-spline knot mesh is shown instead)")
	("snapshotgridspacing", po::value<unsigned int>()->default_value(20), "Grid line spacing in voxels for the deformation grid panel")
	("snapshotlinewidth", po::value<double>()->default_value(0.0),    "Overlay line width in pixels (0 = auto)")
	("version", "Print version and exit")
	("overlappadding", po::value<unsigned int>()->default_value(1),
		"Number of B-spline control points outside the image domain per side "
		"(min = spline order = 3). Higher values give more deformation support at image borders.")
	("metricpadding", po::value<unsigned int>()->default_value(0), "Metric overlap padding in voxels (default 20)")
	("modality", po::value<std::string>()->default_value("custom"),
		"Preset modality: 'multimodal' (MI+NGF), 'singlemodal' (MSE+NC), or 'custom' (manual weights)")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	// ── Handle --help and --version BEFORE any file I/O ────────────────────────
	if (vm.count("version")) {
		std::cout << "3DRegBsplines v" << VERSION_MAJOR << "." << VERSION_MINOR << std::endl;
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

	MetricType::Pointer metric = MetricType::New();
	OptimizerType::Pointer optimizer = OptimizerType::New();
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	RegistrationType::Pointer registration = RegistrationType::New();

	auto s = vm["ngfspacing"].as<std::string>();
	std::replace(s.begin(), s.end(), ',', ' ');
	std::istringstream iss(s);
	std::vector<double> tmp((std::istream_iterator<double>(iss)),
							std::istream_iterator<double>());
	if (tmp.size() != ImageDimension)
	{
		std::cerr << "Error: --ngfspacing must have exactly " << ImageDimension
		          << " comma-separated values, got " << tmp.size() << std::endl;
		return EXIT_FAILURE;
	}
	ImageType::SpacingType ngf;
	for (unsigned i = 0; i < ImageDimension; ++i)
		ngf[i] = tmp[i];

	ImageType::SpacingType workingSpacing;
	bool useWorkingResolution = false;
	if (!RegCommon::ParseOptionalSpacing<ImageType::SpacingType, ImageDimension>(
	        vm["workingresolution"].as<std::string>(),
	        "workingresolution",
	        workingSpacing,
	        useWorkingResolution))
	{
		return EXIT_FAILURE;
	}

	if (vm["verbose"].as<bool>())
		RegCommon::PrintOptions(vm);

	std::string fixedImageFN = vm["fixedimage"].as<std::string>();
	std::string movingImageFN = vm["movingimage"].as<std::string>();

	std::string ou = vm["outputimage"].as<std::string>();
	std::string VOUT = vm["vfout"].as<std::string>();

	int NT = vm["numberofthreads"].as<int>();
	int NB = vm["mattesnumberofbins"].as<int>();
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
	if (NGFevaluator < 0 || NGFevaluator > 4)
	{
		std::cerr << "Error: NGFevaluator must be between 0 and 4" << std::endl;
		return EXIT_FAILURE;
	}
	double GRIDRESOLUTION = vm["gridresolution"].as<double>();
	int NI = vm["maxnumberofiterations"].as<int>();
	double CFCF = vm["costfunctionconvergencefactor"].as<double>();
	double PGT = vm["projectedgradienttolerance"].as<double>();
	int NE = vm["numberofevaluations"].as<int>();
	int NC = vm["numberofcorrections"].as<int>();
	double TR = vm["fixedimagethreshold"].as<double>();
	bool TB = vm["bsplinecaching"].as<bool>();
	bool EPDF = vm["explicitPDFderivatives"].as<bool>();
	int BOUND = vm["bound"].as<int>();
	double LBOUND = vm["lbound"].as<double>();
	double UBOUND = vm["ubound"].as<double>();
	std::string TOUT = vm["transformout"].as<std::string>();
	std::string TIN = vm["transformin"].as<std::string>();
	double DFLTPIXELVALUE = vm["dfltpixelvalue"].as<double>();
	bool V = vm["verbose"].as<bool>();
	std::string GRIDPOSITION = vm["gridposition"].as<std::string>();
	// double RHO = vm["rho"].as<double>();
	// double RHODERIVATIVE = vm["rhoderivative"].as<double>();
	double RHO             = vm["rho"].as<double>();
	double RHODERIVATIVE   = vm["rhoderivative"].as<double>();
	double SIGMA           = vm["sigma"].as<double>();
	double SIGMADERIVATIVE = vm["sigmaderivative"].as<double>();
	int    NMIBINS         = vm["nmibins"].as<int>();
	bool METRICOVERLAP = vm["metricoverlap"].as<bool>();
	int  DERIVMODE = vm["derivativemode"].as<int>();
	int  MAINMETRIC = vm["mainmetric"].as<int>();
	if (DERIVMODE == 1) {
		std::cerr << "ERROR: derivativemode=1 (normalized) is not compatible with "
		             "LBFGS-B optimizer used by B-splines. Use 0 or 2." << std::endl;
		return EXIT_FAILURE;
	}
	double meshMargin = vm["meshmarginsize"].as<double>();

	// ── label map options ───────────────────────────────────────────────────────
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

	// Read label maps
	typedef itk::Image<short, ImageDimension> LabelImageType;
	LabelImageType::ConstPointer fixedLabelMap;
	LabelImageType::ConstPointer movingLabelMap;
	if (FIXEDLABELMAP != "N" && MOVINGLABELMAP != "N")
	{
		typedef itk::ImageFileReader<LabelImageType> LabelReaderType;
		auto flr = LabelReaderType::New(); flr->SetFileName(FIXEDLABELMAP); flr->Update();
		{ LabelImageType::Pointer _tmp = flr->GetOutput(); _tmp->DisconnectPipeline(); fixedLabelMap = _tmp; }
		auto mlr = LabelReaderType::New(); mlr->SetFileName(MOVINGLABELMAP); mlr->Update();
		{ LabelImageType::Pointer _tmp = mlr->GetOutput(); _tmp->DisconnectPipeline(); movingLabelMap = _tmp; }
		std::cout << "[Label] Fixed label map:  " << FIXEDLABELMAP  << "\n"
		          << "[Label] Moving label map: " << MOVINGLABELMAP << std::endl;
	}

	registration->SetMetric(metric);
	registration->SetOptimizer(optimizer);
	registration->SetInterpolator(interpolator);
	registration->SetNumberOfThreads(NT);
	registration->SetFixedImageRegion(ImageType::RegionType());

	TransformType::Pointer transform = TransformType::New();
	registration->SetTransform(transform);

	typedef itk::ImageFileReader<ImageType> FixedImageReaderType;
	typedef itk::ImageFileReader<ImageType> MovingImageReaderType;

	FixedImageReaderType::Pointer fixedImageReader = FixedImageReaderType::New();
	MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

	fixedImageReader->SetFileName(fixedImageFN);
	movingImageReader->SetFileName(movingImageFN);

	fixedImageReader->Update();
	movingImageReader->Update();

	ImageType::ConstPointer originalFixedImage = fixedImageReader->GetOutput();
	ImageType::ConstPointer originalMovingImage = movingImageReader->GetOutput();
	ImageType::ConstPointer fixedImage = originalFixedImage;
	ImageType::Pointer movingImage = const_cast<ImageType*>(originalMovingImage.GetPointer());

	// ── Input validation ──────────────────────────────────────────────────────
	if (!originalFixedImage || originalFixedImage->GetLargestPossibleRegion().GetNumberOfPixels() == 0) {
		std::cerr << "Error: Failed to load fixed image or image is empty." << std::endl;
		return EXIT_FAILURE;
	}
	if (!originalMovingImage || originalMovingImage->GetLargestPossibleRegion().GetNumberOfPixels() == 0) {
		std::cerr << "Error: Failed to load moving image or image is empty." << std::endl;
		return EXIT_FAILURE;
	}

	if (useWorkingResolution)
	{
		std::cout << "[WorkingResolution] Resampling registration inputs to spacing "
		          << workingSpacing << std::endl;
		if (!RegCommon::SpacingEquals(originalFixedImage->GetSpacing(), workingSpacing))
		{
			fixedImage = RegCommon::ResampleScalarImageToSpacing<ImageType>(
			    originalFixedImage, workingSpacing);
		}
		if (!RegCommon::SpacingEquals(originalMovingImage->GetSpacing(), workingSpacing))
		{
			movingImage = RegCommon::ResampleScalarImageToSpacing<ImageType>(
			    originalMovingImage, workingSpacing, DFLTPIXELVALUE);
		}
		if (fixedLabelMap && !RegCommon::SpacingEquals(fixedLabelMap->GetSpacing(), workingSpacing))
		{
			fixedLabelMap = RegCommon::ResampleNearestNeighborImageToSpacing<LabelImageType>(
			    fixedLabelMap, workingSpacing);
		}
		if (movingLabelMap && !RegCommon::SpacingEquals(movingLabelMap->GetSpacing(), workingSpacing))
		{
			movingLabelMap = RegCommon::ResampleNearestNeighborImageToSpacing<LabelImageType>(
			    movingLabelMap, workingSpacing);
		}
	}

	registration->SetFixedImage(fixedImage);
	registration->SetMovingImage(movingImage);

	ImageType::RegionType fixedRegion = fixedImage->GetBufferedRegion();
	registration->SetFixedImageRegion(fixedRegion);

	TransformType::PhysicalDimensionsType fixedPhysicalDimensions;
	TransformType::MeshSizeType meshSize;
	TransformType::OriginType fixedOrigin;

	ImageType::SpacingType meshspacing = fixedImage->GetSpacing();
	ImageType::PointType meshorigin = fixedImage->GetOrigin();
	ImageType::DirectionType meshdirection = fixedImage->GetDirection();
	ImageType::SizeType meshsize = fixedImage->GetLargestPossibleRegion().GetSize();

	// #let's fix a few things
	if ((LAMBDA != 0) || (LAMBDADERIVATIVE != 0))
	{

		if ((ETAF == -1) || (ETAM == -1))
		{
			metric->SetAutoEstimateEta(true);
		}
	}

	if (GRIDPOSITION != "N")
	{
		// read the image that specifies the position of the grid
		FixedImageReaderType::Pointer meshImageReader = FixedImageReaderType::New();
		meshImageReader->SetFileName(GRIDPOSITION);
		meshImageReader->Update();
		ImageType::ConstPointer meshImage = meshImageReader->GetOutput();
		if (useWorkingResolution &&
		    !RegCommon::SpacingEquals(meshImage->GetSpacing(), workingSpacing))
		{
			meshImage = RegCommon::ResampleNearestNeighborImageToSpacing<ImageType>(
			    meshImage, workingSpacing);
		}

		// If the input is a mask/object image, shrink the mesh domain to the
		// non-zero voxel bounding box. If it is empty, preserve legacy behavior
		// and use the full image geometry.
		ImageType::IndexType bboxMinIndex;
		ImageType::IndexType bboxMaxIndex;
		bool hasActiveMaskVoxel = false;
		using MeshIteratorType = itk::ImageRegionConstIteratorWithIndex<ImageType>;
		for (MeshIteratorType it(meshImage, meshImage->GetLargestPossibleRegion());
			 !it.IsAtEnd(); ++it)
		{
			if (it.Get() == PixelType{})
			{
				continue;
			}

			const ImageType::IndexType currentIndex = it.GetIndex();
			if (!hasActiveMaskVoxel)
			{
				bboxMinIndex = currentIndex;
				bboxMaxIndex = currentIndex;
				hasActiveMaskVoxel = true;
				continue;
			}

			for (unsigned int d = 0; d < ImageDimension; ++d)
			{
				bboxMinIndex[d] = std::min(bboxMinIndex[d], currentIndex[d]);
				bboxMaxIndex[d] = std::max(bboxMaxIndex[d], currentIndex[d]);
			}
		}

		// get the information of the image that specifies the position of the grid and overwrite the information of the fixed image
		meshspacing = meshImage->GetSpacing();
		meshdirection = meshImage->GetDirection();
		meshorigin = meshImage->GetOrigin();
		meshsize = meshImage->GetLargestPossibleRegion().GetSize();

		if (hasActiveMaskVoxel)
		{
			ImageType::PointType bboxOrigin;
			meshImage->TransformIndexToPhysicalPoint(bboxMinIndex, bboxOrigin);
			meshorigin = bboxOrigin;

			ImageType::SizeType bboxSize;
			for (unsigned int d = 0; d < ImageDimension; ++d)
			{
				bboxSize[d] = static_cast<ImageType::SizeType::SizeValueType>(
					bboxMaxIndex[d] - bboxMinIndex[d] + 1);
			}
			meshsize = bboxSize;

			if (vm["verbose"].as<bool>())
			{
				std::cout << "[GridPosition] Using non-zero mask bounding box";
				std::cout << " minIndex=[";
				for (unsigned int d = 0; d < ImageDimension; ++d)
				{
					if (d != 0) std::cout << ", ";
					std::cout << bboxMinIndex[d];
				}
				std::cout << "] maxIndex=[";
				for (unsigned int d = 0; d < ImageDimension; ++d)
				{
					if (d != 0) std::cout << ", ";
					std::cout << bboxMaxIndex[d];
				}
				std::cout << "] origin=[";
				for (unsigned int d = 0; d < ImageDimension; ++d)
				{
					if (d != 0) std::cout << ", ";
					std::cout << meshorigin[d];
				}
				std::cout << "] size=[";
				for (unsigned int d = 0; d < ImageDimension; ++d)
				{
					if (d != 0) std::cout << ", ";
					std::cout << meshsize[d];
				}
				std::cout << "]" << std::endl;
			}
		}
		else if (vm["verbose"].as<bool>())
		{
			std::cout << "[GridPosition] No non-zero voxels found in " << GRIDPOSITION
			          << "; falling back to full image geometry." << std::endl;
		}

		// resample the meshimage on the fixedimage space in casse the two images have different size
		typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
		ResampleFilterType::Pointer resampler = ResampleFilterType::New();
		resampler->SetInput(meshImage);
		resampler->SetOutputParametersFromImage(fixedImage);
		resampler->Update();
		ImageType::RegionType meshregionresampled = resampler->GetOutput()->GetLargestPossibleRegion();

		fixedImage->SetRequestedRegion(meshregionresampled);
	}

	itk::Vector<double, SpaceDimension> axisPadding;
	axisPadding.Fill(0.0);
	for (unsigned int i = 0; i < SpaceDimension; ++i)
	{
			// Number of extra B-spline control points outside the image domain
			// per side. ITK internally handles the SplineOrder border
			// coefficients; this setting only controls how far the domain
			// extends beyond the anatomy for better edge deformation.
			const unsigned int borderNodesPerSide =
				vm["overlappadding"].as<unsigned int>();
			const double extension = borderNodesPerSide * GRIDRESOLUTION;
			const double totalPad = meshMargin + extension;
			axisPadding[i] = totalPad;

			// Grow the physical size along each transform-domain axis.
			fixedPhysicalDimensions[i] =
				meshspacing[i] * (meshsize[i] - 1)
				+ 2.0 * totalPad;

			// Recompute how many grid-nodes you need.
			const unsigned int totalGridNodes =
				static_cast<unsigned int>(fixedPhysicalDimensions[i] / GRIDRESOLUTION) + 1;

			// Subtract the spline order to get the final meshSize.
			meshSize[i] = totalGridNodes > SplineOrder
						? totalGridNodes - SplineOrder
						: 1;  // guard against too small

			if (vm["verbose"].as<bool>())
			{
				std::cout
				<< "Dim " << i
				<< ", axisPad = " << totalPad
				<< ", physSize = " << fixedPhysicalDimensions[i]
				<< ", meshSize = " << meshSize[i]
				<< std::endl;
			}
	}

	// Shift the mesh origin "before" the anatomy along the full transform
	// domain basis, not just by the sign of the diagonal. This keeps the
	// domain placement correct for axis permutations and oblique directions.
	for (unsigned int row = 0; row < SpaceDimension; ++row)
	{
		fixedOrigin[row] = meshorigin[row];
		for (unsigned int col = 0; col < SpaceDimension; ++col)
		{
			fixedOrigin[row] -= meshdirection[row][col] * axisPadding[col];
		}
	}

	if (vm["verbose"].as<bool>())
	{
		std::cout << "[GridPosition] Transform domain origin = [";
		for (unsigned int d = 0; d < SpaceDimension; ++d)
		{
			if (d != 0) std::cout << ", ";
			std::cout << fixedOrigin[d];
		}
		std::cout << "] direction = [";
		for (unsigned int row = 0; row < SpaceDimension; ++row)
		{
			if (row != 0) std::cout << "; ";
			for (unsigned int col = 0; col < SpaceDimension; ++col)
			{
				if (col != 0) std::cout << ", ";
				std::cout << meshdirection[row][col];
			}
		}
		std::cout << "]" << std::endl;
	}

	transform->SetTransformDomainOrigin(fixedOrigin);
	transform->SetTransformDomainPhysicalDimensions(
		fixedPhysicalDimensions);
	transform->SetTransformDomainMeshSize(meshSize);
	transform->SetTransformDomainDirection(meshdirection);
	transform->SetIdentity();

	const unsigned int numberOfGridNodes = transform->GetNumberOfParameters() / SpaceDimension;

	// metric is negative
	optimizer->MinimizeOn();

	registration->SetInitialTransformParameters(transform->GetParameters());

	const unsigned int numberOfPixels = fixedImage->GetLargestPossibleRegion().GetNumberOfPixels();


	double NGFPERCENTAGE = vm["ngfpercentage"].as<double>();
	double MSEPERCENTAGE = vm["msepercentage"].as<double>();
	double NCPERCENTAGE = vm["ncpercentage"].as<double>();
	// double CHPERCENTAGE = vm["chpercentage"].as<double>();
	
	const unsigned int numberOfSamplesMA = static_cast<unsigned int>(numberOfPixels * MAPERCENTAGE);
	// const unsigned int numberOfSamplesCH =static_cast<unsigned int>(numberOfPixels * CHPERCENTAGE);
	const unsigned int numberOfSamplesNGF = static_cast<unsigned int>(numberOfPixels * NGFPERCENTAGE);
	const unsigned int numberOfSamplesMSE = static_cast<unsigned int>(numberOfPixels * MSEPERCENTAGE);
	const unsigned int numberOfSamplesNC = static_cast<unsigned int>(numberOfPixels * NCPERCENTAGE);

	metric->SetUseCachingOfBSplineWeights(TB);
	metric->SetUseExplicitPDFDerivatives(EPDF);
	metric->SetNumberOfThreads(NT);
	metric->SetDerivativeMode(DERIVMODE);
	metric->SetMainMetricIndex(MAINMETRIC);
	metric->SetComputeOverlap(METRICOVERLAP);
	metric->SetOverlapPadding(vm["metricpadding"].as<unsigned int>());
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


	if (TR != -99999999)
	{
		metric->SetFixedImageThreshold(TR);
	}

	// ── label metric wiring ─────────────────────────────────────────────────────
	if (fixedLabelMap && movingLabelMap)
	{
		metric->SetFixedLabelMap(fixedLabelMap);
		metric->SetMovingLabelMap(movingLabelMap);
		metric->SetLabelKappa(LABELKAPPA);
		metric->SetLabelKappaDerivative(LABELKAPPADERIV);
		metric->SetLabelNumberOfSamples(LABELSAMPLES);
		if (!LABELKAPPAVEC.empty())      metric->SetLabelKappaWeights(LABELKAPPAVEC);
		if (!LABELKAPPADERIVVEC.empty()) metric->SetLabelKappaDerivativeWeights(LABELKAPPADERIVVEC);
	}

	const unsigned int numParameters = transform->GetNumberOfParameters();
	OptimizerType::BoundSelectionType boundSelect(numParameters);
	OptimizerType::BoundValueType upperBound(numParameters);
	OptimizerType::BoundValueType lowerBound(numParameters);

	boundSelect.Fill(BOUND);
	std::cout << "Optimization Bound" << BOUND << std::endl;
	switch (BOUND)
	{
	case 1:
		lowerBound.Fill(LBOUND);
		std::cout << "Lower Bound" << LBOUND << std::endl;
		break;
	case 2:
		lowerBound.Fill(LBOUND);
		upperBound.Fill(UBOUND);
		std::cout << "Lower Bound" << LBOUND << " and " << "Upper Bound" << UBOUND << std::endl;
		break;
	case 3:
		upperBound.Fill(UBOUND);
		std::cout << "Upper Bound" << UBOUND << std::endl;
		break;
	default:
		std::cout << "Ubounded " << UBOUND << std::endl;
		break;
	}

	optimizer->SetBoundSelection(boundSelect);
	optimizer->SetUpperBound(upperBound);
	optimizer->SetLowerBound(lowerBound);
	// CostFunctionConvergenceFactor 1e+12 for low accuracy; 1e+7 for moderate accuracy and 1e+1 for extremely high accuracy.
	optimizer->SetCostFunctionConvergenceFactor(CFCF);

	if (V)
	{
		optimizer->TraceOn();
	}

	optimizer->SetProjectedGradientTolerance(PGT);
	optimizer->SetMaximumNumberOfIterations(NI);
	optimizer->SetMaximumNumberOfEvaluations(NE);
	optimizer->SetMaximumNumberOfCorrections(NC);

	// Optional initial linear transform (set when --transformin contains an
	// affine/rigid/similarity that is used to pre-warp the moving image).
	typedef itk::MatrixOffsetTransformBase<double, ImageDimension, ImageDimension> GenericLinearTransformType;
	GenericLinearTransformType::Pointer initialLinearTransform;

	if (TIN != "N")
	{
		// Try reading as a B-spline transform first
		TransformType::Pointer bsplineIn = ReadTransform<TransformType>(TIN);
		if (bsplineIn)
		{
			transform = bsplineIn;
			registration->SetInitialTransformParameters(transform->GetParameters());
			std::cout << "[TransformIn] Loaded B-spline transform from " << TIN << std::endl;
		}
		else
		{
			// Not a B-spline — read as generic transform and pre-warp the moving image
			itk::TransformBase::Pointer genericTransform = ReadTransformGeneric(TIN);
			if (!genericTransform)
			{
				std::cerr << "Error: could not read transform from " << TIN << std::endl;
				return EXIT_FAILURE;
			}

			initialLinearTransform =
			    dynamic_cast<GenericLinearTransformType*>(genericTransform.GetPointer());
			if (!initialLinearTransform)
			{
				std::cerr << "Error: transform in " << TIN
				          << " is neither a BSplineTransform nor a linear transform (Affine/Rigid/Similarity).\n"
				          << "  Actual type: " << genericTransform->GetTransformTypeAsString() << std::endl;
				return EXIT_FAILURE;
			}

			std::cout << "[TransformIn] Loaded linear transform (" << genericTransform->GetTransformTypeAsString()
			          << ") from " << TIN << "\n"
			          << "  Pre-warping moving image and label map before B-spline registration." << std::endl;

			// Pre-warp moving image
			typedef itk::ResampleImageFilter<ImageType, ImageType> PreWarpFilterType;
			PreWarpFilterType::Pointer prewarp = PreWarpFilterType::New();
			prewarp->SetInput(movingImageReader->GetOutput());
			prewarp->SetTransform(initialLinearTransform);
			prewarp->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
			prewarp->SetOutputOrigin(fixedImage->GetOrigin());
			prewarp->SetOutputSpacing(fixedImage->GetSpacing());
			prewarp->SetOutputDirection(fixedImage->GetDirection());
			prewarp->SetDefaultPixelValue(DFLTPIXELVALUE);
			prewarp->Update();

			ImageType::Pointer prewarpedMoving = prewarp->GetOutput();
			prewarpedMoving->DisconnectPipeline();
			movingImage = prewarpedMoving;
			registration->SetMovingImage(movingImage);

			// Pre-warp moving label map if provided
			if (movingLabelMap)
			{
				typedef itk::ResampleImageFilter<LabelImageType, LabelImageType> LabelPreWarpType;
				typedef itk::NearestNeighborInterpolateImageFunction<LabelImageType, double> NNInterpType;
				LabelPreWarpType::Pointer labelPrewarp = LabelPreWarpType::New();
				NNInterpType::Pointer nnInterp = NNInterpType::New();
				labelPrewarp->SetInput(movingLabelMap);
				labelPrewarp->SetTransform(initialLinearTransform);
				labelPrewarp->SetInterpolator(nnInterp);
				labelPrewarp->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
				labelPrewarp->SetOutputOrigin(fixedImage->GetOrigin());
				labelPrewarp->SetOutputSpacing(fixedImage->GetSpacing());
				labelPrewarp->SetOutputDirection(fixedImage->GetDirection());
				labelPrewarp->SetDefaultPixelValue(0);
				labelPrewarp->Update();
				LabelImageType::Pointer prewarpedLabel = labelPrewarp->GetOutput();
				prewarpedLabel->DisconnectPipeline();
				movingLabelMap = prewarpedLabel;
			}

			std::cout << "  Pre-warping complete." << std::endl;
		}
	}

	// if TODO
	// if (VIN!='N')
	// {
	// 	DeformationTransformImageType::Pointer td  = DeformationTransformImageType::New();
	// 	td=ReadDeformationField<DeformationTransformImageType>(VIN);
	// 	transform=DeformationFieldToTransform<TransformType,ImageType,DeformationTransformImageType>(td, fixedImage);
	// 	registration->SetInitialTransformParameters(transform->GetParameters());

	// }
	std::cout << "\n B-spline transform using itkMplus,	Threads: " << metric->GetNumberOfThreads() << ", Variables: " << transform->GetNumberOfParameters() << ", Grid Nodes " << numberOfGridNodes << "\nTransform Domain meshes: " << transform->GetTransformDomainMeshSize() << "\nMontin, E., et al. A multi-metric registration strategy for the alignment of longitudinal brain images in pediatric oncology. Med Biol Eng Comput 58, 843-855 (2020). https://doi.org/10.1007/s11517-019-02109-4" << std::endl;

	LBFGSBOptimizeCommandIterationUpdate::Pointer observer = LBFGSBOptimizeCommandIterationUpdate::New();
	optimizer->AddObserver(itk::IterationEvent(), observer);

	// ── label Dice monitoring ───────────────────────────────────────────────────
	if (fixedLabelMap && movingLabelMap && LABELREPORT > 0)
	{
		LabelMapDiceObserver<TransformType, LabelImageType>::Pointer labelObs =
		    LabelMapDiceObserver<TransformType, LabelImageType>::New();
		labelObs->SetFixedLabelMap(fixedLabelMap);
		labelObs->SetMovingLabelMap(movingLabelMap);
		labelObs->SetTransform(transform);
		labelObs->SetEvaluateEveryNIterations(static_cast<unsigned int>(LABELREPORT));
		optimizer->AddObserver(itk::IterationEvent(), labelObs);
	}
	using SnapObsType = IterationSnapshotObserver<TransformType, ImageType>;
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
		// When --snapshotgrid is off (the default for B-splines), show the
		// real B-spline control-point lattice instead of a pixel grid.
		snapObs->SetShowBSplineMesh(!SNAPSHOTGRID);
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
	itk::TimeProbesCollectorBase chronometer;
	itk::MemoryProbesCollectorBase memorymeter;

	try
	{
		memorymeter.Start("Registration");
		chronometer.Start("Registration");

		registration->Update();

		chronometer.Stop("Registration");
		memorymeter.Stop("Registration");
		if (vm["verbose"].as<bool>())
		{
			std::cout << "Registration completed successfully." << std::endl;
			std::cout << "Optimizer stop condition = "
					  << registration->GetOptimizer()->GetStopConditionDescription()
					  << std::endl;
		}
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}

	if (snapObs) snapObs->FinalizeConvergencePlot();

	// Report the time and memory taken by the registration
	chronometer.Report(std::cout);
	memorymeter.Report(std::cout);

	transform->SetParameters(registration->GetLastTransformParameters());

	// ── Build the final output transform ──────────────────────────────────────
	// If a linear pre-warp was used, compose: linear → B-spline
	typedef itk::CompositeTransform<double, ImageDimension> CompositeTransformType;
	typedef itk::Transform<double, ImageDimension, ImageDimension> GenericTransformType;
	typename CompositeTransformType::Pointer compositeTransform;
	const GenericTransformType* outputTransform = transform.GetPointer();

	if (initialLinearTransform)
	{
		compositeTransform = CompositeTransformType::New();
		compositeTransform->AddTransform(initialLinearTransform);
		compositeTransform->AddTransform(transform);
		outputTransform = compositeTransform.GetPointer();
		std::cout << "[Output] Composing initial linear + B-spline transform." << std::endl;
	}

	typedef itk::ResampleImageFilter<
		ImageType,
		ImageType>
		ResampleFilterType;

	ResampleFilterType::Pointer resample = ResampleFilterType::New();

	resample->SetTransform(outputTransform);
	resample->SetInput(movingImageReader->GetOutput());

	resample->SetSize(originalFixedImage->GetLargestPossibleRegion().GetSize());
	resample->SetOutputOrigin(originalFixedImage->GetOrigin());
	resample->SetOutputSpacing(originalFixedImage->GetSpacing());
	resample->SetOutputDirection(originalFixedImage->GetDirection());
	resample->SetDefaultPixelValue(DFLTPIXELVALUE);

	typedef itk::ImageFileWriter<ImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();

	writer->SetFileName(ou);
	writer->SetInput(resample->GetOutput());

	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}

	if (VOUT != "N")
	{
		DeformationTransformImageType::Pointer td = DeformationTransformImageType::New();
		if (initialLinearTransform)
			td = TransformToDeformationField<CompositeTransformType, ImageType, DeformationTransformImageType>(compositeTransform, movingImageReader->GetOutput());
		else
			td = TransformToDeformationField<TransformType, ImageType, DeformationTransformImageType>(transform, movingImageReader->GetOutput());
		WriteDeformationField<DeformationTransformImageType>(VOUT, td);
	};

	if (TOUT != "N")
	{
#if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 5) || ITK_VERSION_MAJOR > 4
		itk::TransformFileWriterTemplate<double>::Pointer twriter =
			itk::TransformFileWriterTemplate<double>::New();
#else
		itk::TransformFileWriter::Pointer twriter = itk::TransformFileWriter::New();
#endif
		if (initialLinearTransform)
		{
			// Write both transforms so the composite can be reconstructed
			twriter->SetInput(initialLinearTransform);
			twriter->AddTransform(registration->GetOutput()->Get());
		}
		else
		{
			twriter->SetInput(registration->GetOutput()->Get());
		}
		twriter->SetFileName(TOUT);
		twriter->Update();
	};

	return EXIT_SUCCESS;
}
