// DeformableRegistration6
#include "itkImageRegistrationMethod.h"
#include "../../../Metrics/Mplus/itkMplus.h"
#include "itkBSplineTransform.h"
#include "itkLBFGSBOptimizer.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkBSplineResampleImageFunction.h"
#include "itkIdentityTransform.h"
#include "itkBSplineDecompositionImageFilter.h"
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
    ("gridposition,G", po::value<std::string>()->default_value("N"), "Read the position of the grid from a file")
    ("dfltpixelvalue", po::value<double>()->default_value(0), "Default pixel value")
    ("verbose,V", po::value<bool>()->default_value(false), "verbose")
    ("yota,y", po::value<double>(&YOTA)->default_value(0), "Yota value NC (Normalized Correlation) weight")
    ("yotaderivative,Y", po::value<double>(&YOTADERIVATIVE)->default_value(0), "Yota derivative NC, 0 = no derivatives")
    ("ngfpercentage", po::value<double>()->default_value(0.1), "NGF percentage of pixels used (0.1 = 10%)")
    ("msepercentage", po::value<double>()->default_value(0.1), "MSE percentage of pixels used (0.1 = 10%)")
    ("ncpercentage", po::value<double>()->default_value(0.1), "NC percentage of pixels used (0.1 = 10%)")
    ("rho", po::value<double>()->default_value(0.0), "Rho weight for Gradient Difference (GD)")
    ("rhoderivative", po::value<double>()->default_value(0.0), "Rho derivative for GD")
    ("sigma", po::value<double>()->default_value(0.0), "Sigma weight for Normalized Mutual Information (NMI)")
    ("sigmaderivative", po::value<double>()->default_value(0.0), "Sigma derivative for NMI")
    ("nmibins", po::value<int>()->default_value(64), "Number of histogram bins for NMI")
    ("derivativemode", po::value<int>()->default_value(0), "Derivative merge mode: 0=consistent, 2=main-metric adaptive (mode 1 not allowed for bsplines)")
    ("mainmetric", po::value<int>()->default_value(0), "Main metric index for mode 2: 0=MI, 1=NGF, 2=MSE, 3=NC, 4=Label, 5=GD, 6=NMI")
    ("ngfspacing", po::value<std::string>()->default_value("4,4,4"), "NGF spacing per dimension (x,y,z)")
    ("meshmarginsize", po::value<double>()->default_value(0.0), "Margin (mm) to extend mesh domain")
	("metricoverlap", po::value<bool>()->default_value(true), "Compute overlap between fixed and moving image (default true)")
	("fixedlabelmap",  po::value<std::string>()->default_value("N"), "Fixed label map filename (N = none)")
	("movinglabelmap", po::value<std::string>()->default_value("N"), "Moving label map filename (N = none)")
	("labelkappa",     po::value<double>()->default_value(0.0),       "Global kappa weight for label-map distance metric (0 = off)")
	("labelkappaderiv",po::value<double>()->default_value(0.0),       "Global kappa weight for label-map derivative")
	("labelkappavec",  po::value<std::string>()->default_value(""),   "Per-label kappa (value) weights: 'L1:w1,L2:w2,...'")
	("labelkappaderivvec", po::value<std::string>()->default_value(""),"Per-label kappa (derivative) weights: 'L1:w1,L2:w2,...'")
	("labelsamples",   po::value<unsigned int>()->default_value(20000),"Samples for label metric")
	("labelreport",    po::value<int>()->default_value(1),            "Report Dice every N iterations (0 = off)")
	("snapshotdir",    po::value<std::string>()->default_value("N"), "Directory for iteration snapshots (N = off)")
	("snapshotevery",  po::value<int>()->default_value(1),            "Save snapshot every N iterations")
	("snapshotstack",  po::value<bool>()->default_value(false),       "Save full 3D .nii.gz instead of mid-slice PNG")
	("version", "Print version and exit")
	("overlappadding", po::value<unsigned int>()->default_value(20), "Overlap padding in voxels")
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
	const unsigned int LABELSAMPLES   = vm["labelsamples"].as<unsigned int>();
	const int         LABELREPORT     = vm["labelreport"].as<int>();
	const std::string SNAPSHOTDIR     = vm["snapshotdir"].as<std::string>();
	const int         SNAPSHOTEVERY   = vm["snapshotevery"].as<int>();
	const bool        SNAPSHOTSTACK   = vm["snapshotstack"].as<bool>();

	const auto LABELKAPPAVEC      = RegCommon::ParseLabelWeights(vm["labelkappavec"].as<std::string>());
	const auto LABELKAPPADERIVVEC = RegCommon::ParseLabelWeights(vm["labelkappaderivvec"].as<std::string>());

	// Read label maps
	typedef itk::Image<short, ImageDimension> LabelImageType;
	LabelImageType::ConstPointer fixedLabelMap, movingLabelMap;
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

	ImageType::ConstPointer fixedImage = fixedImageReader->GetOutput();
	ImageType::ConstPointer movingImage = movingImageReader->GetOutput();

	// ── Input validation ──────────────────────────────────────────────────────
	if (!fixedImage || fixedImage->GetLargestPossibleRegion().GetNumberOfPixels() == 0) {
		std::cerr << "Error: Failed to load fixed image or image is empty." << std::endl;
		return EXIT_FAILURE;
	}
	if (!movingImage || movingImage->GetLargestPossibleRegion().GetNumberOfPixels() == 0) {
		std::cerr << "Error: Failed to load moving image or image is empty." << std::endl;
		return EXIT_FAILURE;
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
		// get the information of the image that specifies the position of the grid and overwrite the information of the fixed image
		meshspacing = meshImage->GetSpacing();
		meshorigin = meshImage->GetOrigin();
		meshdirection = meshImage->GetDirection();
		meshsize = meshImage->GetLargestPossibleRegion().GetSize();

		// resample the meshimage on the fixedimage space in casse the two images have different size
		typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
		ResampleFilterType::Pointer resampler = ResampleFilterType::New();
		resampler->SetInput(meshImage);
		resampler->SetOutputParametersFromImage(fixedImage);
		resampler->Update();
		ImageType::RegionType meshregionresampled = resampler->GetOutput()->GetLargestPossibleRegion();

		fixedImage->SetRequestedRegion(meshregionresampled);
	}

	for (unsigned int i = 0; i < SpaceDimension; ++i)
	{
				// correlate the number of extra nodes to the spline order:
				constexpr unsigned int borderNodesPerSide = SplineOrder;
			const double extension = borderNodesPerSide * GRIDRESOLUTION;

			// 1) shift the origin back by “extension”
			fixedOrigin[i] = meshorigin[i] - meshMargin - extension;

			// 2) grow the physical size by 2*extension
			fixedPhysicalDimensions[i] =
				meshspacing[i] * (meshsize[i] - 1)
				+ 2.0 * meshMargin
				+ 2.0 * extension;

			// 3) now recompute how many grid‐nodes you need
			const unsigned int totalGridNodes =
				static_cast<unsigned int>(fixedPhysicalDimensions[i] / GRIDRESOLUTION) + 1;

			// 4) subtract the spline order to get the final meshSize
			meshSize[i] = totalGridNodes > SplineOrder
						? totalGridNodes - SplineOrder
						: 1;  // guard against too small

			if (vm["verbose"].as<bool>())
			{
				std::cout
				<< "Dim " << i
				<< " origin = " << fixedOrigin[i]
				<< ", physSize = " << fixedPhysicalDimensions[i]
				<< ", meshSize = " << meshSize[i]
				<< std::endl;
			}
	}

	transform->SetTransformDomainOrigin(fixedOrigin);
	transform->SetTransformDomainPhysicalDimensions(
		fixedPhysicalDimensions);
	transform->SetTransformDomainMeshSize(meshSize);
	transform->SetTransformDomainDirection(fixedImage->GetDirection());
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
	metric->SetOverlapPadding(vm["overlappadding"].as<unsigned int>());
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

	if (TIN != "N")
	{

		transform = ReadTransform<TransformType>(TIN);
		registration->SetInitialTransformParameters(transform->GetParameters());
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
	if (SNAPSHOTDIR != "N") {
		auto snapObs = IterationSnapshotObserver<TransformType, ImageType>::New();
		snapObs->SetFixedImage(fixedImage);
		snapObs->SetMovingImage(movingImage);
		snapObs->SetTransform(transform);
		snapObs->SetOutputDirectory(SNAPSHOTDIR);
		snapObs->SetSaveEveryNIterations(static_cast<unsigned int>(SNAPSHOTEVERY));
		snapObs->SetSaveStack(SNAPSHOTSTACK);
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

	// Report the time and memory taken by the registration
	chronometer.Report(std::cout);
	memorymeter.Report(std::cout);

	transform->SetParameters(registration->GetLastTransformParameters());

	typedef itk::ResampleImageFilter<
		ImageType,
		ImageType>
		ResampleFilterType;

	ResampleFilterType::Pointer resample = ResampleFilterType::New();

	resample->SetTransform(transform);
	resample->SetInput(movingImageReader->GetOutput());

	resample->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
	resample->SetOutputOrigin(fixedImage->GetOrigin());
	resample->SetOutputSpacing(fixedImage->GetSpacing());
	resample->SetOutputDirection(fixedImage->GetDirection());
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
		td = TransformToDeformationField<TransformType, ImageType, DeformationTransformImageType>(transform, movingImageReader->GetOutput());
		WriteDeformationField<DeformationTransformImageType>(VOUT, td);
	};

	if (TOUT != "N")
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
