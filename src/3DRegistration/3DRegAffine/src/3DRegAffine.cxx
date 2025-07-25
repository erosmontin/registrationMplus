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

#include "../../../Metrics/NGF/NGFImageMetric/NGFImageToImageMetric/Code/itkGetImageNoiseFunction.h"
#include "../../../includes/imageUtils.h"
#include "../../../includes/registrationUtils.h"

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"

#include "itkCenteredTransformInitializer.h"
#include "itkRegularStepGradientDescentOptimizer.h"

#include <boost/program_options.hpp>
#include <algorithm>    // std::replace
#include <iterator>     // std::istream_iterator
#include <sstream>      // std::istringstream
namespace po = boost::program_options;

const unsigned int ImageDimension = 3;

	typedef  float          PixelType;
	typedef itk::Image< PixelType, ImageDimension >  FixedImageType;
	typedef itk::Image< PixelType, ImageDimension >  MovingImageType;

	const unsigned int SpaceDimension = ImageDimension;
	const unsigned int SplineOrder = 3;
	typedef double CoordinateRepType;
	
	typedef itk::AffineTransform<double, 3> TransformType;


	typedef itk::Mplus<
			FixedImageType,
			MovingImageType >    MetricType;

	typedef itk:: LinearInterpolateImageFunction<
			MovingImageType,
			double          >    InterpolatorType;

	typedef itk::ImageRegistrationMethod<
			FixedImageType,
			MovingImageType >    RegistrationType;

	MetricType::Pointer         metric        = MetricType::New();
	InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
	RegistrationType::Pointer   registration  = RegistrationType::New();

typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
int main( int argc, char *argv[] )
{
    po::options_description desc("B-spline Registration\n"
	"Dr. Eros Montin Ph.D., 2014\n"
	"eros.montin@gmail.com\n\n"
	"cite us:\n\nMontin, E., Belfatto, A., Bologna, M., Meroni, S., Cavatorta, C., Pecori, E., Diletto, B., Massimino, M., Oprandi, M. C., Poggi, G., Arrigoni, F., Peruzzo, D., Pignoli, E., Gandola, L., Cerveri, P., & Mainardi, L. (2020). A multi-metric registration strategy for the alignment of longitudinal brain images in pediatric oncology. Medical & biological engineering & computing, 58(4), 843–855. https://doi.org/10.1007/s11517-019-02109-4\n\n"
	"Allowed options for alpha MI + lambda NGF +  nu MSE +yota NMI\n\n");
    std::string method;
	bool ND=false;
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
		("normalizederivatives,Z", po::value<bool>(&ND)->default_value(false), "Normalize derivatives")
        ("yota,y", po::value<double>(&YOTA)->default_value(0.1), "Yota value NMI 0.1")
        ("yotaderivative,Y", po::value<double>(&YOTADERIVATIVE)->default_value(0), "Yota derivative NMI 0 no derivatives")
        ("msepercentage",   po::value<double>()->default_value(0.1), "MSE percentage of pixels used (0.1 = 10%)")
        ("ngfpercentage",   po::value<double>()->default_value(0.1), "NGF percentage of pixels used (0.1 = 10%)")
        ("ncpercentage,p",  po::value<double>()->default_value(0.1), "Histogram‐Correlation (CH) percentage of pixels used (0.1 = 10%)")
        // ("rho",             po::value<double>()->default_value(0.0), "rho weight for CH")
        // ("rhoderivative",   po::value<double>()->default_value(0.0), "rho derivative for CH")
        ("ngfspacing",      po::value<std::string>()->default_value("4,4,4"), "NGF spacing per dimension (x,y,z)")
 ;
	


    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

	bool NORMALIZE_DERIVATIVES = vm["normalizederivatives"].as<bool>();

	if (vm.count("help") || !vm.count("fixedimage") || !vm.count("movingimage") || !vm.count("outputimage") || !vm.count("vfout")){
    std::cout << desc << "\n";
    return 1;
	}


	auto s = vm["ngfspacing"].as<std::string>();
	std::replace(s.begin(), s.end(), ',', ' ');
	std::istringstream iss(s);
	std::vector<double> tmp((std::istream_iterator<double>(iss)),
							std::istream_iterator<double>());
	if (tmp.size() != ImageDimension)
	{ /* error */
	}
	FixedImageType::SpacingType ngf;
	for (unsigned i = 0; i < ImageDimension; ++i)
		ngf[i] = tmp[i];
		if (vm["verbose"].as<bool>())
		{
for(const auto& it : vm) {
	std::cout << "Option: " << it.first.c_str() << "\nValue: ";
	
	auto& value = it.second.value();
	
	if(value.type() == typeid(int)) {
		auto v = boost::any_cast<int>(&value);
		std::cout << *v;
	} else if(value.type() == typeid(double)) {
		auto v = boost::any_cast<double>(&value);
		std::cout << *v;
	} else if(value.type() == typeid(std::string)) {
		auto v = boost::any_cast<std::string>(&value);
		std::cout << *v;
	} else if(value.type() == typeid(bool)) {
		auto v = boost::any_cast<bool>(&value);
		std::cout << std::boolalpha << *v;
	} else {
		std::cout << "Unknown type";
	}
	std::cout << "\n-------------------\n";
}
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

    double MSEPERCENTAGE   = vm["msepercentage"].as<double>();
    double NGFPERCENTAGE   = vm["ngfpercentage"].as<double>();
    double NMIPERCENTAGE   = vm["nmipercentage"].as<double>();
    double NCPERCENTAGE   = vm["ncpercentage"].as<double>();
    // double RHO             = vm["rho"].as<double>();
    // double RHODERIVATIVE   = vm["rhoderivative"].as<double>();

    // parse ngfspacing → SpacingType
    {
        auto s = vm["ngfspacing"].as<std::string>();
        std::replace(s.begin(), s.end(), ',', ' ');
        std::istringstream iss(s);
        std::vector<double> v{
            std::istream_iterator<double>(iss),
            std::istream_iterator<double>()};
        if(v.size() != ImageDimension)
        {
            std::cerr << "Error: ngfspacing must have "
                      << ImageDimension << " comma‐separated values\n";
            return EXIT_FAILURE;
        }
        MovingImageType::SpacingType ngf;
        for(unsigned i=0; i<ImageDimension; ++i) ngf[i] = v[i];
        vm.insert({ "parsed_ngfspacing",
            po::variable_value(boost::any(ngf),false) });
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
	typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;

	FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
	MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

	fixedImageReader->SetFileName(  fixedImageFN );
	movingImageReader->SetFileName( movingImageFN );

	fixedImageReader->Update();
	movingImageReader->Update();

	FixedImageType::ConstPointer fixedImage = fixedImageReader->GetOutput();

	MovingImageType::ConstPointer movingImage = movingImageReader->GetOutput();
	

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
	optimizerScales.Fill(1.0);
  optimizerScales[9] = GMT;
  optimizerScales[10] = GMT;
  optimizerScales[11] = GMT;
  optimizer->SetNumberOfIterations(NI);
  optimizer->SetMinimumStepLength(MINSTEP);
  optimizer->SetRelaxationFactor(RF);
  optimizer->SetGradientMagnitudeTolerance(GMT);
  optimizer->SetMaximumStepLength(MAXSTEP);

	

if (method == "translation") {
    // Only optimize translation parameters
	optimizerScales.Fill(1.0);
    for (int i = 0; i < 9; ++i) {
        optimizerScales[i] = GMT;
    }
} else if (method == "rotation") {
    // Only optimize rotation parameters
		optimizerScales.Fill(1.0);
    for (int i = 3; i < 12; ++i) {
        if (i < 9) {
            optimizerScales[i] = 1.0;
        } else {
            optimizerScales[i] = GMT;
        }
    }
} else if (method == "scaling") {
	optimizerScales.Fill(1.0);

    // Only optimize scaling parameters along the diagonal of the rotation matrix
    for (int i = 0; i < 12; ++i) {
        if (i == 0 || i == 4 || i == 8) {
            optimizerScales[i] = 1.0;
        } else {
            optimizerScales[i] = GMT;
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


	// itk::Point<double, 3> center;
	// itk::Index<3> centerIndex;

	// for (int i = 0; i < 3; ++i) {
	// 	centerIndex[i] = sourceImage->GetLargestPossibleRegion().GetSize()[i] / 2;
	// }
	// sourceImage->TransformIndexToPhysicalPoint(centerIndex, center);
	// transform->SetCenter(center);



	registration->SetInitialTransformParameters( transform->GetParameters() );

	const unsigned int numberOfPixels = fixedImage->GetLargestPossibleRegion().GetNumberOfPixels();
	const unsigned int numberOfSamplesMA  = static_cast<unsigned int>(numberOfPixels * MAPERCENTAGE);
    const unsigned int numberOfSamplesMSE = static_cast<unsigned int>(numberOfPixels * MSEPERCENTAGE);
    const unsigned int numberOfSamplesNGF = static_cast<unsigned int>(numberOfPixels * NGFPERCENTAGE);
    const unsigned int numberOfSamplesNC  = static_cast<unsigned int>(numberOfPixels * NCPERCENTAGE);

	metric->SetUseExplicitPDFDerivatives(EPDF);
	metric->SetNumberOfThreads(NT);
	metric->SetNormalizeDerivatives(NORMALIZE_DERIVATIVES);

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



    // apply parsed NGF‐spacing
    {
        auto ngf = boost::any_cast<MovingImageType::SpacingType>(
            vm["parsed_ngfspacing"].value());
        metric->SetNGFSpacing(ngf);
    }
	
	if (TR!=-99999999)
	{
		metric->SetUseFixedImageSamplesIntensityThreshold(TR);
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

	resample->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
	resample->SetOutputOrigin(  fixedImage->GetOrigin() );
	resample->SetOutputSpacing( fixedImage->GetSpacing() );
	resample->SetOutputDirection( fixedImage->GetDirection() );
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
