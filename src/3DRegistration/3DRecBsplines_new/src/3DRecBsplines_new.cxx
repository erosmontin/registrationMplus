


//DeformableRegistration6
#include "itkImageRegistrationMethod.h"
#include "../../../Metrics/Mplus/itkMANGFMSE.h"
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
#include "../../Version.h"


#include "../../../includes/imageUtils.h"
#include "../../../includes/registrationUtils.h"

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;



const unsigned int ImageDimension = 3;

int main( int argc, char *argv[] )
{
    po::options_description desc("B-spline Registration\n"
	"Dr. Eros Montin Ph.D., 2014\n"
	"eros.montin@gmail.com\n\n"
	"cite us:\n\nMontin, E., Belfatto, A., Bologna, M., Meroni, S., Cavatorta, C., Pecori, E., Diletto, B., Massimino, M., Oprandi, M. C., Poggi, G., Arrigoni, F., Peruzzo, D., Pignoli, E., Gandola, L., Cerveri, P., & Mainardi, L. (2020). A multi-metric registration strategy for the alignment of longitudinal brain images in pediatric oncology. Medical & biological engineering & computing, 58(4), 843–855. https://doi.org/10.1007/s11517-019-02109-4\n\n"
	"Allowed options for mu MI + lambda NGF +  nu MSE");
    desc.add_options()
	    ("help,h", "produce help message")
        ("fixedimage,f", po::value<std::string>(), "Fixed image filename")
        ("movingimage,m", po::value<std::string>(), "Moving image filename")
        ("outputimage,o", po::value<std::string>(), "Output registered imagefilename")
        ("vfout,v", po::value<std::string>()->default_value("N"), "VF output filename")
        ("numberofthreads", po::value<int>()->default_value(2), "Number of threads 2")
	    ("mu", po::value<double>()->default_value(1.0), "mu value MI 1.0")
		("muderivative,M", po::value<double>()->default_value(1.0), "mu derivative 1.0")
		("mattespercentage,p", po::value<double>()->default_value(0.1), "Mattes percentage 0.1")
        ("mattesnumberofbins,b", po::value<int>()->default_value(64), "Mattes number of bins 64")
	    ("bsplinecaching,B", po::value<bool>()->default_value(true), "Bspline caching, 1 for true")
		("explicitPDFderivatives", po::value<bool>()->default_value(false), "Explicit PDF derivatives, 0 for false")
		("lambda,l", po::value<double>()->default_value(1.0), "lambda value MI + lambda NGF 1.0")
        ("lambdaderivative,L", po::value<double>()->default_value(0), "Lambda derivative 0 no derivatives")
        ("etavaluefixed,r", po::value<double>()->default_value(2), "Eta value fixed image(NGF noise) 2")
        ("etavaluemoving,s", po::value<double>()->default_value(2), "Eta value moving image (NGF noise) 2")
        ("NGFevaluator", po::value<int>()->default_value(0), "NGF Evaluator (0 scalar,1cross,2scdelta,3Delta,4Delta2)")
	    ("nu,n", po::value<double>()->default_value(1.0), "yota value MSE 1.0")
		("nuderivative,N", po::value<double>()->default_value(1.0), "nu derivative 1.0")
        ("gridresolution,g", po::value<double>()->default_value(2), "Mesh resolution (mm) 2")
        ("maxnumberofiterations,I", po::value<int>()->default_value(1000), "Max number of Iterations 1000")
        ("costfunctionconvergencefactor,F", po::value<double>()->default_value(1.e12), "CostFunctionConvergenceFactor 1e+12 for low accuracy; 1e+7 for moderate accuracy and 1e+1 for extremely high accuracy.")
        ("projectedgradienttolerance,P", po::value<double>()->default_value(1.e-5), "ProjectedGradientTolerance. Algorithm terminates when the project gradient is below the tolerance. Default value is 1e-5.")
        ("numberofevaluations,E", po::value<int>()->default_value(500), "Number of Evaluations")
        ("numberofcorrections,C", po::value<int>()->default_value(5), "Number of Corrections")
		("fixedimagethreshold,t", po::value<double>()->default_value(-99999999), "Fixed image threshold")
    
        ("bound", po::value<int>()->default_value(0), "Set the boundary condition for each variable, where = 0 if x[i] is unbounded, 1 if x[i] has only a lower bound,2 if x[i] has both lower and upper bounds and 3 if x[1] has only an upper bound")
        ("lbound", po::value<double>()->default_value(0), "Lower bound")
        ("ubound", po::value<double>()->default_value(0), "Upper bound")
        ("transformout,T", po::value<std::string>()->default_value("N"), "Output for transform")
		("gridposition,G", po::value<std::string>()->default_value("N"), "Read the position of the grid from a file")
		("dfltpixelvalue,P", po::value<double>()->default_value(0), "Default pixel value")
		("verbose,V", po::value<bool>()->default_value(false), "verbose")

    ;
	


    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);


	if (vm.count("help") || !vm.count("fixedimage") || !vm.count("movingimage") || !vm.count("outputimage") || !vm.count("vfout")){
    std::cout << desc << "\n";
    return 1;
	}

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

    std::string fixedImageFN = vm["fixedimage"].as<std::string>();
    std::string movingImageFN = vm["movingimage"].as<std::string>();

	std::string ou = vm["outputimage"].as<std::string>();
	std::string ouvf = vm["vfout"].as<std::string>();

	int NT=vm["numberofthreads"].as<int>();
	int NB=vm["mattesnumberofbins"].as<int>();
    double MAPERCENTAGE = vm["mattespercentage"].as<double>();
	double MU = vm["mu"].as<double>();
	double MUDERIVATIVE = vm["muderivative"].as<double>();
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
	double GRIDRESOLUTION=vm["gridresolution"].as<double>();
	int NI=vm["maxnumberofiterations"].as<int>();
	double CFCF=vm["costfunctionconvergencefactor"].as<double>();
	double PGT=vm["projectedgradienttolerance"].as<double>();
	int NE=vm["numberofevaluations"].as<int>();
	int NC=vm["numberofcorrections"].as<int>();
	double TR=vm["fixedimagethreshold"].as<double>();
	bool TB=vm["bsplinecaching"].as<bool>();
	bool EPDF=vm["explicitPDFderivatives"].as<bool>();
	int BOUND=vm["bound"].as<int>();
	double LBOUND=vm["lbound"].as<double>();
	double UBOUND=vm["ubound"].as<double>();
	std::string TOUT=vm["transformout"].as<std::string>();
	double DFLTPIXELVALUE=vm["dfltpixelvalue"].as<double>();
	bool V=vm["verbose"].as<bool>();
	std::string GRIDPOSITION=vm["gridposition"].as<std::string>();

	typedef  float          PixelType;
	typedef itk::Image< PixelType, ImageDimension >  FixedImageType;
	typedef itk::Image< PixelType, ImageDimension >  MovingImageType;

	const unsigned int SpaceDimension = ImageDimension;
	const unsigned int SplineOrder = 3;
	typedef double CoordinateRepType;

	typedef itk::BSplineTransform<
			CoordinateRepType,
			SpaceDimension,
			SplineOrder >     TransformType;


	typedef itk::LBFGSBOptimizer       OptimizerType;


	typedef itk::MANGFMSE<
			FixedImageType,
			MovingImageType >    MetricType;

	typedef itk:: LinearInterpolateImageFunction<
			MovingImageType,
			double          >    InterpolatorType;

	typedef itk::ImageRegistrationMethod<
			FixedImageType,
			MovingImageType >    RegistrationType;

	MetricType::Pointer         metric        = MetricType::New();
	OptimizerType::Pointer      optimizer     = OptimizerType::New();
	InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
	RegistrationType::Pointer   registration  = RegistrationType::New();

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

	TransformType::PhysicalDimensionsType   fixedPhysicalDimensions;
	TransformType::MeshSizeType             meshSize;
	TransformType::OriginType               fixedOrigin;


	FixedImageType::SpacingType meshspacing = fixedImage->GetSpacing();
	FixedImageType::PointType meshorigin = fixedImage->GetOrigin();
	FixedImageType::DirectionType meshdirection = fixedImage->GetDirection();
	FixedImageType::SizeType meshsize = fixedImage->GetLargestPossibleRegion().GetSize();


	if (GRIDPOSITION!="N"){
		//read the image that specifies the position of the grid
		FixedImageReaderType::Pointer  meshImageReader  = FixedImageReaderType::New();
		meshImageReader->SetFileName(  GRIDPOSITION );
		meshImageReader->Update();
		FixedImageType::ConstPointer meshImage = meshImageReader->GetOutput();
		//get the information of the image that specifies the position of the grid and overwrite the information of the fixed image
		meshspacing=meshImage->GetSpacing();
		meshorigin=meshImage->GetOrigin();
		meshdirection=meshImage->GetDirection();
		meshsize=meshImage->GetLargestPossibleRegion().GetSize();

		// resample the meshimage on the fixedimage space in casse the two images have different size

		typedef itk::ResampleImageFilter<FixedImageType, FixedImageType> ResampleFilterType;
		ResampleFilterType::Pointer resampler = ResampleFilterType::New();
		resampler->SetInput(meshImage);
		resampler->SetOutputParametersFromImage(fixedImage);
		resampler->Update();
		FixedImageType::RegionType meshregionresampled=resampler->GetOutput()->GetLargestPossibleRegion();
		
		fixedImage->SetRequestedRegion(meshregionresampled);
	}

	
	unsigned int numberOfGridNodes=0;
	for( unsigned int i=0; i< SpaceDimension; i++ )
	{
		fixedOrigin[i] = meshorigin[i];
		fixedPhysicalDimensions[i] = meshspacing[i] *
				static_cast<double>(
						meshsize[i] - 1 );
		numberOfGridNodes=static_cast<int>((fixedPhysicalDimensions[i]-fixedOrigin[i])/GRIDRESOLUTION);
		    if (numberOfGridNodes <= SplineOrder) {
        std::cerr << "Error: numberOfGridNodes must be greater than 0" << std::endl;
        return EXIT_FAILURE;
    	}

		meshSize[i] =  numberOfGridNodes - SplineOrder;
    std::cout << "Dimension " << i << ": numberOfGridNodes = " << numberOfGridNodes << ", meshSize = " << meshSize[i] << std::endl;

	}

	transform->SetTransformDomainOrigin( fixedOrigin );
	transform->SetTransformDomainPhysicalDimensions(
			fixedPhysicalDimensions );
	transform->SetTransformDomainMeshSize( meshSize );
	transform->SetTransformDomainDirection( fixedImage->GetDirection() );
	transform->SetIdentity();


	//metric is negative
	optimizer->MinimizeOn();
	typedef TransformType::ParametersType     ParametersType;

	const unsigned int numberOfParameters =
			transform->GetNumberOfParameters();

	ParametersType parametersLow( numberOfParameters );

	parametersLow.Fill( 0.0 );

	transform->SetParameters( parametersLow );

	registration->SetInitialTransformParameters( transform->GetParameters() );

	const unsigned int numberOfPixels = fixedImage->GetLargestPossibleRegion().GetNumberOfPixels();
	const unsigned int numberOfSamplesMA =static_cast< unsigned int >( numberOfPixels * MAPERCENTAGE );

	metric->SetMu(MU);
	metric->SetMuDerivative(MUDERIVATIVE);
	metric->SetMSENumberOfSamples(numberOfSamplesMA);
	metric->SetBinNumbers(NB);
	metric->SetMANumberOfSamples(numberOfSamplesMA);
	metric->SetUseCachingOfBSplineWeights(TB);
	metric->SetUseExplicitPDFDerivatives(EPDF);
	metric->SetNu(NU);
	metric->SetNuDerivative(NUDERIVATIVE);
	metric->SetLambda(LAMBDA);
	metric->SetLambdaDerivative(LAMBDADERIVATIVE);
	metric->SetNGFNumberOfSamples(numberOfSamplesMA);

	metric->SetNumberOfThreads(NT);
	metric->SetFixedEta(ETAF);
	metric->SetMovingEta(ETAM);
	metric->SetEvaluator(NGFevaluator);
	
	if (TR!=-99999999)
	{
		metric->SetUseFixedImageSamplesIntensityThreshold(TR);
	}

	const unsigned int numParameters = transform->GetNumberOfParameters();
	OptimizerType::BoundSelectionType boundSelect( numParameters );
	OptimizerType::BoundValueType upperBound( numParameters );
	OptimizerType::BoundValueType lowerBound( numParameters );

	boundSelect.Fill( BOUND );
	upperBound.Fill( UBOUND );
	lowerBound.Fill(  LBOUND);

	optimizer->SetBoundSelection( boundSelect );
	optimizer->SetUpperBound( upperBound );
	optimizer->SetLowerBound( lowerBound );
	//CostFunctionConvergenceFactor 1e+12 for low accuracy; 1e+7 for moderate accuracy and 1e+1 for extremely high accuracy.
	optimizer->SetCostFunctionConvergenceFactor( CFCF );

	if (V){
		optimizer->TraceOn();
	}

	optimizer->SetProjectedGradientTolerance( PGT );
	optimizer->SetMaximumNumberOfIterations( NI );
	optimizer->SetMaximumNumberOfEvaluations( NE );
	optimizer->SetMaximumNumberOfCorrections( NC);

	std::cout<<"\nUpper"<<upperBound[0]<<std::endl;

	std::cout<<"\nLower"<<lowerBound[0]<<std::endl;

	optimizer->SetMinimize(true);


	std::cout << "Starting Registration "
			<< std::endl;


	std::cout<< "\n\n\n\n B-spline transform usin itkMANGFMSE	\n\tThread: "<< metric->GetNumberOfThreads() <<
			"\nVariables: " <<transform->GetNumberOfParameters() <<
			"\nVariables: " << transform->GetTransformDomainMeshSize() <<
			
			std::endl;


	
	// size_t number_of_points = meshSize[0] * meshSize[1] * meshSize[2];
	// int number_of_cells = transform->GetNumberOfParameters();

	// // Estimate memory usage
	// using CellType = MeshType::CellType;
	// using CellAutoPointer = CellType::CellAutoPointer;

	// CellAutoPointer cellPointer;
	// size_t number_of_points_per_cell = cellPointer->GetNumberOfPoints();
	// long int  memory_for_points = number_of_points * ImageDimension * sizeof(float);
	// long int memory_for_cells = number_of_cells * (1 + number_of_points_per_cell) * sizeof(int);

	// size_t total_memory = memory_for_points + memory_for_cells;

	// std::cout<<totak_memory<<std::endl;


	LBFGSBOptimizeCommandIterationUpdate::Pointer observer = LBFGSBOptimizeCommandIterationUpdate::New();
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


	if (ouvf!="N")
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
		rtd->SetFileName(ouvf);
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
