


//DeformableRegistration6
#include "itkImageRegistrationMethod.h"
#include "../../../Metrics/MANGF2ImageToImageMetric/itkMANGF2.h"
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
//#include "itkRescaleIntensityImageFilter.h"

//#include "itkShrinkImageFilter.h"

#include "itkTransformToDeformationFieldSource.h"

#include "../../Version.h"

#include "../../../include/funzioniIO.h"
#include "itkTransformFileWriter.h"

class CommandIterationUpdate : public itk::Command
{
public:
	typedef  CommandIterationUpdate   Self;
	typedef  itk::Command             Superclass;
	typedef itk::SmartPointer<Self>   Pointer;
	itkNewMacro( Self );
protected:
	CommandIterationUpdate() {};
public:
	typedef itk::LBFGSBOptimizer    OptimizerType;
	typedef   const OptimizerType * OptimizerPointer;
	void Execute(itk::Object *caller, const itk::EventObject & event)
	{
		Execute( (const itk::Object *)caller, event);
	}
	void Execute(const itk::Object * object, const itk::EventObject & event)
	{
		OptimizerPointer optimizer = static_cast< OptimizerPointer >( object );
		if( !(itk::IterationEvent().CheckEvent( &event )) )
		{
			return;
		}
		std::cout << optimizer->GetCurrentIteration() << "   ";
		std::cout << optimizer->GetCachedValue() << "   ";
		std::cout << optimizer->GetInfinityNormOfProjectedGradient() << std::endl;
	}
};

int main( int argc, char *argv[] )
{
	/**help function and parmaeters */
	std::vector<std::string>Vparams;
	Vparams.push_back("program name");
	int T=1;Vparams.push_back("fixedImageFile");
	int S=2;Vparams.push_back("movingImageFile");
	int O=3;Vparams.push_back("output Image");
	int OVF=4;Vparams.push_back("VF out");
	int NT=5; int nt=2; Vparams.push_back("Number of thread" + tostr(nt));readParametersFromStdin(&nt,NT,argv,argc,no);
	int BIN=NT+1; int NB=50;Vparams.push_back("Mattes number of bins");readParametersFromStdin(&NB,BIN,argv,argc,no);
	int MAPERCENTAGE=BIN+1; double nPCTMA=0.1; Vparams.push_back("Mattes percecntage" + tostr(nPCTMA));readParametersFromStdin(&nPCTMA,MAPERCENTAGE,argv,argc,no);
	int LAMBDA=MAPERCENTAGE+1; double lambda=1.0;Vparams.push_back("lambda values MI + lambda NGF" + tostr(lambda));readParametersFromStdin(&lambda,LAMBDA,argv,argc,no);
	int LAMBDADERIVATIVE=LAMBDA+1;double lambdaderivative=1e-4;Vparams.push_back("Lambda derivatuive " + tostr(lambdaderivative));readParametersFromStdin(&lambdaderivative,LAMBDADERIVATIVE,argv,argc,no);
	int ETAF=LAMBDADERIVATIVE+1;double etaF=2; Vparams.push_back("Eta values Reference (NGF noise)");readParametersFromStdin(&etaF,ETAF,argv,argc,no);
	int ETAM=ETAF+1;	double etaM=2; Vparams.push_back("Eta values Source (NGF noise)");readParametersFromStdin(&etaM,ETAM,argv,argc,no);
	int NGFEV=ETAM+1; int NGFEvaluator=0;Vparams.push_back("NGF Evaluator (0 scalar,1cross,2scdelta,3Delta,4Delta2)");	readParametersFromStdin(&NGFEvaluator,NGFEV,argv,argc,no);
	int GRID=NGFEV+1; double numberOfGridNodes=8;Vparams.push_back("Grid nodes: " + tostr(numberOfGridNodes) );	readParametersFromStdin(&numberOfGridNodes,GRID,argv,argc,no);
	int GRADCONV=GRID+1; double gradconv=0.005;Vparams.push_back("Grad conv: " + tostr(gradconv));	readParametersFromStdin(&gradconv,GRADCONV,argv,argc,no);
	int LINESEARCHACCURACY=GRADCONV+1;double linesearchaccuracy=0.9;Vparams.push_back("linesearchaccuracy: " + tostr(linesearchaccuracy));	readParametersFromStdin(&linesearchaccuracy,LINESEARCHACCURACY,argv,argc,no);
	int DEFAULTSTEPLENGTH=LINESEARCHACCURACY+1;double defaultStepLength=1.5;Vparams.push_back("default step legngth : " + tostr(defaultStepLength));	readParametersFromStdin(&defaultStepLength,DEFAULTSTEPLENGTH,argv,argc,no);
	int NI=DEFAULTSTEPLENGTH+1;int ni=1000;Vparams.push_back("Max number of Iterations: " + tostr(ni)); 	readParametersFromStdin(&ni,NI,argv,argc,no);
	int CFCF=NI+1;double cfcf=1.e12;Vparams.push_back("CostFunctionConvergenceFactor 1e+12 for low accuracy; 1e+7 for moderate accuracy and 1e+1 for extremely high accuracy." + tostr(cfcf)); 	readParametersFromStdin(&cfcf,CFCF,argv,argc,no);
	int PGT=CFCF+1;double pgt=1.e-5;Vparams.push_back("ProjectedGradientTolerance. Algorithm terminates when the project gradient is below the tolerance. Default value is 1e-5." + tostr(pgt)); 	readParametersFromStdin(&pgt,PGT,argv,argc,no);
	int NE=PGT+1;int ne=500;Vparams.push_back("Number of evaluations" + tostr(ne)); 	readParametersFromStdin(&ne,NE,argv,argc,no);
	int NC=NE+1;int nc=5;Vparams.push_back("Number of corrections" + tostr(nc)); 	readParametersFromStdin(&nc,NC,argv,argc,no);
	int TR=NC+1;int tr=1;Vparams.push_back("fixed Image threshold: " + tostr(tr)); 	readParametersFromStdin(&tr,TR,argv,argc,no);
	int BB=TR+1;bool TB=true;Vparams.push_back("Bspleine chaching 1 true ");readParametersFromStdin(&TB,BB,argv,argc,no);
	int BOUND=BB+1;int bound=0;Vparams.push_back(" Set the boundary condition for each variable, where = 0 if x[i] is unbounded, 1 if x[i] has only a lower bound,2 if x[i] has both lower and upper bounds and 3 if x[1] has only an upper bound: " + tostr(bound)); 	readParametersFromStdin(&bound,BOUND,argv,argc,no);
	int LBOUND=BOUND+1; double lbound=0;Vparams.push_back("Lower bound " + tostr(lbound));readParametersFromStdin(&lbound,LBOUND,argv,argc,no);
	int HBOUND=LBOUND+1; double hbound=0;Vparams.push_back("higher bound " + tostr(hbound));readParametersFromStdin(&hbound,HBOUND,argv,argc,no);
	int TOUT=HBOUND+1; Vparams.push_back("transforout.txt");
	int FMASK=TOUT+1;Vparams.push_back("Mask Reference");
	int MMASK=FMASK+1; Vparams.push_back("Mask source");


	if( argc < 4 )
	{
		std::cout<<argv[0]<<" ";
			printVersion();
			printParam(Vparams);
			std::cout << "Skip feature\t"<<no<< std::endl;
			std::cout << "Output and input type unsigned char \t"<< std::endl;


			firma();

		return EXIT_FAILURE;
	}


	

	const    unsigned int    ImageDimension = 3;
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


	typedef itk::MANGF2<
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
	registration->SetNumberOfThreads(nt);


	TransformType::Pointer  transform = TransformType::New();
	registration->SetTransform( transform );

	typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
	typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;

	FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
	MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

	fixedImageReader->SetFileName(  argv[T] );
	movingImageReader->SetFileName( argv[S] );

	fixedImageReader->Update();
	movingImageReader->Update();

	FixedImageType::ConstPointer fixedImage = fixedImageReader->GetOutput();
	





	registration->SetFixedImage(  fixedImage   );
	registration->SetMovingImage(   movingImageReader->GetOutput()   );



	FixedImageType::RegionType fixedRegion = fixedImage->GetBufferedRegion();

	registration->SetFixedImageRegion( fixedRegion );



	TransformType::PhysicalDimensionsType   fixedPhysicalDimensions;
	TransformType::MeshSizeType             meshSize;
	TransformType::OriginType               fixedOrigin;

	for( unsigned int i=0; i< SpaceDimension; i++ )
	{
		fixedOrigin[i] = fixedImage->GetOrigin()[i];
		fixedPhysicalDimensions[i] = fixedImage->GetSpacing()[i] *
				static_cast<double>(
						fixedImage->GetLargestPossibleRegion().GetSize()[i] - 1 );
	}
	meshSize.Fill( numberOfGridNodes - SplineOrder );

	transform->SetTransformDomainOrigin( fixedOrigin );
	transform->SetTransformDomainPhysicalDimensions(
			fixedPhysicalDimensions );
	transform->SetTransformDomainMeshSize( meshSize );
	transform->SetTransformDomainDirection( fixedImage->GetDirection() );


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
	const unsigned int numberOfSamplesMA =static_cast< unsigned int >( numberOfPixels * nPCTMA);

	metric->Setlambda(lambda);
	metric->SetlambdaDerivative(lambdaderivative);
	metric->SetBinNumbers(NB);
	metric->SetMAnumberOfSamples(numberOfSamplesMA);
	metric->SetNumberOfThreads(nt);
	metric->SetFixedEta(etaF);
	metric->SetMovingEta(etaM);
	metric->SetEvaluator(NGFEvaluator);
	metric->SetUseCachingOfBSplineWeights(TB);
	metric->SetUseFixedImageSamplesIntensityThreshold(tr);



	/*	optimizer->SetGradientConvergenceTolerance( gradconv1 );
	optimizer->SetLineSearchAccuracy( linesearch1 );
	optimizer->SetDefaultStepLength( defaultStepLength );
	optimizer->TraceOn();
	optimizer->SetMaximumNumberOfFunctionEvaluations( ni );
	 */


	const unsigned int numParameters = transform->GetNumberOfParameters();
	OptimizerType::BoundSelectionType boundSelect( numParameters );
	OptimizerType::BoundValueType upperBound( numParameters );
	OptimizerType::BoundValueType lowerBound( numParameters );

	boundSelect.Fill( bound );
	upperBound.Fill( hbound );
	lowerBound.Fill( lbound );

	optimizer->SetBoundSelection( boundSelect );
	optimizer->SetUpperBound( upperBound );
	optimizer->SetLowerBound( lowerBound );
	//CostFunctionConvergenceFactor 1e+12 for low accuracy; 1e+7 for moderate accuracy and 1e+1 for extremely high accuracy.
	optimizer->SetCostFunctionConvergenceFactor( cfcf );
	optimizer->SetProjectedGradientTolerance( pgt );
	optimizer->SetMaximumNumberOfIterations( ni );
	optimizer->SetMaximumNumberOfEvaluations( ne );
	optimizer->SetMaximumNumberOfCorrections( nc);

	std::cout<<"\nUpper"<<upperBound<<std::endl;

	std::cout<<"\nLower"<<lowerBound<<std::endl;

	optimizer->SetMinimize(true);




	std::cout << "Starting Registration "
			<< std::endl;

	optimizer->SetMinimize((bool)1);
	registration->SetNumberOfThreads(nt);



	if(argc>FMASK)
		if(atof(argv[FMASK])!=no)
		{
			typedef itk::ImageMaskSpatialObject< ImageDimension > MaskType;
			MaskType::Pointer spatialObjectMask = MaskType::New();
			typedef itk::Image< unsigned char, ImageDimension > ImageMaskType;
			typedef itk::ImageFileReader< ImageMaskType > MaskReaderType;
			MaskReaderType::Pointer maskReader = MaskReaderType::New();
			maskReader->SetFileName( argv[FMASK] );
			try
			{
				maskReader->Update();
			}
			catch( itk::ExceptionObject & err )
			{
				std::cerr << "ExceptionObject caught !" << std::endl;
				std::cerr << err << std::endl;
				return EXIT_FAILURE;
			}
			spatialObjectMask->SetImage( maskReader->GetOutput() );
			metric->SetFixedImageMask( spatialObjectMask );
			std::cerr<<"FixedMASK: "<< argv[FMASK]<<std::endl;
		}




	if(argc>MMASK)
		if (atof(argv[MMASK])!=no)
		{
			typedef itk::ImageMaskSpatialObject< ImageDimension > MaskType;
			MaskType::Pointer spatialObjectMask = MaskType::New();
			typedef itk::Image< unsigned char, ImageDimension > ImageMaskType;
			typedef itk::ImageFileReader< ImageMaskType > MaskReaderType;
			MaskReaderType::Pointer maskReader = MaskReaderType::New();
			maskReader->SetFileName( argv[MMASK] );
			try
			{
				maskReader->Update();
			}
			catch( itk::ExceptionObject & err )
			{
				std::cerr << "ExceptionObject caught !" << std::endl;
				std::cerr << err << std::endl;
				return EXIT_FAILURE;
			}
			spatialObjectMask->SetImage( maskReader->GetOutput() );
			metric->SetMovingImageMask( spatialObjectMask );
			std::cerr<<"MovingMASK: "<< argv[MMASK]<<std::endl;
		}




	std::cout<< "\n\n\n\n FFD MANGF	\n\tThread: "<< metric->GetNumberOfThreads() <<
			"\n\tlambda: " << metric->Getlambda() <<
			"\n\tlambda derivative: " << metric->GetlambdaDerivative() <<
			"\n\tBin: " <<metric->GetBinNumbers() <<
			"\n\tETAF: " <<metric->GetFixedEta() <<
			"\n\t ETAM: " <<metric->GetMovingEta() <<
			"\n\tMA Pixels Number: "<<metric->GetMAnumberOfSamples()<<
			"\n\tNGF Pixels Number: "<<metric->GetNGFnumberOfSamples()<<
			"\n\tNGF Evaluator: " << (int) metric->GetEvaluator()<<
			"\n\tpixels in the image: " << fixedImageReader->GetOutput()->GetLargestPossibleRegion().GetNumberOfPixels()<<
			"\nVAriables: " <<transform->GetNumberOfParameters() <<
			std::endl;




	CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
	optimizer->AddObserver( itk::IterationEvent(), observer );



	try
	{
		registration->Update();
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
	resample->SetDefaultPixelValue( 0 );

	typedef short  OutputPixelType;

	typedef itk::Image< OutputPixelType, ImageDimension > OutputImageType;

	typedef itk::CastImageFilter<
			FixedImageType,
			OutputImageType > CastFilterType;

	typedef itk::ImageFileWriter< OutputImageType >  WriterType;


	WriterType::Pointer      writer =  WriterType::New();
	CastFilterType::Pointer  caster =  CastFilterType::New();


	writer->SetFileName( argv[O] );


	caster->SetInput( resample->GetOutput() );
	writer->SetInput( caster->GetOutput()   );


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


	if (argc>OVF)
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
		rtd->SetFileName(argv[OVF]);
		rtd->Update();
	};



	return EXIT_SUCCESS;
}
