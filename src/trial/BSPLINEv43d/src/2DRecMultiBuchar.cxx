


//DeformableRegistration6
#include "itkImageRegistrationMethod.h"
#include "../../../Metrics/MANGFImageToImageMetric/itkMANGF.h"
#include "itkBSplineTransform.h"
#include  "itkSPSAOptimizer.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkBSplineResampleImageFunction.h"
#include "itkIdentityTransform.h"
#include "itkBSplineDecompositionImageFilter.h"
#include "itkImageMaskSpatialObject.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkShrinkImageFilter.h"

#include "itkTransformToDeformationFieldSource.h"


int main( int argc, char *argv[] )
{
int T=1;
	int S=2;
	int O=3;
	int OVF=4;
	int THREAD=5; int NT=4;
	int BIN=6; int NB=32;
	int MAPERCENTAGE=7; double nPCTMA=0.1;
	int LAMBDA=MAPERCENTAGE+1; double lambda=1;
	int LAMBDADERIVATIVE=LAMBDA+1; double lambdaderivative=0.001;
	int ETAF=LAMBDADERIVATIVE+1; double etaF= 5;
	int ETAM=ETAF+1; double etaM= 5;
	int NGFEV=ETAM+1;	int NGFEvaluator=0;//scalar
	int GRID1=NGFEV+1; double numberOfGridNodesCoarse=8;
	int GRADIENTCONVERGENCETOL1=GRID1+1; double gradconv1=0.005;
	int LINEARSEARCHACURACY1=GRADIENTCONVERGENCETOL1+1;double linesearch1=0.9;
	int DEFAULTSTEPLENGTH1=LINEARSEARCHACURACY1+1;double stepl1=1.5;
	int MAXITERATIONS1=DEFAULTSTEPLENGTH1+1;int NI1=1000;
	int BD=MAXITERATIONS1+1;bool TD=1;
	int BB=MAXITERATIONS1+2;bool TB=1;
	int FMASK=BB+1;
	int MMASK=FMASK+1;
	int no=-1;


	if( argc < 4 )
	{
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << "\n01 fixedImageFile\n02 movingImageFile\n03 outputImagefile  ";
		std::cerr << " \n04 deformationField \n05Thread ";
		std::cerr << " \n06 BIN \n07 MApercentage1"
					"\n08 Lambda \n09 Derivative lambda"
					 "\n10 ETAF \n11 ETAM \n12NGF Evaluator (0 scalar,1cross,2scdelta,3Delta,4Delta2)"
				"\n13 numberof meshpoint per axis "
				"\n14 gradient convergencetol1 0.005"
				"\n15 linearsearchaccuracy1 0.9 \n16 step length 1.5 "
				"\n17 maxiteraiontio1 \n18 BD \n19BB \n20FMASK \n21MMASK \n eros.montin@polimi.it";
		return EXIT_FAILURE;
	}


	if(argc>THREAD)
		if(atoi(argv[THREAD])!=no){NT=atoi(argv[THREAD]);};
	if(argc>GRID1)
		if(atoi(argv[GRID1])!=no){numberOfGridNodesCoarse=atoi(argv[GRID1]);};
	if(argc>GRADIENTCONVERGENCETOL1)
		if(atoi(argv[GRADIENTCONVERGENCETOL1])!=no){gradconv1=atof(argv[GRADIENTCONVERGENCETOL1]);};
	if(argc>LINEARSEARCHACURACY1)
		if(atoi(argv[LINEARSEARCHACURACY1])!=no){linesearch1=atof(argv[LINEARSEARCHACURACY1]);};
	if(argc>DEFAULTSTEPLENGTH1)
		if(atoi(argv[DEFAULTSTEPLENGTH1])!=no){stepl1=atof(argv[DEFAULTSTEPLENGTH1]);};
	if(argc>MAXITERATIONS1)
		if(atoi(argv[MAXITERATIONS1])!=no){NI1=atoi(argv[MAXITERATIONS1]);};

		if (argc>BIN)
		if(atof(argv[BIN])!=no){NB=atoi(argv[BIN]);}

	if (argc>MAPERCENTAGE)
		if( atof(argv[MAPERCENTAGE])!=no){nPCTMA=atof(argv[MAPERCENTAGE]);}
		if (argc>LAMBDA)
		if( atof(argv[LAMBDA])!=no){lambda=atof(argv[LAMBDA]);}

	if (argc>LAMBDADERIVATIVE)
			if( atof(argv[LAMBDADERIVATIVE])!=no){lambdaderivative=atof(argv[LAMBDADERIVATIVE]);}

	if (argc>ETAF)
			if( atof(argv[ETAF])!=no){etaF=atof(argv[ETAF]);}

	if (argc>ETAM)
			if( atof(argv[ETAM])!=no){etaM=atof(argv[ETAM]);}

	if (argc>NGFEV)
		if( atof(argv[NGFEV])!=no){NGFEvaluator=atof(argv[NGFEV]);}






	const    unsigned int    ImageDimension = 3;
	typedef  unsigned int          PixelType;

	typedef itk::Image< PixelType, ImageDimension >  FixedImageType;
	typedef itk::Image< PixelType, ImageDimension >  MovingImageType;

	const unsigned int SpaceDimension = ImageDimension;
	const unsigned int SplineOrder = 3;
	typedef double CoordinateRepType;

	typedef itk::BSplineTransform<
			CoordinateRepType,
			SpaceDimension,
			SplineOrder >     TransformType;


	 typedef itk::SPSAOptimizer       OptimizerType;


	typedef itk::MANGF<
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
	meshSize.Fill( numberOfGridNodesCoarse - SplineOrder );

	transform->SetTransformDomainOrigin( fixedOrigin );
	transform->SetTransformDomainPhysicalDimensions(
			fixedPhysicalDimensions );
	transform->SetTransformDomainMeshSize( meshSize );
	transform->SetTransformDomainDirection( fixedImage->GetDirection() );


		//metric is negative
		//optimizer->MinimizeOn();
	typedef TransformType::ParametersType     ParametersType;

	const unsigned int numberOfParameters =
			transform->GetNumberOfParameters();

	ParametersType parametersLow( numberOfParameters );

	parametersLow.Fill( 0.0 );

	transform->SetParameters( parametersLow );


	registration->SetInitialTransformParameters( transform->GetParameters() );

	const unsigned int numberOfPixels = fixedImage->GetLargestPossibleRegion().GetNumberOfPixels();




const unsigned int numberOfSamplesMA =static_cast< unsigned int >( numberOfPixels * nPCTMA);

	const unsigned int numberOfSamplesNGF =static_cast< unsigned int >( numberOfPixels);



	metric->Setlambda(lambda);
	metric->SetlambdaDerivative(lambdaderivative);
	metric->SetBinNumbers(NB);
	metric->SetMAnumberOfSamples(numberOfSamplesMA);
	metric->SetNGFnumberOfSamples(numberOfSamplesNGF);
	metric->SetNumberOfThreads(NT);
	metric->SetFixedEta(etaF);
	metric->SetMovingEta(etaM);
	metric->SetEvaluator(NGFEvaluator);

	optimizer->SetMinimize(1);

/*	optimizer->SetGradientConvergenceTolerance( gradconv1 );
	optimizer->SetLineSearchAccuracy( linesearch1 );
	optimizer->SetDefaultStepLength( stepl1 );
	optimizer->TraceOn();
	optimizer->SetMaximumNumberOfFunctionEvaluations( NI1 );
	*/


	/* const unsigned int numParameters = transform->GetNumberOfParameters();
		  OptimizerType::BoundSelectionType boundSelect( numParameters );
		  OptimizerType::BoundValueType upperBound( numParameters );
		  OptimizerType::BoundValueType lowerBound( numParameters );
		  boundSelect.Fill( 0 );
		  upperBound.Fill( 0.0 );
		  lowerBound.Fill( 0.0 );
		  optimizer->SetBoundSelection( boundSelect );
		  optimizer->SetUpperBound( upperBound );
		  optimizer->SetLowerBound( lowerBound );
		  //1e+12 for low accuracy; 1e+7 for moderate accuracy and 1e+1 for extremely high accuracy.
		  optimizer->SetCostFunctionConvergenceFactor( 1.e7 );
		  optimizer->SetProjectedGradientTolerance( 1e-6 );
		  optimizer->SetMaximumNumberOfIterations( NI1 );
		  optimizer->SetMaximumNumberOfEvaluations( 30 );
		  optimizer->SetMaximumNumberOfCorrections( 5 );


*/

	std::cout << "Starting Registration "
			<< std::endl;

	//optimizer->SetMinimize((bool)1);
	registration->SetNumberOfThreads(NT);



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

	typedef  unsigned short  OutputPixelType;

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
			td->SetOutputParametersFromImage(fixedImageReader->GetOutput());
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
