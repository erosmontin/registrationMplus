//DeformableRegistration6
#include "itkImageRegistrationMethod.h"
#include "/DATA/Dropbox/itkSRC/DONE/MANGFImageToImageMetric/itkMANGF.h"
#include "itkBSplineTransform.h"
#include "itkLBFGSOptimizer.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkBSplineResampleImageFunction.h"
#include "itkIdentityTransform.h"
#include "itkBSplineDecompositionImageFilter.h"
#include "itkImageMaskSpatialObject.h"


#include "itkTransformToDeformationFieldSource.h"

// NOTE: the LBFGSOptimizer does not invoke events
int main( int argc, char *argv[] )
{
int T=1;
	int S=2;
	int O=3;
	int OVF=4;
	int THREAD=5; int NT=1;
	int BIN1=THREAD+1; int BN1=50;
	int MAPERCENTAGE1=BIN1+1; double nPCT1MA=0.01;
	int NGFPERCENTAGE1=MAPERCENTAGE1+1; double nPCT1NGF=0.01;
	int NOISE= NGFPERCENTAGE1+1;double noise=0.1;
	int LAMBDA=NGFPERCENTAGE1+2; double lambda=1;
	int LAMBDAD=LAMBDA+1; double lambdaDerivative=1;
	int NGFEV=LAMBDAD+1; char NGFEvaluator=0;//scalar
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
		std::cerr << " \n04 deformationField \nThread ";
		std::cerr << " \n BIN \n MApercentage1 \n NGFpercentage"
					"\n eta \n Lambda \n Derivative lambda\n NGF Evaluator (0 scalar,1cross,2scdelta,3Delta,4Delta2)"
				"\n numberof meshpoint per axis "
				"\n gradient convergencetol1 0.005"
				"\n linearsearchaccuracy1 0.9 \n step length 1.5 "
				"\n maxiteraiontio1 \n BD \nBB \nFMASK \nMMASK \n eros.montin@polimi.it";
		return EXIT_FAILURE;
	}


	if(argc>THREAD)
		if(atoi(argv[THREAD])!=no){NT=atoi(argv[THREAD]);};
	if(argc>BIN1)
		if(atoi(argv[BIN1])!=no){BN1=atoi(argv[BIN1]);};
	if(argc>MAPERCENTAGE1)
		if(atof(argv[MAPERCENTAGE1])!=no){nPCT1MA=atof(argv[MAPERCENTAGE1]);};
	if(argc>NGFPERCENTAGE1)
		if(atof(argv[NGFPERCENTAGE1])!=no){nPCT1NGF=atof(argv[NGFPERCENTAGE1]);};
	if (argc>LAMBDA)
		if( atof(argv[LAMBDA])!=no){lambda=atof(argv[LAMBDA]);}
	if (argc>LAMBDAD)
			if( atof(argv[LAMBDAD])!=no){lambdaDerivative=atof(argv[LAMBDAD]);}
	if (argc>NGFEV)
		if( atof(argv[NGFEV])!=no){NGFEvaluator=atof(argv[NGFEV]);}
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
	if(argc>NOISE)
			if(atoi(argv[NOISE])!=no){noise=atof(argv[NOISE]);};

	const    unsigned int    ImageDimension = 2;
	typedef  float           PixelType;

	typedef itk::Image< PixelType, ImageDimension >  FixedImageType;
	typedef itk::Image< PixelType, ImageDimension >  MovingImageType;

	const unsigned int SpaceDimension = ImageDimension;
	const unsigned int SplineOrder = 3;
	typedef double CoordinateRepType;

	typedef itk::BSplineTransform<
			CoordinateRepType,
			SpaceDimension,
			SplineOrder >     TransformType;


	typedef itk::LBFGSOptimizer       OptimizerType;


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


	TransformType::Pointer  transformLow = TransformType::New();
	registration->SetTransform( transformLow );

	typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
	typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;

	FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
	MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

	fixedImageReader->SetFileName(  argv[T] );
	movingImageReader->SetFileName( argv[S] );

	FixedImageType::ConstPointer fixedImage = fixedImageReader->GetOutput();

	registration->SetFixedImage(  fixedImage   );
	registration->SetMovingImage(   movingImageReader->GetOutput()   );

	fixedImageReader->Update();

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

	transformLow->SetTransformDomainOrigin( fixedOrigin );
	transformLow->SetTransformDomainPhysicalDimensions(
			fixedPhysicalDimensions );
	transformLow->SetTransformDomainMeshSize( meshSize );
	transformLow->SetTransformDomainDirection( fixedImage->GetDirection() );


		//metric is negative
		optimizer->MinimizeOn();
	typedef TransformType::ParametersType     ParametersType;

	const unsigned int numberOfParameters =
			transformLow->GetNumberOfParameters();

	ParametersType parametersLow( numberOfParameters );

	parametersLow.Fill( 0.0 );

	transformLow->SetParameters( parametersLow );


	registration->SetInitialTransformParameters( transformLow->GetParameters() );

	const unsigned int TnumberOfPixels = fixedImage->GetLargestPossibleRegion().GetNumberOfPixels();

	const unsigned int TnumberOfSamples1MA =
			static_cast< unsigned int >( TnumberOfPixels * nPCT1MA);

	const unsigned int TnumberOfSamples1NGF =
			static_cast< unsigned int >( TnumberOfPixels * nPCT1NGF);





	metric->SetS(lambda);
	metric->SetMABinNumbers(BN1);
	metric->SetMAnumberOfSamples(TnumberOfSamples1MA);
	metric->SetNGFnumberOfSamples(TnumberOfSamples1NGF);
	metric->SetNumberOfThreads(NT);
	metric->SetEvaluator(NGFEvaluator);
	metric->SetNoise(noise);
	metric->SetSDerivative(lambdaDerivative);


	optimizer->SetGradientConvergenceTolerance( gradconv1 );
	optimizer->SetLineSearchAccuracy( linesearch1 );
	optimizer->SetDefaultStepLength( stepl1 );
	optimizer->TraceOn();
	optimizer->SetMaximumNumberOfFunctionEvaluations( NI1 );
	std::cout << "Starting Registration with low resolution transform"
			<< std::endl;

	optimizer->SetMinimize((bool)1);
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




	std::cout<< "\n\n\n\nMattes MI"
			"\n\tThread: "<< metric->GetNumberOfThreads() <<
								  "\n\tlambda: " << metric->GetS() <<
								  "\n\tlambda derivative: " << metric->GetSDerivative() <<
								  "\n\tNoise: " << metric->GetNoise() <<
								  "\n\tBin: " <<metric->GetMABinNumbers() <<
								   "\n\tMA Pixels Number: "<<metric->GetMAnumberOfSamples()<<
								  "\n\tNGF Pixels Number: "<<metric->GetNGFnumberOfSamples()<<
								  "\n\tNGF Evaluator: " << (int) metric->GetEvaluator()<<
								  "\n\tpixels in the image: " << fixedImageReader->GetOutput()->GetLargestPossibleRegion().GetNumberOfPixels()<<
			"\nVAriables: " <<transformLow->GetNumberOfParameters() <<
			"\nOptimizer conv tol: "<< optimizer->GetGradientConvergenceTolerance()<<
			"\noptimzer linear search accuracy" <<optimizer->GetLineSearchAccuracy()<<
			"\nstep length: " <<optimizer->GetDefaultStepLength()<<
			"\nMax iterations: " <<optimizer->GetMaximumNumberOfFunctionEvaluations()<<
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


	transformLow->SetParameters( registration->GetLastTransformParameters() );

	typedef itk::ResampleImageFilter<
			MovingImageType,
			FixedImageType >    ResampleFilterType;

	ResampleFilterType::Pointer resample = ResampleFilterType::New();

	resample->SetTransform( transformLow );
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
