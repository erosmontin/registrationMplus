#include "itkImageRegistrationMethod.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "../../../Metrics/MANGFImageToImageMetric/Code/itkNormalizedGradientFieldImageToImageMetric.h"
#include "../../../Metrics/MANGFImageToImageMetric/Code/itkNGFMetricKernel.h"
#include "../../../Metrics/MANGFImageToImageMetric/itkMANGF.h"
#include "../../../Metrics/MANGF2ImageToImageMetric/itkMANGF2.h"
#include "../../../Metrics/MANGFMSEImageToImageMetric/itkMANGFMSE.h"
#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"
#include "itkBSplineTransform.h"
#include "itkLBFGSBOptimizer.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkTransformToDeformationFieldSource.h"
#include "itkMeanSquaresImageToImageMetric.h"




#include "itkBSplineTransformInitializer.h"



#include "itkCommand.h"
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
  typedef itk::LBFGSBOptimizer         OptimizerType;
  typedef   const OptimizerType *      OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
    Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
    OptimizerPointer optimizer =
      dynamic_cast< OptimizerPointer >( object );
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
	int T=1;
		int S=2;
		int O=3;
		int OVF=4;
	 	int MMODE=5; int mmode=0;
		int BIN=MMODE+1; int NB=50;
		int MAPERCENTAGE=BIN+1; double nPCTMA=0.1;
		int LAMBDA=MAPERCENTAGE+1; double lambda=1.0;
		int LAMBDADERIVATIVE=LAMBDA+1;double lambdaderivative=1e-4;
		int ETAF=LAMBDADERIVATIVE+1;double etaF=2;
		int ETAM=ETAF+1;	double etaM=2;
		int NGFEV=ETAM+1; int NGFEvaluator=0;//scalar
		int CP=NGFEV+1; 	int numberOfLevels=1;
		int IT=CP+1; int cp=10;
		int NL=IT+1; 		int NI=200;

		int no=-1;


			int NT=2;

  if( argc < 4 )
    {
		std::cerr << "Missing Parameters " << std::endl;
			std::cerr << "Usage: " << argv[0];
			std::cerr << " \n01 fixedImageFile  \n02 movingImageFile ";
			std::cerr << " \n03 outputImagefile \n04 outputVF"
					"\n05 metric mode (0 MATTES, 1 NGF, 2 MANGF+, 3 MANG*,4 MANGFMSE,5MSE) "
					"\n06 MA number of bins \n07 MA percentage"
					"\n08 lambda term"
					"\n09 lambda derivative"
					"\n10 EtaF \n11EtaM"
					"\n12 NGF Evaluator (0 scalar,1cross,2scdelta,3Delta,4Delta2)"
					"\n13 control points \n14 iterations\n15 Levels " << no << " to skip feature\n\neros.montin@polimi.it" << std::endl;
			return EXIT_FAILURE;
    }


  if (argc>MMODE )
  		if(atoi(argv[MMODE])!=no){mmode=atof(argv[MMODE]);}


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

    if (argc>CP )
    		if(atoi(argv[CP])!=no){cp=atof(argv[CP]);}
	if (argc>IT)
			if( atoi(argv[IT])!=no){NI=atof(argv[IT]);}



  const    unsigned int    ImageDimension = 2;
  typedef  unsigned short   PixelType;

  typedef itk::Image< PixelType, ImageDimension >  ImageType;



  const unsigned int SpaceDimension = ImageDimension;
  const unsigned int SplineOrder = 3;
  typedef double CoordinateRepType;

  typedef itk::BSplineTransform<
                            CoordinateRepType,
                            SpaceDimension,
                            SplineOrder >     TransformType;


  typedef itk::LBFGSBOptimizer       OptimizerType;


  typedef itk::MANGF<
                                    ImageType,
                                    ImageType >    MetricType;


  typedef itk:: LinearInterpolateImageFunction<
                                    ImageType,
                                    double          >    InterpolatorType;

  typedef itk::ImageRegistrationMethod<
                                    ImageType,
                                    ImageType >    RegistrationType;

  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();



  registration->SetOptimizer(     optimizer     );
  registration->SetInterpolator(  interpolator  );


 
  TransformType::Pointer  transform = TransformType::New();
  registration->SetTransform( transform );
 
  typedef itk::ImageFileReader< ImageType  > FixedImageReaderType;
  typedef itk::ImageFileReader< ImageType > MovingImageReaderType;

  FixedImageReaderType::Pointer  fixedImageReader = FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader =
    MovingImageReaderType::New();

  fixedImageReader->SetFileName(  argv[T] );
  movingImageReader->SetFileName( argv[S] );

  ImageType::ConstPointer fixedImage = fixedImageReader->GetOutput();

  registration->SetFixedImage(  fixedImage   );
  registration->SetMovingImage(   movingImageReader->GetOutput()   );

  fixedImageReader->Update();

  ImageType::RegionType fixedRegion = fixedImage->GetBufferedRegion();

 registration->SetFixedImageRegion( fixedRegion );

  // Software Guide : BeginCodeSnippet

  unsigned int numberOfGridNodesInOneDimension = cp;

  TransformType::PhysicalDimensionsType   fixedPhysicalDimensions;
  TransformType::MeshSizeType             meshSize;
  TransformType::OriginType               fixedOrigin;

  for( unsigned int i=0; i< SpaceDimension; i++ )
    {
    fixedOrigin = fixedImage->GetOrigin()[i];
    fixedPhysicalDimensions[i] = fixedImage->GetSpacing()[i] *
      static_cast<double>(
      fixedImage->GetLargestPossibleRegion().GetSize()[i] - 1 );
    }
  meshSize.Fill( numberOfGridNodesInOneDimension - SplineOrder );

  transform->SetTransformDomainOrigin( fixedOrigin );
  transform->SetTransformDomainPhysicalDimensions(
    fixedPhysicalDimensions );
  transform->SetTransformDomainMeshSize( meshSize );
  transform->SetTransformDomainDirection( fixedImage->GetDirection() );

  typedef TransformType::ParametersType     ParametersType;

  const unsigned int numberOfParameters =
               transform->GetNumberOfParameters();

  ParametersType parameters( numberOfParameters );

  parameters.Fill( 0.0 );

  transform->SetParameters( parameters );

  registration->SetInitialTransformParameters( transform->GetParameters() );
  OptimizerType::BoundSelectionType boundSelect(
    transform->GetNumberOfParameters() );
  OptimizerType::BoundValueType upperBound( transform->GetNumberOfParameters() );
  OptimizerType::BoundValueType lowerBound( transform->GetNumberOfParameters() );

  boundSelect.Fill( 0 );
  upperBound.Fill( 0.0 );
  lowerBound.Fill( 0.0 );

  optimizer->SetBoundSelection( boundSelect );
  optimizer->SetUpperBound( upperBound );
  optimizer->SetLowerBound( lowerBound );

//Set/Get the CostFunctionConvergenceFactor. Algorithm terminates when the reduction in cost function is less than factor * epsmcj where epsmch is the machine precision. Typical values for factor: 1e+12 for low accuracy; 1e+7 for moderate accuracy and 1e+1 for extremely high accuracy. 
  optimizer->SetCostFunctionConvergenceFactor( 1e+7 );
//Set/Get the ProjectedGradientTolerance. Algorithm terminates when the project gradient is below the tolerance. Default value is 1e-5.
  optimizer->SetProjectedGradientTolerance( 1e-5);
  optimizer->SetMaximumNumberOfIterations( NI );
  optimizer->SetMaximumNumberOfEvaluations(NI );
  optimizer->SetMaximumNumberOfCorrections( 7 );
  // Software Guide : EndCodeSnippet

  // Create the Command observer and register it with the optimizer.
  //
  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );

 

  const unsigned int numberOfPixels = fixedImageReader->GetOutput()->GetLargestPossibleRegion().GetNumberOfPixels();
  const unsigned int numberOfSamplesMA =static_cast< unsigned int >( numberOfPixels * nPCTMA);



  if (mmode==0){
  	typedef itk::MattesMutualInformationImageToImageMetric<
  					ImageType,
  					ImageType >   MetricType;
  		MetricType::Pointer         metric        = MetricType::New();
  			registration->SetMetric( metric  );
  			metric->SetNumberOfHistogramBins(NB);
  			metric->SetNumberOfSpatialSamples(numberOfSamplesMA);
  std::cout<<"MATTES\n"<<
  							  "\n\tBin: " <<metric->GetNumberOfHistogramBins() <<
  							  "\n\tSpatial samples: " <<metric->GetNumberOfSpatialSamples()<<
  				"\n";

  }else if (mmode==1) {
  	std::cout<<"NGF\n";

  	typedef itk::NormalizedGradientFieldImageToImageMetric < ImageType , ImageType >MetricType;
  	typedef MetricType::MovingNGFType MovingNGFType;
  	typedef MetricType::FixedNGFType  FixedNGFType;
  	  MetricType::Pointer metric = MetricType::New();


  		  if (NGFEvaluator==0) {metric->SetEvaluator(new itk::NGFScalarKernel<MovingNGFType, FixedNGFType>());
  		  }else if (NGFEvaluator==1) {metric->SetEvaluator(new itk::NGFCrossKernel<MovingNGFType, FixedNGFType>());
  		  }	else if (NGFEvaluator==2){metric->SetEvaluator(new itk::NGFScaledDeltaKernel<MovingNGFType, FixedNGFType>());
  		  } else if (NGFEvaluator==3){metric->SetEvaluator(new itk::NGFDeltaKernel<MovingNGFType, FixedNGFType>());
  		  } else if (NGFEvaluator==4){metric->SetEvaluator(new itk::NGFDelta2Kernel<MovingNGFType, FixedNGFType>());}

  	registration->SetMetric( metric  );
  }else if (mmode==2) {
  		std::cout<<"MANGF+\n";
  		typedef itk::MANGF<
  						ImageType,
  						ImageType >   MetricType;
  		MetricType::Pointer         metric        = MetricType::New();


  	registration->SetMetric( metric  );
  	metric->Setlambda(lambda);
  	metric->SetlambdaDerivative(lambdaderivative);
  	metric->SetBinNumbers(NB);
  	metric->SetMAnumberOfSamples(numberOfSamplesMA);
  	metric->SetNumberOfThreads(NT);
  	metric->SetFixedEta(etaF);
  	metric->SetMovingEta(etaM);
  	metric->SetEvaluator(NGFEvaluator);


  	  std::cout<<"\n\tlambda: " << metric->Getlambda() <<
  								  "\n\tlambda derivative: " << metric->GetlambdaDerivative() <<
  								  "\n\tBin: " <<metric->GetBinNumbers() <<
  								  "\n\tETAF: " <<metric->GetFixedEta() <<
  								  "\n\t ETAM: " <<metric->GetMovingEta() <<
  								  "\n\tMA Pixels Number: "<<metric->GetMAnumberOfSamples()<<
  								  "\n\tNGF Pixels Number: "<<metric->GetNGFnumberOfSamples()<<
  								  "\n\tNGF Evaluator: " << (int) metric->GetEvaluator()<<"\n";

  }else if (mmode==3) {
  	std::cout<<"MANGF*\n";
  	typedef itk::MANGF2<
  					ImageType,
  					ImageType >   MetricType;
  	MetricType::Pointer         metric        = MetricType::New();


  registration->SetMetric( metric  );
  metric->Setlambda(lambda);
  metric->SetlambdaDerivative(lambdaderivative);
  metric->SetBinNumbers(NB);
  metric->SetMAnumberOfSamples(numberOfSamplesMA);
  metric->SetNumberOfThreads(NT);
  metric->SetFixedEta(etaF);
  metric->SetMovingEta(etaM);
  metric->SetEvaluator(NGFEvaluator);


    std::cout<<"\n\tlambda: " << metric->Getlambda() <<
  							  "\n\tlambda derivative: " << metric->GetlambdaDerivative() <<
  							  "\n\tBin: " <<metric->GetBinNumbers() <<
  							  "\n\tETAF: " <<metric->GetFixedEta() <<
  							  "\n\t ETAM: " <<metric->GetMovingEta() <<
  							  "\n\tMA Pixels Number: "<<metric->GetMAnumberOfSamples()<<
  							  "\n\tNGF Pixels Number: "<<metric->GetNGFnumberOfSamples()<<
  							  "\n\tNGF Evaluator: " << (int) metric->GetEvaluator()<<"\n";

  }else if (mmode==4) {
	std::cout<<"MA+NGF+MSE\n";
	typedef itk::MANGFMSE<
					ImageType,
					ImageType >   MetricType;
	MetricType::Pointer         metric        = MetricType::New();


double mu=0.01;
double muderivative=0.1;

registration->SetMetric( metric  );
metric->Setlambda(lambda);
metric->SetlambdaDerivative(lambdaderivative);
metric->SetBinNumbers(NB);
metric->SetMAnumberOfSamples(numberOfSamplesMA);
metric->SetNumberOfThreads(NT);
metric->SetFixedEta(etaF);
metric->SetMovingEta(etaM);
metric->SetEvaluator(NGFEvaluator);
metric->Setmu(mu);
metric->SetmuDerivative(muderivative);


  std::cout<<"\n\tlambda: " << metric->Getlambda() <<
							  "\n\tlambda derivative: " << metric->GetlambdaDerivative() <<
							  "\n\tBin: " <<metric->GetBinNumbers() <<
							  "\n\tETAF: " <<metric->GetFixedEta() <<
							  "\n\t ETAM: " <<metric->GetMovingEta() <<
							  "\n\tMA Pixels Number: "<<metric->GetMAnumberOfSamples()<<
							  "\n\tNGF Pixels Number: "<<metric->GetNGFnumberOfSamples()<<
							  "\n\tNGF Evaluator: " << (int) metric->GetEvaluator()<<"\n";

}else if (mmode==5){
typedef itk::MeanSquaresImageToImageMetric<
  					ImageType,
  					ImageType >   MetricType;
  		MetricType::Pointer         metric        = MetricType::New();
  			registration->SetMetric( metric  );
  			metric->SetNumberOfSpatialSamples(numberOfSamplesMA);
  std::cout<<"MeansSquared\n"<<
  							  "\n\tSpatial samples: " <<metric->GetNumberOfSpatialSamples()<<
  				"\n";

  }







	  // Add a time probe
  itk::TimeProbesCollectorBase chronometer;
  itk::MemoryProbesCollectorBase memorymeter;

  std::cout << std::endl << "Starting Registration" << std::endl;

  try
    {
    memorymeter.Start( "Registration" );
    chronometer.Start( "Registration" );

    registration->Update();

    chronometer.Stop( "Registration" );
    memorymeter.Stop( "Registration" );

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

  OptimizerType::ParametersType finalParameters =
                    registration->GetLastTransformParameters();


  // Report the time and memory taken by the registration
  chronometer.Report( std::cout );
  memorymeter.Report( std::cout );

  // Software Guide : BeginCodeSnippet
  transform->SetParameters( finalParameters );
  // Software Guide : EndCodeSnippet


  typedef itk::ResampleImageFilter<
                            ImageType,
                            ImageType >    ResampleFilterType;

  ResampleFilterType::Pointer resample = ResampleFilterType::New();

  resample->SetTransform( transform );
  resample->SetInput( movingImageReader->GetOutput() );

  resample->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
  resample->SetOutputOrigin(  fixedImage->GetOrigin() );
  resample->SetOutputSpacing( fixedImage->GetSpacing() );
  resample->SetOutputDirection( fixedImage->GetDirection() );
  resample->SetDefaultPixelValue( 0 );
  resample->SetNumberOfThreads(NT);

  typedef  unsigned char  OutputPixelType;

  typedef itk::Image< OutputPixelType, ImageDimension > OutputImageType;

  typedef itk::CastImageFilter<
                        ImageType,
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
		if (atoi(argv[OVF])!=-1)
	{
		typedef itk::Vector< float,  ImageDimension >  VectorType;
		typedef itk::Image< VectorType,  ImageDimension >   OutputTransformationImageType;
		typedef itk::TransformToDeformationFieldSource< OutputTransformationImageType, double >TransformToDeformationFieldSourceType;
		TransformToDeformationFieldSourceType::Pointer td = TransformToDeformationFieldSourceType::New();
		td->SetOutputParametersFromImage(fixedImageReader->GetOutput());
		td->SetTransform( registration->GetOutput()->Get() );
		//std::cout<<registration->GetOutput()->Get()<<std::endl;

		typedef itk::ImageFileWriter< OutputTransformationImageType>TransformToDeformationFieldSourceWriterType;
		TransformToDeformationFieldSourceWriterType::Pointer rtd = TransformToDeformationFieldSourceWriterType::New();
		rtd->SetInput(td->GetOutput());
		rtd->SetFileName(argv[OVF]);
		rtd->Update();
	};
 
//  // Generate the explicit deformation field resulting from
//  // the registration.
//  if( argc > 8)
//    {
//
//    typedef itk::Vector< float, ImageDimension >      VectorType;
//    typedef itk::Image< VectorType, ImageDimension >  DisplacementFieldType;
//
//    DisplacementFieldType::Pointer field = DisplacementFieldType::New();
//    field->SetRegions( fixedRegion );
//    field->SetOrigin( fixedImage->GetOrigin() );
//    field->SetSpacing( fixedImage->GetSpacing() );
//    field->SetDirection( fixedImage->GetDirection() );
//    field->Allocate();
//
//    typedef itk::ImageRegionIterator< DisplacementFieldType > FieldIterator;
//    FieldIterator fi( field, fixedRegion );
//
//    fi.GoToBegin();
//
//    TransformType::InputPointType  fixedPoint;
//    TransformType::OutputPointType movingPoint;
//    DisplacementFieldType::IndexType index;
//
//    VectorType displacement;
//
//    while( ! fi.IsAtEnd() )
//      {
//      index = fi.GetIndex();
//      field->TransformIndexToPhysicalPoint( index, fixedPoint );
//      movingPoint = transform->TransformPoint( fixedPoint );
//      displacement = movingPoint - fixedPoint;
//      fi.Set( displacement );
//      ++fi;
//      }
//
//    typedef itk::ImageFileWriter< DisplacementFieldType >  FieldWriterType;
//    FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
//
//    fieldWriter->SetInput( field );
//
//    fieldWriter->SetFileName( argv[6] );
//    try
//      {
//      fieldWriter->Update();
//      }
//    catch( itk::ExceptionObject & excp )
//      {
//      std::cerr << "Exception thrown " << std::endl;
//      std::cerr << excp << std::endl;
//      return EXIT_FAILURE;
//      }
//    }
//
//  // Optionally, save the transform parameters in a file
//  if( argc > 9 )
//    {
//    std::ofstream parametersFile;
//    parametersFile.open( argv[9] );
//    parametersFile << finalParameters << std::endl;
//    parametersFile.close();
//    }

  return EXIT_SUCCESS;
}
