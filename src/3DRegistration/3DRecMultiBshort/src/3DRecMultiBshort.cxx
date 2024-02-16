#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "../../../Metrics/MANGFImageToImageMetric/Code/itkNormalizedGradientFieldImageToImageMetric.h"
#include "../../../Metrics/MANGFImageToImageMetric/Code/itkNGFMetricKernel.h"
#include "../../../Metrics/MANGFImageToImageMetric/itkMANGF.h"
#include "../../../Metrics/MANGF2ImageToImageMetric/itkMANGF2.h"
#include "../../../Metrics/MANGFMSEImageToImageMetric/itkMANGFMSE.h"
#include "itkCastImageFilter.h"
#include "itkBSplineTransform.h"
#include "itkResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkLBFGSBOptimizer.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"


#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkImageMaskSpatialObject.h"
#include "itkTransformToDeformationFieldSource.h"

#include "itkVersion.h"
#include "itkTransformFileWriter.h"
#include "itkTransformFileReader.h"

#include <iostream>
#include <algorithm> // for copy
#include <iterator> // for ostream_iterator
#include <vector>
#include<fstream>
#include<sstream>
#include<iostream>
#include "../../Version.h"



#include "itkVersion.h"
#include "itkTransformFileWriter.h"
#include "itkTransformFileReader.h"




const    unsigned int    Dimension = 3;


typedef  short  PixelType;


typedef itk::Image< PixelType, Dimension >  FixedImageType;
typedef itk::Image< PixelType, Dimension >  MovingImageType;
typedef   float                                    InternalPixelType;
typedef itk::Image< InternalPixelType, Dimension > InternalImageType;


typedef itk::LinearInterpolateImageFunction<
		InternalImageType,
		double             > InterpolatorType;

//read
typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
typedef itk::CastImageFilter<
		FixedImageType, InternalImageType > FixedCastFilterType;
typedef itk::CastImageFilter<
		MovingImageType, InternalImageType > MovingCastFilterType;



typedef itk::MultiResolutionImageRegistrationMethod<
		InternalImageType,
		InternalImageType >   RegistrationType;

typedef itk::MultiResolutionPyramidImageFilter<
		InternalImageType,
		InternalImageType >   FixedImagePyramidType;
typedef itk::MultiResolutionPyramidImageFilter<
		InternalImageType,
		InternalImageType >   MovingImagePyramidType;

const unsigned int SpaceDimension = Dimension;
const unsigned int SplineOrder = 3;
typedef double CoordinateRepType;

#if ITK_VERSION_MAJOR < 4
typedef itk::BSplineDeformableTransform<
		CoordinateRepType,
		SpaceDimension,
		SplineOrder >     TransformType;
#else
typedef itk::BSplineTransform<
		CoordinateRepType,
		SpaceDimension,
		SplineOrder >     TransformType;
#endif
typedef itk::LBFGSBOptimizer       OptimizerType;

//output
typedef itk::ResampleImageFilter<
		InternalImageType,
		FixedImageType >    ResampleFilterType;

typedef unsigned char OutputPixelType;
typedef itk::Image< PixelType, Dimension >  OutputImageType;

typedef itk::CastImageFilter<FixedImageType,
		OutputImageType > CastFilterType;

typedef itk::ImageFileWriter< OutputImageType >  WriterType;


typedef itk::Vector< float,  Dimension >  VectorType;
typedef itk::Image< VectorType,  Dimension >   OutputTransformationImageType;

typedef RegistrationType::ParametersType ParametersType;

#include "../../../include/funzioniIO.h"
#include "../../../include/funzioniREC.h"
#include "../../../include/funzioniIOtranformationTXT.h"


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



	std::vector<std::string>Vparams;
	Vparams.push_back("program name");
	int T=1;Vparams.push_back("fixedImageFile");
	int S=2;Vparams.push_back("movingImageFile");
	int O=3;Vparams.push_back("output Image");
	int OVF=4;Vparams.push_back("VF out");
	int MMODE=5; int mmode=0;Vparams.push_back(	availableMetrics()); readParametersFromStdin(&mmode,MMODE,argv,argc,no);
	int BIN=MMODE+1; int NB=50;Vparams.push_back("Mattes number of bins");readParametersFromStdin(&NB,BIN,argv,argc,no);
	int MAPERCENTAGE=BIN+1; double nPCTMA=0.1; Vparams.push_back("Mattes percecntage");readParametersFromStdin(&nPCTMA,MAPERCENTAGE,argv,argc,no);
	int LAMBDA=MAPERCENTAGE+1; double lambda=1.0;Vparams.push_back("lambda values MI + lambda NGF");readParametersFromStdin(&lambda,LAMBDA,argv,argc,no);
	int LAMBDADERIVATIVE=LAMBDA+1;double lambdaderivative=1e-4;Vparams.push_back("Lambda derivatuive ");readParametersFromStdin(&lambdaderivative,LAMBDADERIVATIVE,argv,argc,no);
	int ETAF=LAMBDADERIVATIVE+1;double etaF=2; Vparams.push_back("Eta values Reference (NGF noise)");readParametersFromStdin(&etaF,ETAF,argv,argc,no);
	int ETAM=ETAF+1;	double etaM=2; Vparams.push_back("Eta values Source (NGF noise)");readParametersFromStdin(&etaM,ETAM,argv,argc,no);
	int NGFEV=ETAM+1; int NGFEvaluator=0;Vparams.push_back("NGF Evaluator (0 scalar,1cross,2scdelta,3Delta,4Delta2)");readParametersFromStdin(&NGFEvaluator,NGFEV,argv,argc,no);
	int MU=NGFEV+1; double mu=0.1;Vparams.push_back("Tau MSE"); readParametersFromStdin(&mu,MU,argv,argc,no);
	int MUDERIVATIVES=MU+1; double muderivatives=0.1;Vparams.push_back("Tau Derivative MSE ");readParametersFromStdin(&muderivatives,MUDERIVATIVES,argv,argc,no);
	int NL=MUDERIVATIVES+1; 	int numberOfLevels=1; Vparams.push_back("Levels");readParametersFromStdin(&numberOfLevels,NL,argv,argc,no);
	int GRID=NL+1; double numberOfGridNodes=8;Vparams.push_back("Grid nodes: " + tostr(numberOfGridNodes) );	readParametersFromStdin(&numberOfGridNodes,GRID,argv,argc,no);
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
	int TIN=HBOUND+1; Vparams.push_back("transformin.txt");
	int TOUT=TIN+1; Vparams.push_back("transforout.txt");
	int NT=TOUT; int nt=2; Vparams.push_back("Number of thread" + tostr(nt));readParametersFromStdin(&nt,NT,argv,argc,no);
	int FMASK=NT+1;Vparams.push_back("Mask Reference");
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
	checkIn(Vparams,argv,argc);
	TransformType::Pointer      transform     = TransformType::New();
	OptimizerType::Pointer      optimizer     = OptimizerType::New();
	InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
	RegistrationType::Pointer   registration  = RegistrationType::New();

	FixedImagePyramidType::Pointer fixedImagePyramid =
			FixedImagePyramidType::New();
	MovingImagePyramidType::Pointer movingImagePyramid =
			MovingImagePyramidType::New();

	registration->SetOptimizer(     optimizer     );
	registration->SetTransform(     transform     );
	registration->SetInterpolator(  interpolator  );

	registration->SetFixedImagePyramid( fixedImagePyramid );
	registration->SetMovingImagePyramid( movingImagePyramid );
	registration->SetNumberOfThreads(nt);

	fixedImagePyramid->SetNumberOfThreads(nt);
	movingImagePyramid->SetNumberOfThreads(nt);


	FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
	MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

	fixedImageReader->SetFileName(  argv[T] );
	movingImageReader->SetFileName( argv[S] );




	FixedCastFilterType::Pointer fixedCaster   = FixedCastFilterType::New();
	MovingCastFilterType::Pointer movingCaster = MovingCastFilterType::New();
	fixedCaster->SetNumberOfThreads(nt);
	movingCaster->SetNumberOfThreads(nt);

	fixedCaster->SetInput(  fixedImageReader->GetOutput() );
	movingCaster->SetInput( movingImageReader->GetOutput() );

	registration->SetFixedImage(    fixedCaster->GetOutput()    );
	registration->SetMovingImage(   movingCaster->GetOutput()   );


	fixedCaster->Update();

	registration->SetFixedImageRegion(
			fixedCaster->GetOutput()->GetBufferedRegion() );


	//  The reader should note that the BSpline computation requires a
		//  finite support region ( 1 grid node at the lower borders and 2
		//  grid nodes at upper borders). place the grid origin such that
		//  grid node (1,1) coincides with the first pixel in the fixed image.


		TransformType::PhysicalDimensionsType   fixedPhysicalDimensions;
		TransformType::MeshSizeType             meshSize;
		for( unsigned int i=0; i < Dimension; i++ )
		{
			fixedPhysicalDimensions[i] = registration->GetFixedImage()->GetSpacing()[i] *
					static_cast<double>(
							registration->GetFixedImage()->GetLargestPossibleRegion().GetSize()[i] - 1 );
		}
		unsigned int numberOfGridNodesInOneDimension = numberOfGridNodes;
		meshSize.Fill( numberOfGridNodesInOneDimension - SplineOrder );
		transform->SetTransformDomainOrigin( registration->GetFixedImage()->GetOrigin() );
		transform->SetTransformDomainPhysicalDimensions( fixedPhysicalDimensions );
		transform->SetTransformDomainMeshSize( meshSize );
		transform->SetTransformDomainDirection( registration->GetFixedImage()->GetDirection() );

		typedef TransformType::ParametersType     ParametersType;



	typedef RegistrationType::ParametersType ParametersType;
	ParametersType initialParameters( transform->GetNumberOfParameters() );



	//intializeit!!
	if (argc>TIN)
		if(*argv[TIN]!=no){

#if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 5) || ITK_VERSION_MAJOR > 4
			itk::TransformFileReaderTemplate<double>::Pointer reader =
					itk::TransformFileReaderTemplate<double>::New();
#else
			itk::TransformFileReader::Pointer reader = itk::TransformFileReader::New();
#endif


			reader->SetFileName(argv[TIN]);
			reader->Update();
			initialParameters= reader->GetTransformList()->begin()->GetPointer()->GetParameters();
			std::cout<<"intial transform: " << initialParameters <<"\n";
		} else {
			transform->SetIdentity();

			initialParameters= transform->GetParameters();

		};

	registration->SetInitialTransformParameters( initialParameters );


	//with this i deete everything before
	typedef TransformType::ParametersType     ParametersType;

		const unsigned int numberOfParameters =
				transform->GetNumberOfParameters();

	ParametersType parametersLow( numberOfParameters );

		parametersLow.Fill( 0.0 );

		transform->SetParameters( parametersLow );


		registration->SetInitialTransformParameters( transform->GetParameters() );


		//  Next we set the parameters of the LBFGS Optimizer.


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


	optimizer->SetMinimize(true);






	const unsigned int numberOfPixels = fixedImageReader->GetOutput()->GetLargestPossibleRegion().GetNumberOfPixels();
	const unsigned int numberOfSamplesMA =static_cast< unsigned int >( numberOfPixels * nPCTMA);
	///////

	//parameter metric to print
	std::vector<std::string> PMETRIC;
	std::vector<std::string> VMETRIC;


	qualemetrica(registration, mmode,NB, numberOfSamplesMA,
			lambda, lambdaderivative, etaF,
			etaM, NGFEvaluator, mu, muderivatives,
			nt,PMETRIC,VMETRIC);

	pMetric("metric",PMETRIC,VMETRIC);
	//MASCHERE
	maschere(registration,FMASK,MMASK, argc, argv,no);


	std::vector<std::string> PR;
	std::vector<std::string> VR;

	registration->SetNumberOfLevels(numberOfLevels);PR.push_back("Registration number of levels ");VR.push_back(tostr(registration->GetNumberOfLevels()));
	pMetric("registration",PR,VR);


	std::cout << std::endl << "Starting Registration" << std::endl;

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



	OptimizerType::ParametersType finalParameters =
					registration->GetLastTransformParameters();



			TransformType::Pointer finalTransform = TransformType::New();

			finalTransform->SetParameters( finalParameters );
			finalTransform->SetFixedParameters( transform->GetFixedParameters() );

			ResampleFilterType::Pointer resample = ResampleFilterType::New();
			WriterType::Pointer      writer =  WriterType::New();
			CastFilterType::Pointer  caster =  CastFilterType::New();


			transformedImageWrite(fixedCaster->GetOutput(),movingCaster->GetOutput(),finalTransform,resample,writer,caster,O,nt,argc,argv);

			if (argc>OVF)
					if (*argv[OVF]!=no)
					{
			transformationWrite(registration,OVF,argc,argv,no);
					};



//			if (argc>TOUT)
//				if (*argv[TOUT]!=no)
//				{
//					writeTranformationTxt(finalTransform, TOUT, argv);
//				};

	return EXIT_SUCCESS;
}







