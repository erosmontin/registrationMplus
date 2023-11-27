

#include "itkMultiResolutionImageRegistrationMethod.h"
//#include "itkCenteredEuler3DTransform.h"
#include "itkRigid3DTransform.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkImageMaskSpatialObject.h"
#include "itkTransformToDeformationFieldSource.h"

#include "itkVersion.h"
#include "itkTransformFileWriter.h"
#include "itkTransformFileReader.h"


#include "Version.h"


const    unsigned int    Dimension = 3;
typedef  unsigned char  PixelType;
typedef  unsigned char  OutputPixelType;

typedef itk::Image< PixelType, Dimension >  FixedImageType;
typedef itk::Image< PixelType, Dimension >  MovingImageType;


typedef   float                                    InternalPixelType;
typedef itk::Image< InternalPixelType, Dimension > InternalImageType;


typedef itk::CastImageFilter<
		FixedImageType, InternalImageType > FixedCastFilterType;
typedef itk::CastImageFilter<
		MovingImageType, InternalImageType > MovingCastFilterType;

//typedef itk::CenteredEuler3DTransform< double > TransformType;
typedef itk::Rigid3DTransform< double > TransformType;
typedef itk::RegularStepGradientDescentOptimizer       OptimizerType;
typedef itk::LinearInterpolateImageFunction<
		InternalImageType,
		double             > InterpolatorType;

typedef itk::MultiResolutionImageRegistrationMethod<
		InternalImageType,
		InternalImageType >   RegistrationType;

typedef itk::MultiResolutionPyramidImageFilter<
		InternalImageType,
		InternalImageType >   FixedImagePyramidType;
typedef itk::MultiResolutionPyramidImageFilter<
		InternalImageType,
		InternalImageType >   MovingImagePyramidType;

typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;

typedef RegistrationType::ParametersType ParametersType;






//OUPTUT
typedef itk::ResampleImageFilter<InternalImageType,InternalImageType >    ResampleFilterType;

typedef itk::Image< OutputPixelType, Dimension > OutputImageType;

	typedef itk::CastImageFilter<
			InternalImageType,
			OutputImageType > CastFilterType;

	typedef itk::ImageFileWriter< OutputImageType >  WriterType;

	typedef itk::Vector< float,  Dimension >  VectorType;
					typedef itk::Image< VectorType,  Dimension >   OutputTransformationImageType;


#include "../../../include/funzioniIO.h"	
#include "../../../include/funzioniREC.h"
float maxsteplengthdivisor =4;
float minimumsteplengthdivisor =10;
int nidivisor =1;
#include "../../../include/commandRSG.h"


#include "itkArray.h"

#include "../../../include/funzioniIOtranformationTXT.h"








int main( int argc, char *argv[] )
{


	std::vector<std::string>Vparams;
	Vparams.push_back("program name");
	int T=1;Vparams.push_back("fixedImageFile");
	int S=2;Vparams.push_back("movingImageFile");
	int O=3;Vparams.push_back("output Image");
	int OVF=4;Vparams.push_back("VF out");
	int MMODE=5; int mmode=0;Vparams.push_back("metric mode (0 MATTES, 1 NGF, 2 MANGF+, 3 MANG*, 4 MANGFMSE, MSE)"); readParametersFromStdin(&mmode,MMODE,argv,argc,no);
	int MAXSTEPLENGTH=MMODE+1; double maximumStepLength=.3;Vparams.push_back("MAX RSG"); readParametersFromStdin(&maximumStepLength,MAXSTEPLENGTH,argv,argc,no);
	int RSGMAXDIVISOR=MAXSTEPLENGTH+1; Vparams.push_back("MAX RSG divisor"); readParametersFromStdin(&maxsteplengthdivisor,RSGMAXDIVISOR,argv,argc,no);
	int MINSTEPLENGTH=RSGMAXDIVISOR+1; 	double minimumStepLength=.001; Vparams.push_back("MIN RSG ");readParametersFromStdin(&minimumStepLength,MINSTEPLENGTH,argv,argc,no);
	int RSGMINDIVISOR=MINSTEPLENGTH+1; Vparams.push_back("Min RSG divisor");readParametersFromStdin(&minimumsteplengthdivisor,RSGMINDIVISOR,argv,argc,no);
	int BIN=RSGMINDIVISOR+1; int NB=50;Vparams.push_back("Mattes number of bins");readParametersFromStdin(&NB,BIN,argv,argc,no);
	int MAPERCENTAGE=BIN+1; double nPCTMA=0.1; Vparams.push_back("Mattes percecntage");readParametersFromStdin(&nPCTMA,MAPERCENTAGE,argv,argc,no);
	int LAMBDA=MAPERCENTAGE+1; double lambda=1.0;Vparams.push_back("lambda values MI + lambda NGF");readParametersFromStdin(&lambda,LAMBDA,argv,argc,no);
	int LAMBDADERIVATIVE=LAMBDA+1;double lambdaderivative=1e-4;Vparams.push_back("Lambda derivatuive ");readParametersFromStdin(&lambdaderivative,LAMBDADERIVATIVE,argv,argc,no);
	int ETAF=LAMBDADERIVATIVE+1;double etaF=2; Vparams.push_back("Eta values Reference (NGF noise)");readParametersFromStdin(&etaF,ETAF,argv,argc,no);
	int ETAM=ETAF+1;	double etaM=2; Vparams.push_back("Eta values Source (NGF noise)");readParametersFromStdin(&etaM,ETAM,argv,argc,no);
	int NGFEV=ETAM+1; int NGFEvaluator=0;Vparams.push_back("NGF Evaluator (0 scalar,1cross,2scdelta,3Delta,4Delta2)");	readParametersFromStdin(&NGFEvaluator,NGFEV,argv,argc,no);
	int MU=NGFEV+1; double mu=0.1;Vparams.push_back("mu MSE(0.1)");	readParametersFromStdin(&mu,MU,argv,argc,no);
	int MUDERIVATIVES=MU+1; double muderivatives=0.1;Vparams.push_back("Mu Derivative MSE ");	readParametersFromStdin(&muderivatives,MUDERIVATIVES,argv,argc,no);
	int NL=MUDERIVATIVES+1; 	int numberOfLevels=1; Vparams.push_back("Levels");readParametersFromStdin(&numberOfLevels,NL,argv,argc,no);
	int NI=NL+1; 	int ni=200; Vparams.push_back("Max number of Iterations"); 	readParametersFromStdin(&ni,NI,argv,argc,no);
	int NIDIV=NI+1; 	Vparams.push_back("Max number of Iterations divisor for step (1)"); 	readParametersFromStdin(&nidivisor,NIDIV,argv,argc,no);
	int FMASK=NIDIV+1;Vparams.push_back("Mask Reference");
	int MMASK=FMASK+1; Vparams.push_back("Mask source");
	int RLF = MMASK+1; double rlf=0.6; Vparams.push_back("RLF (0.6)"); readParametersFromStdin(&rlf,RLF,argv,argc,no);
	int OPTGRADIENTTOL=RLF+1; double optgradienttol=1e-5; Vparams.push_back("optimizer mag grad tol (1e-5)");readParametersFromStdin(&optgradienttol,OPTGRADIENTTOL,argv,argc,no);
	int NT=OPTGRADIENTTOL+1; 	int nt=2; Vparams.push_back("Number of thread");readParametersFromStdin(&nt,NT,argv,argc,no);
	int TIN=NT+1; Vparams.push_back("transformin.txt");
	int TOUT=TIN+1; Vparams.push_back("transforout.txt");
	int OPTSCALE =TOUT+1; double optimizerScale[12]={1.0, 1.0 , 1.0 ,
													1.0 , 1.0, 1.0,
													1.0, 1.0, 1.0, 1.0/1000,1.0/1000,1.0/1000};
	for (int a=0;a<9;a++)
	{
		Vparams.push_back("oprimizers scales a" + tostr(a));
	}

	for (int a=0;a<3;a++)
		{
			Vparams.push_back("oprimizers scales t" + tostr(a));
		}


	for (int a=0;a<12;a++)
	{

		readParametersFromStdin(&optimizerScale[a],OPTSCALE+a,argv,argc,no);

	}



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

	//checkIn(Vparams,argv, argc);


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
	registration->SetNumberOfThreads(NT);

	fixedImagePyramid->SetNumberOfThreads(NT);
	movingImagePyramid->SetNumberOfThreads(NT);


	FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
	MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

	fixedImageReader->SetFileName(  argv[T] );
	movingImageReader->SetFileName( argv[S] );




	FixedCastFilterType::Pointer fixedCaster   = FixedCastFilterType::New();
	MovingCastFilterType::Pointer movingCaster = MovingCastFilterType::New();
	fixedCaster->SetNumberOfThreads(NT);
	movingCaster->SetNumberOfThreads(NT);

	fixedCaster->SetInput(  fixedImageReader->GetOutput() );
	movingCaster->SetInput( movingImageReader->GetOutput() );

	registration->SetFixedImage(    fixedCaster->GetOutput()    );
	registration->SetMovingImage(   movingCaster->GetOutput()   );


	fixedCaster->Update();

	registration->SetFixedImageRegion(
			fixedCaster->GetOutput()->GetBufferedRegion() );



	ParametersType initialParameters( transform->GetNumberOfParameters() );

	//intializeit!!

	if (argc>TIN)
		if( *argv[TIN]!=no){
			initialParameters=readInitFromTxt(registration,initialParameters,TIN,argv);
		} else {
			//Transformation
	typedef FixedImageType::SpacingType    SpacingType;
	  typedef FixedImageType::PointType      OriginType;
	  typedef FixedImageType::RegionType     RegionType;
	  typedef FixedImageType::SizeType       SizeType;


	  const SpacingType fixedSpacing = fixedImageReader->GetOutput()->GetSpacing();
	  const OriginType  fixedOrigin  = fixedImageReader->GetOutput()->GetOrigin();
	  const RegionType  fixedRegion  = fixedImageReader->GetOutput()->GetLargestPossibleRegion();
	  const SizeType    fixedSize    = fixedRegion.GetSize();

	  TransformType::InputPointType centerFixed;

	  centerFixed[0] = fixedOrigin[0] + fixedSpacing[0] * fixedSize[0] / 2.0;
	  centerFixed[1] = fixedOrigin[1] + fixedSpacing[1] * fixedSize[1] / 2.0;
	  centerFixed[2] = fixedOrigin[2] + fixedSpacing[2] * fixedSize[2] / 2.0;


	  const SpacingType movingSpacing = movingImageReader->GetOutput()->GetSpacing();
	  const OriginType  movingOrigin  = movingImageReader->GetOutput()->GetOrigin();
	  const RegionType  movingRegion  = movingImageReader->GetOutput()->GetLargestPossibleRegion();
	  const SizeType    movingSize    = movingRegion.GetSize();

	  TransformType::InputPointType centerMoving;

	  centerMoving[0] = movingOrigin[0] + movingSpacing[0] * movingSize[0] / 2.0;
	  centerMoving[1] = movingOrigin[1] + movingSpacing[1] * movingSize[1] / 2.0;
	  centerMoving[2] = movingOrigin[2] + movingSpacing[2] * movingSize[2] / 2.0;

	  transform->SetCenter( centerFixed );
	  transform->SetTranslation( centerMoving - centerFixed );
	  


	  //transform->SetRotation(toRad(0.0),toRad(0.0),toRad(0.0));

	  transform->SetIdentity();

	  //transform->SetMatrix(tor1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,0.0);
			
			
		};
	

	transform->SetParameters(initialParameters);
	registration->SetInitialTransformParameters( transform->GetParameters() );

	std::cout<<registration->GetInitialTransformParameters()<<'\n';



	const unsigned int numberOfPixels = fixedImageReader->GetOutput()->GetLargestPossibleRegion().GetNumberOfPixels();
	const unsigned int numberOfSamplesMA =static_cast< unsigned int >( numberOfPixels * nPCTMA);

	//parameter metric to print
	std::vector<std::string> PMETRIC;
	std::vector<std::string> VMETRIC;


	//METRICS
	qualemetrica(registration, mmode,NB, numberOfSamplesMA,
			lambda, lambdaderivative, etaF,
			etaM, NGFEvaluator, mu, muderivatives,
			nt,PMETRIC,VMETRIC);

	pMetric("metric",PMETRIC,VMETRIC);
		//MASCHERE
	maschere(registration,FMASK,MMASK, argc, argv,no);




	std::vector<std::string> POPT;
	std::vector<std::string> VOPT;

	optimizer->SetNumberOfIterations( ni ); POPT.push_back("NI: ");VOPT.push_back(tostr(optimizer->GetNumberOfIterations()));
	optimizer->MinimizeOn();
	optimizer->SetMaximumStepLength( maximumStepLength );POPT.push_back("MAX step length: ");VOPT.push_back(tostr(optimizer->GetMaximumStepLength()));
	optimizer->SetMinimumStepLength( minimumStepLength );POPT.push_back("MIN step length: ");VOPT.push_back(tostr(optimizer->GetMinimumStepLength()));
	optimizer->SetGradientMagnitudeTolerance(optgradienttol);POPT.push_back("Grad Mag tol: ");VOPT.push_back(tostr(optimizer->GetGradientMagnitudeTolerance()));
	optimizer->SetRelaxationFactor(rlf); POPT.push_back("RLF: ");VOPT.push_back(tostr(optimizer->GetRelaxationFactor()));
	optimizer->SetMinimize(1);



	typedef OptimizerType::ScalesType       OptimizerScalesType;
	OptimizerScalesType   scale = OptimizerScalesType( transform->GetNumberOfParameters() );
	int u;
	std::cout<<scale<<std::endl;
	for(u=0;u<transform->GetNumberOfParameters();u++)
	{
		scale.SetElement(u,optimizerScale[u]);
	}

	optimizer->SetScales(scale); POPT.push_back("OPT scales: ");VOPT.push_back(tostr(optimizer->GetScales()));

	pMetric("optimizer",POPT,VOPT);







	// Create the Command observer and register it with the optimizer.
	//
	CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
	optimizer->AddObserver( itk::IterationEvent(), observer );

	typedef RegistrationInterfaceCommand<RegistrationType> CommandType;
	CommandType::Pointer command = CommandType::New();
	registration->AddObserver( itk::IterationEvent(), command );


	std::vector<std::string> PR;
	std::vector<std::string> VR;

	registration->SetNumberOfLevels(numberOfLevels);PR.push_back("Registration number of levels ");VR.push_back(tostr(registration->GetNumberOfLevels()));

	pMetric("registration",PR,VR);


	try
	{
		registration->Update();
		std::cout << "Optimizer stop condition: "
				<< registration->GetOptimizer()->GetStopConditionDescription()
				<< std::endl;
	}
	catch( itk::ExceptionObject & err )
	{
		std::cout << "ExceptionObject caught !" << std::endl;
		std::cout << err << std::endl;
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



		if (argc>TOUT)
			if (*argv[TOUT]!=no)
			{
				writeTranformationTxt(finalTransform, TOUT, argv);
			};




	return EXIT_SUCCESS;
}

