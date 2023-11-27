#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkAffineTransform.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "../../../Metrics/MANGFImageToImageMetric/Code/itkNormalizedGradientFieldImageToImageMetric.h"
#include "../../../Metrics/MANGFImageToImageMetric/Code/itkNGFMetricKernel.h"
#include "../../../Metrics/MANGFImageToImageMetric/itkMANGF.h"
#include "../../../Metrics/MANGF2ImageToImageMetric/itkMANGF2.h"
#include "../../../Metrics/MANGFMSEImageToImageMetric/itkMANGFMSE.h"
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

#include <iostream>
#include <algorithm> // for copy
#include <iterator> // for ostream_iterator
#include <vector>
#include<fstream>
#include<sstream>
#include<iostream>
#include "Version.h"
#include "../../../include/old_funzioni.h"

float maxsteplengthdivisor =4;
float minimumsteplengthdivisor =10;
int nidivisor =1;
#include "../../../include/commandRSG.h"


#include "itkArray.h"

#include "itkCommand.h"







int main( int argc, char *argv[] )
{
	char no='N';

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
	int LAMBDADERIVATIVE=LAMBDA+1;double lambdaderivative=1e-4;Vparams.push_back("Lambda derivatuive ");
	int ETAF=LAMBDADERIVATIVE+1;double etaF=2; Vparams.push_back("Eta values Reference (NGF noise)");
	int ETAM=ETAF+1;	double etaM=2; Vparams.push_back("Eta values Source (NGF noise)");
	int NGFEV=ETAM+1; int NGFEvaluator=0;Vparams.push_back("NGF Evaluator (0 scalar,1cross,2scdelta,3Delta,4Delta2)");
	int MU=NGFEV+1; double mu=0.1;Vparams.push_back("Tau MSE");
	int MUDERIVATIVES=MU+1; double muderivatives=0.1;Vparams.push_back("Tau Derivative MSE ");
	int NL=MUDERIVATIVES+1; 	int numberOfLevels=1; Vparams.push_back("Levels");
	int NI=NL+1; 	int ni=200; Vparams.push_back("Max number of Iterations");
	int NIDIV=NI+1; 	Vparams.push_back("Max number of Iterations divisor for step (1)");
	int FMASK=NIDIV+1;Vparams.push_back("Mask Reference");
	int MMASK=FMASK+1; Vparams.push_back("Mask source");
	int RLF = MMASK+1; double rlf=0.6; Vparams.push_back("RLF (0.6)");
	int OPTGRADIENTTOL=RLF+1; double optgradienttol=1e-5; Vparams.push_back("optimizer mag grad tol (1e-5)");
	int NT=OPTGRADIENTTOL+1; 	int nt=2; Vparams.push_back("Number of thread");
	int TIN=NT+1; Vparams.push_back("transformin.txt");
	int TOUT=TIN+1; Vparams.push_back("transforout.txt");
	int OPTSCALE =TOUT+1; double optimizerScale[6]={1.0, 1.0, 1.0, 1.0,0.001,0.001}; Vparams.push_back("oprimizers scales a1");
	Vparams.push_back("oprimizers scales a2");
	Vparams.push_back("oprimizers scales a3");
	Vparams.push_back("oprimizers scales a4");
	Vparams.push_back("oprimizers scales t1");
	Vparams.push_back("oprimizers scales t2");








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

	readParametersFromStdin(&mmode,MMODE,argv,argc,no);
	readParametersFromStdin(&maximumStepLength,MAXSTEPLENGTH,argv,argc,no);
	readParametersFromStdin(&maxsteplengthdivisor,RSGMAXDIVISOR,argv,argc,no);
	readParametersFromStdin(&minimumStepLength,MINSTEPLENGTH,argv,argc,no);
	readParametersFromStdin(&minimumsteplengthdivisor,RSGMINDIVISOR,argv,argc,no);
	readParametersFromStdin(&NB,BIN,argv,argc,no);
	readParametersFromStdin(&nPCTMA,MAPERCENTAGE,argv,argc,no);
	readParametersFromStdin(&lambda,LAMBDA,argv,argc,no);
	readParametersFromStdin(&lambdaderivative,LAMBDADERIVATIVE,argv,argc,no);
	readParametersFromStdin(&etaF,ETAF,argv,argc,no);
	readParametersFromStdin(&etaM,ETAM,argv,argc,no);
	readParametersFromStdin(&NGFEvaluator,NGFEV,argv,argc,no);
	readParametersFromStdin(&mu,MU,argv,argc,no);
	readParametersFromStdin(&muderivatives,MUDERIVATIVES,argv,argc,no);
	readParametersFromStdin(&numberOfLevels,NL,argv,argc,no);
	readParametersFromStdin(&ni,NI,argv,argc,no);
	readParametersFromStdin(&nidivisor,NIDIV,argv,argc,no);
	readParametersFromStdin(&rlf,RLF,argv,argc,no);
	readParametersFromStdin(&optgradienttol,OPTGRADIENTTOL,argv,argc,no);
	readParametersFromStdin(&nt,NT,argv,argc,no);
	for (int a=0;a<6;a++)
	{

		readParametersFromStdin(&optimizerScale[a],OPTSCALE+a,argv,argc,no);

	}


	const    unsigned int    Dimension = 2;
	typedef  unsigned short  PixelType;
	typedef  unsigned short  OutputPixelType;

	typedef itk::Image< PixelType, Dimension >  FixedImageType;
	typedef itk::Image< PixelType, Dimension >  MovingImageType;
	typedef   float                                    InternalPixelType;
	typedef itk::Image< InternalPixelType, Dimension > InternalImageType;
	typedef itk::AffineTransform< double, Dimension > TransformType;
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
	typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
	typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;

	FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
	MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

	fixedImageReader->SetFileName(  argv[T] );
	movingImageReader->SetFileName( argv[S] );


	typedef itk::CastImageFilter<
			FixedImageType, InternalImageType > FixedCastFilterType;
	typedef itk::CastImageFilter<
			MovingImageType, InternalImageType > MovingCastFilterType;

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


	const unsigned int numberOfPixels = fixedImageReader->GetOutput()->GetLargestPossibleRegion().GetNumberOfPixels();
	const unsigned int numberOfSamplesMA =static_cast< unsigned int >( numberOfPixels * nPCTMA);

	//parameter metric to print
	std::vector<std::string> PMETRIC;
	std::vector<std::string> VMETRIC;

	if (mmode==0){
		typedef itk::MattesMutualInformationImageToImageMetric<
				InternalImageType,
				InternalImageType >   MetricType;
		MetricType::Pointer         metric        = MetricType::New();
		registration->SetMetric( metric  );
		metric->SetNumberOfHistogramBins(NB);
		metric->SetNumberOfSpatialSamples(numberOfSamplesMA);
		PMETRIC.push_back("Metric: ");VMETRIC.push_back("MI");
		PMETRIC.push_back("Number of histogram bin"); VMETRIC.push_back(tostr(metric->GetNumberOfHistogramBins()));
		PMETRIC.push_back("Spatial samples");VMETRIC.push_back(tostr(metric->GetNumberOfSpatialSamples()));
		pMetric("metric",PMETRIC,VMETRIC);
		registration->SetMetric( metric  );
		//MASK
		if(argc>FMASK)
		if(*argv[FMASK]!=no)
		{
			typedef itk::ImageMaskSpatialObject< Dimension > MaskType;
			MaskType::Pointer spatialObjectMask = MaskType::New();
			typedef itk::Image< unsigned char, Dimension > ImageMaskType;
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
		if(*argv[MMASK]!=no)
		{
			typedef itk::ImageMaskSpatialObject< Dimension > MaskType;
			MaskType::Pointer spatialObjectMask = MaskType::New();
			typedef itk::Image< unsigned char, Dimension > ImageMaskType;
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
		
		
	}else if (mmode==1) {

		typedef itk::NormalizedGradientFieldImageToImageMetric < InternalImageType , InternalImageType >MetricType;
		typedef MetricType::MovingNGFType MovingNGFType;
		typedef MetricType::FixedNGFType  FixedNGFType;
		MetricType::Pointer metric = MetricType::New();

		if (NGFEvaluator==0) {metric->SetEvaluator(new itk::NGFScalarKernel<MovingNGFType, FixedNGFType>());
		}else if (NGFEvaluator==1) {metric->SetEvaluator(new itk::NGFCrossKernel<MovingNGFType, FixedNGFType>());
		}	else if (NGFEvaluator==2){metric->SetEvaluator(new itk::NGFScaledDeltaKernel<MovingNGFType, FixedNGFType>());
		} else if (NGFEvaluator==3){metric->SetEvaluator(new itk::NGFDeltaKernel<MovingNGFType, FixedNGFType>());
		} else if (NGFEvaluator==4){metric->SetEvaluator(new itk::NGFDelta2Kernel<MovingNGFType, FixedNGFType>());
		}else {metric->SetEvaluator(new itk::NGFScalarKernel<MovingNGFType, FixedNGFType>());}

		PMETRIC.push_back("Metric: ");VMETRIC.push_back("NGF");
		PMETRIC.push_back("EValuator: ");VMETRIC.push_back(tostr(NGFEvaluator));
		pMetric("metric",PMETRIC,VMETRIC);
		registration->SetMetric( metric  );
//MASK
		if(argc>FMASK)
		if(*argv[FMASK]!=no)
		{
			typedef itk::ImageMaskSpatialObject< Dimension > MaskType;
			MaskType::Pointer spatialObjectMask = MaskType::New();
			typedef itk::Image< unsigned char, Dimension > ImageMaskType;
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
		if(*argv[MMASK]!=no)
		{
			typedef itk::ImageMaskSpatialObject< Dimension > MaskType;
			MaskType::Pointer spatialObjectMask = MaskType::New();
			typedef itk::Image< unsigned char, Dimension > ImageMaskType;
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

	}else if (mmode==2) {
		typedef itk::MANGF<
				InternalImageType,
				InternalImageType >   MetricType;
		MetricType::Pointer         metric        = MetricType::New();


		PMETRIC.push_back("Metric: ");VMETRIC.push_back("MI + lambda NGF");

		metric->Setlambda(lambda);	PMETRIC.push_back("Lambda: ");VMETRIC.push_back(tostr(metric->Getlambda()));
		metric->SetlambdaDerivative(lambdaderivative); PMETRIC.push_back("LambdaDerivative: ");VMETRIC.push_back(tostr(metric->GetlambdaDerivative()));
		metric->SetBinNumbers(NB);PMETRIC.push_back("Mi bin: ");VMETRIC.push_back(tostr(metric->GetBinNumbers()));
		metric->SetMAnumberOfSamples(numberOfSamplesMA);PMETRIC.push_back("MI numbe rof samples: ");VMETRIC.push_back(tostr(metric->GetMAnumberOfSamples()));
		metric->SetNumberOfThreads(nt); PMETRIC.push_back("Metric thread : ");VMETRIC.push_back(tostr(metric->GetNumberOfThreads()));
		metric->SetFixedEta(etaF); PMETRIC.push_back("Reference eta: ");VMETRIC.push_back(tostr(metric->GetFixedEta()));
		metric->SetMovingEta(etaM); PMETRIC.push_back("Source eta: ");VMETRIC.push_back(tostr(metric->GetMovingEta()));
		metric->SetEvaluator(NGFEvaluator);PMETRIC.push_back("Evaluator: ");VMETRIC.push_back(tostr(NGFEvaluator));

		pMetric("metric",PMETRIC,VMETRIC);

		registration->SetMetric( metric  );
	}else if (mmode==3) {

		typedef itk::MANGF2<
				InternalImageType,
				InternalImageType >   MetricType;
		MetricType::Pointer         metric        = MetricType::New();


		PMETRIC.push_back("Metric: ");VMETRIC.push_back("MI * lambda NGF");
		metric->Setlambda(lambda);	PMETRIC.push_back("Lambda: ");VMETRIC.push_back(tostr(metric->Getlambda()));
		metric->SetlambdaDerivative(lambdaderivative); PMETRIC.push_back("LambdaDerivative: ");VMETRIC.push_back(tostr(metric->GetlambdaDerivative()));
		metric->SetBinNumbers(NB);PMETRIC.push_back("Mi bin: ");VMETRIC.push_back(tostr(metric->GetBinNumbers()));
		metric->SetMAnumberOfSamples(numberOfSamplesMA);PMETRIC.push_back("MI numbe rof samples: ");VMETRIC.push_back(tostr(metric->GetMAnumberOfSamples()));
		metric->SetNumberOfThreads(nt); PMETRIC.push_back("Metric thread : ");VMETRIC.push_back(tostr(metric->GetNumberOfThreads()));
		metric->SetFixedEta(etaF); PMETRIC.push_back("Reference eta: ");VMETRIC.push_back(tostr(metric->GetFixedEta()));
		metric->SetMovingEta(etaM); PMETRIC.push_back("Source eta: ");VMETRIC.push_back(tostr(metric->GetMovingEta()));
		metric->SetEvaluator(NGFEvaluator);PMETRIC.push_back("Evaluator: ");VMETRIC.push_back(tostr(NGFEvaluator));

		pMetric("metric",PMETRIC,VMETRIC);
		registration->SetMetric( metric  );

//MASK
		if(argc>FMASK)
		if(*argv[FMASK]!=no)
		{
			typedef itk::ImageMaskSpatialObject< Dimension > MaskType;
			MaskType::Pointer spatialObjectMask = MaskType::New();
			typedef itk::Image< unsigned char, Dimension > ImageMaskType;
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
		if(*argv[MMASK]!=no)
		{
			typedef itk::ImageMaskSpatialObject< Dimension > MaskType;
			MaskType::Pointer spatialObjectMask = MaskType::New();
			typedef itk::Image< unsigned char, Dimension > ImageMaskType;
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
		
	}else if (mmode==4) {
		typedef itk::MANGFMSE<
				InternalImageType,
				InternalImageType >   MetricType;
		MetricType::Pointer         metric        = MetricType::New();


		PMETRIC.push_back("Metric: ");VMETRIC.push_back("MI * lambda NGF + mu MSE");
		metric->Setlambda(lambda);	PMETRIC.push_back("Lambda: ");VMETRIC.push_back(tostr(metric->Getlambda()));
		metric->SetlambdaDerivative(lambdaderivative); PMETRIC.push_back("LambdaDerivative: ");VMETRIC.push_back(tostr(metric->GetlambdaDerivative()));
		metric->SetBinNumbers(NB);PMETRIC.push_back("Mi bin: ");VMETRIC.push_back(tostr(metric->GetBinNumbers()));
		metric->SetMAnumberOfSamples(numberOfSamplesMA);PMETRIC.push_back("MI numbe rof samples: ");VMETRIC.push_back(tostr(metric->GetMAnumberOfSamples()));
		metric->SetNumberOfThreads(nt); PMETRIC.push_back("Metric thread : ");VMETRIC.push_back(tostr(metric->GetNumberOfThreads()));
		metric->SetFixedEta(etaF); PMETRIC.push_back("Reference eta: ");VMETRIC.push_back(tostr(metric->GetFixedEta()));
		metric->SetMovingEta(etaM); PMETRIC.push_back("Source eta: ");VMETRIC.push_back(tostr(metric->GetMovingEta()));
		metric->SetEvaluator(NGFEvaluator);PMETRIC.push_back("Evaluator: ");VMETRIC.push_back(tostr(NGFEvaluator));
		metric->Setmu(mu); PMETRIC.push_back("mu : ");VMETRIC.push_back(tostr(metric->Getmu()));
		metric->SetmuDerivative(muderivatives); PMETRIC.push_back("mu derivatives: ");VMETRIC.push_back(tostr(metric->GetmuDerivative()));
		registration->SetMetric( metric  );
		pMetric("metric",PMETRIC,VMETRIC);
		//MASK
		if(argc>FMASK)
		if(*argv[FMASK]!=no)
		{
			typedef itk::ImageMaskSpatialObject< Dimension > MaskType;
			MaskType::Pointer spatialObjectMask = MaskType::New();
			typedef itk::Image< unsigned char, Dimension > ImageMaskType;
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
		if(*argv[MMASK]!=no)
		{
			typedef itk::ImageMaskSpatialObject< Dimension > MaskType;
			MaskType::Pointer spatialObjectMask = MaskType::New();
			typedef itk::Image< unsigned char, Dimension > ImageMaskType;
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

	}else if (mmode==5){
		typedef itk::MeanSquaresImageToImageMetric<
				InternalImageType,
				InternalImageType >   MetricType;
		MetricType::Pointer         metric        = MetricType::New();
		registration->SetMetric( metric  );
		metric->SetNumberOfSpatialSamples(numberOfSamplesMA); PMETRIC.push_back("Number of samples :");VMETRIC.push_back(tostr(metric->GetNumberOfSpatialSamples()));
		pMetric("metric",PMETRIC,VMETRIC);
		//MASK
		if(argc>FMASK)
		if(*argv[FMASK]!=no)
		{
			typedef itk::ImageMaskSpatialObject< Dimension > MaskType;
			MaskType::Pointer spatialObjectMask = MaskType::New();
			typedef itk::Image< unsigned char, Dimension > ImageMaskType;
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
		if(*argv[MMASK]!=no)
		{
			typedef itk::ImageMaskSpatialObject< Dimension > MaskType;
			MaskType::Pointer spatialObjectMask = MaskType::New();
			typedef itk::Image< unsigned char, Dimension > ImageMaskType;
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
	}


	

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


	ParametersType finalParameters = registration->GetLastTransformParameters();

	unsigned int numberOfIterations = optimizer->GetCurrentIteration();

	double bestValue = optimizer->GetValue();


	// Print out results
	//
	std::cout << "Result = " << std::endl;
	std::cout << " Translation  = " << finalParameters << std::endl;
	std::cout << " Iterations    = " << numberOfIterations << std::endl;
	std::cout << " Metric value  = " << bestValue          << std::endl;

	typedef itk::ResampleImageFilter<
			MovingImageType,
			FixedImageType >    ResampleFilterType;

	TransformType::Pointer finalTransform = TransformType::New();

	finalTransform->SetParameters( finalParameters );
	finalTransform->SetFixedParameters( transform->GetFixedParameters() );

	ResampleFilterType::Pointer resample = ResampleFilterType::New();

	resample->SetTransform( finalTransform );
	resample->SetInput( movingImageReader->GetOutput() );

	FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();

	PixelType backgroundGrayLevel = 0;

	resample->SetOutputParametersFromImage(fixedImageReader->GetOutput());
	resample->SetDefaultPixelValue( backgroundGrayLevel );




	typedef itk::Image< OutputPixelType, Dimension > OutputImageType;

	typedef itk::CastImageFilter<
			FixedImageType,
			OutputImageType > CastFilterType;

	typedef itk::ImageFileWriter< OutputImageType >  WriterType;


	WriterType::Pointer      writer =  WriterType::New();
	CastFilterType::Pointer  caster =  CastFilterType::New();


	writer->SetFileName( argv[O] );


	caster->SetInput( resample->GetOutput() );
	writer->SetInput( caster->GetOutput()   );
	writer->Update();



	if (argc>OVF)
		if (*argv[OVF]!=no)
		{
			typedef itk::Vector< float,  Dimension >  VectorType;
			typedef itk::Image< VectorType,  Dimension >   OutputTransformationImageType;
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


if (argc>TOUT)
			if (*argv[TOUT]!=no)
		{
#if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 5) || ITK_VERSION_MAJOR > 4
  itk::TransformFileWriterTemplate<double>::Pointer writer =
    itk::TransformFileWriterTemplate<double>::New();
#else
  itk::TransformFileWriter::Pointer writer = itk::TransformFileWriter::New();
#endif

  writer->SetInput(registration->GetOutput()->Get());
  writer->SetFileName(argv[TOUT]);
  writer->Update();
		};


	return EXIT_SUCCESS;
}

