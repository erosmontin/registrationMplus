#include "itkMattesMutualInformationImageToImageMetric.h"
#include "../Metrics/MANGFImageToImageMetric/Code/itkNormalizedGradientFieldImageToImageMetric.h"
#include "../Metrics/MANGFImageToImageMetric/Code/itkNGFMetricKernel.h"
#include "../Metrics/MANGFImageToImageMetric/itkMANGF.h"
#include "../Metrics/MANGF2ImageToImageMetric/itkMANGF2.h"
#include "../Metrics/MANGFMSEImageToImageMetric/itkMANGFMSE.h"

#include <iostream>
#include <algorithm> // for copy
#include <iterator> // for ostream_iterator
#include <vector>
#include<fstream>
#include<sstream>
#include<iostream>


std::string availableMetrics(){
std::string S=	"metric mode (0 MATTES, 1 NGF, 2 MANGF+, 3 MANG*, 4 MANGFMSE, MSE)";
	return S;
}
void qualemetrica(RegistrationType::Pointer r, int mmode, int NB, int numberOfSamplesMA, double lambda, double lambdaderivative, double etaF, double etaM, int NGFEvaluator, double mu, double muderivatives,int nt, std::vector<std::string>& PMETRIC, std::vector<std::string>& VMETRIC){
	if (mmode==0){
		typedef itk::MattesMutualInformationImageToImageMetric<
				InternalImageType,
				InternalImageType >   MetricType;
		MetricType::Pointer    metric        = MetricType::New();
		metric->SetNumberOfHistogramBins(NB);
		metric->SetNumberOfSpatialSamples(numberOfSamplesMA);
		metric->ReinitializeSeed();

		PMETRIC.push_back("Metric: ");VMETRIC.push_back("MI");
		PMETRIC.push_back("Number of histogram bin"); VMETRIC.push_back(tostr(metric->GetNumberOfHistogramBins()));
		PMETRIC.push_back("Spatial samples");VMETRIC.push_back(tostr(metric->GetNumberOfSpatialSamples()));

		r->SetMetric( metric  );

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

		r->SetMetric( metric  );

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
		r->SetMetric( metric  );
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

		r->SetMetric( metric  );
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
		r->SetMetric( metric  );
	}else if (mmode==5){
		typedef itk::MeanSquaresImageToImageMetric<
				InternalImageType,
				InternalImageType >   MetricType;
		MetricType::Pointer         metric        = MetricType::New();
		metric->SetNumberOfSpatialSamples(numberOfSamplesMA); PMETRIC.push_back("Number of samples :");VMETRIC.push_back(tostr(metric->GetNumberOfSpatialSamples()));
		r->SetMetric( metric  );
	}

	};

	void maschere(RegistrationType::Pointer r, int FMASK,int MMASK, int argc, char* argv[], char no){
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

				}
				spatialObjectMask->SetImage( maskReader->GetOutput() );
				r->GetMetric()->SetFixedImageMask( spatialObjectMask );
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

				}
				spatialObjectMask->SetImage( maskReader->GetOutput() );
				r->GetMetric()->SetMovingImageMask( spatialObjectMask );
				std::cerr<<"MovingMASK: "<< argv[MMASK]<<std::endl;
			}
	};




	void transformationWrite(RegistrationType::Pointer r, int OVF, int argc, char* argv[], char no){



				typedef itk::TransformToDeformationFieldSource< OutputTransformationImageType, double >TransformToDeformationFieldSourceType;
				TransformToDeformationFieldSourceType::Pointer td = TransformToDeformationFieldSourceType::New();
				td->SetOutputParametersFromImage(r->GetMovingImage());
				td->SetTransform( r->GetOutput()->Get() );

				typedef itk::ImageFileWriter< OutputTransformationImageType>TransformToDeformationFieldSourceWriterType;
				TransformToDeformationFieldSourceWriterType::Pointer rtd = TransformToDeformationFieldSourceWriterType::New();
				rtd->SetInput(td->GetOutput());
				rtd->SetFileName(argv[OVF]);
				rtd->Update();

};





	void transformedImageWrite(InternalImageType::Pointer R,InternalImageType::Pointer M,
			TransformType::Pointer transform,
			ResampleFilterType::Pointer resample,
			WriterType::Pointer writer,
			CastFilterType::Pointer caster, int O,int nt, int argc, char* argv[]){



		resample->SetTransform( transform );
		resample->SetInput( M);

		resample->SetSize(   R->GetLargestPossibleRegion().GetSize() );
		resample->SetOutputOrigin(  R->GetOrigin() );
		resample->SetOutputSpacing( R->GetSpacing() );
		resample->SetOutputDirection( R->GetDirection() );
		resample->SetDefaultPixelValue( 0 );
		resample->SetNumberOfThreads(nt);



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

		}

};

	/*void readAndCast(RegistrationType::Pointer registration,int T,int S,int nt,char * argv[]){
		FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
			MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

			fixedImageReader->SetFileName(  argv[T] );
			movingImageReader->SetFileName( argv[S] );


			fixedImageReader->Update();
			movingImageReader->Update();

			FixedCastFilterType::Pointer fixedCaster   = FixedCastFilterType::New();
			MovingCastFilterType::Pointer movingCaster = MovingCastFilterType::New();
			fixedCaster->SetNumberOfThreads(nt);
			movingCaster->SetNumberOfThreads(nt);

			fixedCaster->SetInput(  fixedImageReader->GetOutput() );
			movingCaster->SetInput( movingImageReader->GetOutput() );

			registration->SetFixedImage(    fixedCaster->GetOutput()    );
			registration->SetMovingImage(   movingCaster->GetOutput()   );


			fixedCaster->Update();
	};*/





