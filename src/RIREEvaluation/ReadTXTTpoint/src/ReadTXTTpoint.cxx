
#include <iostream>
#include <fstream>
#include <cstring>

#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkTranslationTransform.h"

#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkImage.h"
#include "itkResampleImageFilter.h"

#include "itkVersion.h"
#include "itkTransformFileReader.h"


#include "Version.h"

#include "../../../include/funzioniIO.h"
#include "../../../include/mymath.h"

const    unsigned int    Dimension = 3;
typedef  short  PixelType;
typedef  short  OutputPixelType;

typedef itk::Image< PixelType, Dimension >  ImageType;

typedef itk::TranslationTransform< double, Dimension > TransformType;
typedef itk::RegularStepGradientDescentOptimizer       OptimizerType;
typedef itk::LinearInterpolateImageFunction<
		ImageType,
		double             > InterpolatorType;

typedef itk::MultiResolutionImageRegistrationMethod<
		ImageType,
		ImageType >   RegistrationType;

typedef itk::MultiResolutionPyramidImageFilter<
		ImageType,
		ImageType >   FixedImagePyramidType;
typedef itk::MultiResolutionPyramidImageFilter<
		ImageType,
		ImageType >   MovingImagePyramidType;

typedef RegistrationType::ParametersType ParametersType;





int main( int argc, char *argv[] )
{


	std::vector<std::string>Vparams;
	Vparams.push_back("program name");
	int T=1;Vparams.push_back("tfm txt");
	int TXTTESTPOISITIONPOINT=T+1;Vparams.push_back("test points txt");
	int TXTREALPOSITIONPOINT=TXTTESTPOISITIONPOINT+1;Vparams.push_back("Real point txt");


	if( argc < 2 )
	{
		std::cout<<argv[0]<<" ";
		printVersion();
		printParam(Vparams);

		firma();
		return EXIT_FAILURE;
	}



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

	ParametersType Parameters( transform->GetNumberOfParameters() );




#if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 5) || ITK_VERSION_MAJOR > 4
	itk::TransformFileReaderTemplate<double>::Pointer reader =
			itk::TransformFileReaderTemplate<double>::New();
#else
	itk::TransformFileReader::Pointer reader = itk::TransformFileReader::New();
#endif


	reader->SetFileName(argv[T]);
	reader->Update();
	Parameters= reader->GetTransformList()->begin()->GetPointer()->GetParameters();

	transform->SetParameters(Parameters);


	ImageType::PointType C;
	ImageType::PointType R;
	TransformType::OutputPointType N;


	int numberofpixels;
	std::string Load;
	std::ifstream fin;
	fin.open(argv[TXTTESTPOISITIONPOINT]);
	double RE;
	std::vector<double>ARE;

	while (fin >> numberofpixels >> C[0] >> C[1]>>C[2]>>R[0]>>R[1]>>R[2])	{

				N = transform->TransformPoint(C);
				RE=0;
				for (int u=0;u<3;u++){
					RE+=std::pow(N[u]-R[u],2);
				};
				ARE.push_back(std::sqrt(RE));
				//std::cout<<C<< "went: " <<N<<" and should be: " <<R<<"error: "<< std::sqrt(RE)<<std::endl;



	}
	fin.close();

	std::cout<<CALC::max(ARE)<<" "<< CALC::mean(ARE)<<" "<<CALC::median(ARE) <<std::endl;




		//for(int t=0;t<Dimension;t++){C[t]=atof(argv[P+t]);};
		//
		//TransformType::OutputPointType trackerPointNewPosition;
		//  trackerPointNewPosition = transform->TransformPoint(C);
		//
		//
		//  std::cout<<trackerPointNewPosition[0]<<
		//		  "\t" <<trackerPointNewPosition[1]<<
		//		  "\t" <<trackerPointNewPosition[2]<<std::endl;



		return EXIT_SUCCESS;
	}

