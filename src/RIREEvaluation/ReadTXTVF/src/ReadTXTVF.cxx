
#include <iostream>
#include <fstream>
#include <cstring>

#include "itkImage.h"
#include "itkResampleImageFilter.h"


#include "itkVector.h"
#include "itkImageFileReader.h"

#include "Version.h"
#include "../../../include/funzioniIO.h"
#include "../../../include/mymath.h"

const    unsigned int    Dimension = 3;
typedef  double  PixelType;
typedef itk::Image< PixelType, Dimension >  ScalarImageType;
typedef itk::ImageFileReader<ScalarImageType> ScalarImageTypeReader;

typedef itk::Vector<double, Dimension> VectorPixelType;
typedef itk::Image<VectorPixelType, Dimension>  VectorImageType;
typedef itk::ImageFileReader<VectorImageType> VectorImageTypeReader;

int main( int argc, char *argv[] )
{


	std::vector<std::string>Vparams;
	Vparams.push_back("program name");
	int T=1;Vparams.push_back("reference image");
	int VF=T+1; Vparams.push_back("VF");
	int TXTTESTPOISITIONPOINT=VF+1;Vparams.push_back("test points txt");
	int TXTOUT = TXTTESTPOISITIONPOINT + 1; Vparams.push_back("newpointpositions txt");


	if( argc < 2 )
	{
		std::cout<<argv[0]<<" ";
		printVersion();
		printParam(Vparams);
		std::cout << "min, mean, median,std,max error" << std::endl;
		firma();
		return EXIT_FAILURE;
	}

	ScalarImageTypeReader::Pointer s=ScalarImageTypeReader::New();

	s->SetFileName(argv[T]);
	s->Update();
	VectorImageTypeReader::Pointer v=VectorImageTypeReader::New();
	v->SetFileName(argv[VF]);
	v->Update();

	ScalarImageType::Pointer im=s->GetOutput();
	VectorImageType::Pointer vf=v->GetOutput();



	ScalarImageType::PointType C;
	ScalarImageType::PointType R;
	ScalarImageType::PointType N;

	ScalarImageType::IndexType c;
	VectorImageType::PixelType r;

	std::ofstream f;
	f.open(argv[TXTOUT]);
	f<<"X,Y,Z"<<std::endl;
	int numberofpixels;
	std::string Load;
	std::ifstream fin;
	fin.open(argv[TXTTESTPOISITIONPOINT]);
	double RE;
	std::vector<double>ARE;


	//std::cout<<im->GetDirection()<<std::endl;
	while (fin >> numberofpixels >> C[0] >> C[1]>>C[2]>>R[0]>>R[1]>>R[2])	{

		im->TransformPhysicalPointToIndex(C,c);
		r=vf->GetPixel(c);
		//std::cout<<"Pixel: "<<c <<"has value"<<r<<std::endl;


		N=C-r;
//		N[2]=C [2]-r[2];
		f<<N[0]<<","<<N[1]<<","<<N[2]<<std::endl;
				RE=0;
				for (int u=0;u<3;u++){
					RE+=std::pow(N[u]-R[u],2);
				};

				ARE.push_back(std::sqrt(RE));
				//std::cout<<C<< "went: " <<N<<" and should be: " <<R<<"error: "<< std::sqrt(RE)<<std::endl;
			//	std::cout<<"index"<<c<< "is: " <<C<<std::endl;


	}
	fin.close();
f.close();
std::cout << CALC::min(ARE) << " " << CALC::mean(ARE) << " " << CALC::median(ARE) << " "<<CALC::std(ARE) << " " << CALC::max(ARE) << std::endl;




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

