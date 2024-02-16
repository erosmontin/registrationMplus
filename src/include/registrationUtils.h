#include <iostream>
#include <algorithm> // for copy
#include <iterator> // for ostream_iterator
#include <vector>
#include<fstream>
#include<sstream>
#include<iostream>

template <class RInternalImageType2,class MInternalImageType2,class TransformType2,class ResamplerType2,class WriterType2>
	void transformAndWriteImage(typename RInternalImageType2::Pointer R,typename MInternalImageType2::Pointer M,
			typename TransformType2::Pointer transform,
			typename ResamplerType2::Pointer resample,
			typename WriterType2::Pointer writer
			,char*  fname,int nt){

		resample->SetTransform( transform );
		resample->SetInput(M);
		resample->SetSize(R->GetLargestPossibleRegion().GetSize() );
		resample->SetOutputOrigin(  R->GetOrigin() );
		resample->SetOutputSpacing( R->GetSpacing() );
		resample->SetOutputDirection( R->GetDirection() );
		resample->SetDefaultPixelValue( 0 );
		resample->SetNumberOfThreads(nt);
		writer->SetFileName(fname);
		writer->SetInput(resample->GetOutput());
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


