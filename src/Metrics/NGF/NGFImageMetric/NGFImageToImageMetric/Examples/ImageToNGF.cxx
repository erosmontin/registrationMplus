/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ImageRegistration5.cxx,v $
  Language:  C++
  Date:      $Date: 2007-11-22 00:30:16 $
  Version:   $Revision: 1.44 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif


#include <itkImage.h>
#include <itkImageToNGFFilter.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkCastImageFilter.h>

using namespace std; 
using namespace itk; 

/*
  This Program evaluates a gradient strength map of the metric between two images
*/

int main( int argc, char *argv[] )
{
  if( argc < 3 )
    {
    cerr << "Missing Parameters " << endl;
    cerr << "Usage: " << argv[0];
    cerr << " ImageFile1   ";
    cerr << " ImageFile2   ";
    return EXIT_FAILURE;
    }
  
  const   unsigned int    Dimension = 2;
  typedef unsigned char   PixelType;

  typedef Image< PixelType, Dimension >  TImage;
  typedef Image< RGBPixel<float> >         Gradient3DType; 
  typedef Gradient3DType::PixelType         Gradient3DPixelType;        
  typedef Image< RGBPixel<unsigned char> > ColorImageType; 
  typedef ImageFileWriter<ColorImageType> ImageFileWriterType; 
  

  typedef ImageFileReader< TImage  > ImageReaderType;

  ImageReaderType::Pointer  fixedImageReader  = ImageReaderType::New();

  fixedImageReader->SetFileName(  argv[1] );

  ImageToNGFFilter<TImage>::Pointer image2NGFFilter = ImageToNGFFilter<TImage>::New(); 
  image2NGFFilter->SetInput(fixedImageReader->GetOutput()); 

  typedef ImageToNGFFilter<TImage>::OutputImageType GradientType; 
  typedef GradientType::PixelType GradientPixelType; 
  GradientType *imageGradient = image2NGFFilter->GetOutput(); 
  imageGradient->Update(); 
  
  
    // make a three channel image out of the two gradient channels
  Gradient3DType::Pointer conversion = Gradient3DType::New(); 
  Gradient3DType::RegionType gradientregion = imageGradient->GetLargestPossibleRegion(); 
  conversion->SetRegions(gradientregion.GetSize());
  conversion->Allocate(); 
  
  ImageRegionConstIterator< GradientType > cgi(imageGradient, gradientregion); 
  ImageRegionIterator< Gradient3DType >    oci(conversion, gradientregion); 
  
  cgi.GoToBegin(); 
  oci.GoToBegin(); 

  Gradient3DPixelType outpixel; 
  GradientPixelType   inpixel; 
  float minp= 0.0; 
  float maxp= 0.0; 
  float maxn= 0.0; 
  float sum = 0.0; 

  while (!cgi.IsAtEnd()) {
	  const GradientPixelType&   inpixel = cgi.Value(); 
	  outpixel[0] = inpixel[0]; 
	  outpixel[1] = inpixel[1];
	  outpixel[2] = sqrt (inpixel[0] * inpixel[0] + inpixel[1] * inpixel[1]); 
	  
	  if (outpixel[2] > 0.1) 
		  sum += 1.0; 

	  if (maxn < outpixel[2])
		  maxn = outpixel[2];

 	  if (minp > outpixel[0])
		  minp=  outpixel[0]; 
	  if (minp > outpixel[1])
		  minp=  outpixel[1]; 

	  if (maxp < outpixel[0])
		  maxp=  outpixel[0]; 
	  if (maxp < outpixel[1])
		  maxp=  outpixel[1]; 

	  oci.Set(outpixel);
	  ++oci; 
	  ++cgi; 
  }
  cout << "sum gradient norms:" << sum << "\n"; 
  
  if (-minp < maxp) 
	  minp =  - maxp; 
  else 
	  maxp = - minp; 

  oci.GoToBegin(); 
  
  float scale = 255.0 / (maxp - minp); 

  
  while (!oci.IsAtEnd()) {
	  Gradient3DPixelType outpixel;
	  const Gradient3DPixelType& inpixel = oci.Value(); 
	  outpixel[1] = inpixel[2] / maxn; 

	  outpixel[0] = outpixel[1] * scale * ( inpixel[0] - minp ); 
	  outpixel[2] = outpixel[1] * scale * ( inpixel[1] - minp ); 

	  outpixel[1] *= 255; 
	  oci.Set(outpixel); 
	  ++oci; 
  }

  // resample the 3-channel image to become a nice colored image of the gradient
  typedef CastImageFilter<Gradient3DType, ColorImageType> ColorCastFilter; 

  ColorCastFilter::Pointer crf = ColorCastFilter::New(); 
  crf->SetInput(conversion); 

  ImageFileWriterType::Pointer ifwp = ImageFileWriterType::New();
  
  ifwp->SetInput(crf->GetOutput()); 
  ifwp->SetFileName(argv[2]); 
  ifwp->Update(); 
  
  return EXIT_SUCCESS;
}

