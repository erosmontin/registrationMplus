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

#include <itkNormalizedGradientFieldImageToImageMetric.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkImage.h>

#include <itkCenteredRigid2DTransform.h>
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
  if( argc < 5 )
    {
    cerr << "Missing Parameters " << endl;
    cerr << "Usage: " << argv[0];
    cerr << " ImageFile1   ";
    cerr << " ImageFile2   ";
    cerr << " metric=(scalar|cross|scdelta|delta|delta2) ";
    cerr << " OutputImage\n";
    return EXIT_FAILURE;
    }
  
  const   unsigned int    Dimension = 2;
  typedef unsigned char   PixelType;

  typedef Image< PixelType, Dimension >  TImage;
  typedef CenteredRigid2DTransform< double > TransformType;

  typedef Image< RGBPixel<float> >         Gradient3DType; 
  typedef Gradient3DType::PixelType         Gradient3DPixelType;        
  typedef Image< RGBPixel<unsigned char> > ColorImageType; 
  typedef ImageFileWriter<ColorImageType> ImageFileWriterType; 
  

  typedef NormalizedGradientFieldImageToImageMetric<TImage, TImage>    NGFMetricType;
  typedef NGFMetricType::MovingNGFType MovingNGFType; 
  typedef NGFMetricType::FixedNGFType FixedNGFType; 

  typedef itk:: LinearInterpolateImageFunction<TImage, double >    InterpolatorType;
  typedef NGFMetricType::MovingNGFType GradientType;  
  typedef GradientType::PixelType GradientPixelType;  
  typedef NGFMetricType::MovingNGFType::Pointer GradientPointerType;  

  NGFMetricType::Pointer      ngfmetric        = NGFMetricType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  TransformType::Pointer  transform = TransformType::New();
  
  typedef ImageFileReader< TImage  > ImageReaderType;

  ImageReaderType::Pointer  fixedImageReader  = ImageReaderType::New();
  ImageReaderType::Pointer  movingImageReader  = ImageReaderType::New();

  fixedImageReader->SetFileName(  argv[1] );
  fixedImageReader->Update();
  TImage::Pointer fixedImage = fixedImageReader->GetOutput();

  movingImageReader->SetFileName(  argv[2] );
  movingImageReader->Update();
  TImage::Pointer movingImage = movingImageReader->GetOutput();

  if (!strcmp("scalar",argv[3])) {
	  ngfmetric->SetEvaluator(new NGFScalarKernel<MovingNGFType, FixedNGFType>()); 
  }else if (!strcmp("cross",argv[3])) {
	  ngfmetric->SetEvaluator(new NGFCrossKernel<MovingNGFType, FixedNGFType>()); 
  }else if (!strcmp("delta2",argv[3])) {
	  ngfmetric->SetEvaluator(new NGFDelta2Kernel<MovingNGFType, FixedNGFType>()); 
  }else if (!strcmp("deltasq",argv[3])) {
	  ngfmetric->SetEvaluator(new NGFDeltaKernel<MovingNGFType, FixedNGFType>()); 
  }else if (strcmp("scdelta",argv[3])) {
	  cerr << "unsupported evaluator '"<< argv[3] << "' specified"; 
	  return -1; 
  } // else use standard = "scdelta" 

  ngfmetric->Initialize(fixedImage, movingImage, 
			transform, interpolator, 
			fixedImage->GetLargestPossibleRegion()); 

  typedef TImage::SpacingType    SpacingType;
  typedef TImage::PointType      OriginType;
  typedef TImage::RegionType     RegionType;
  typedef TImage::SizeType       SizeType;

  const SpacingType fixedSpacing = fixedImage->GetSpacing();
  const OriginType  fixedOrigin  = fixedImage->GetOrigin();
  const RegionType  fixedRegion  = fixedImage->GetLargestPossibleRegion(); 
  const SizeType    fixedSize    = fixedRegion.GetSize();
  
  TransformType::InputPointType centerFixed;
  
  centerFixed[0] = fixedOrigin[0] + fixedSpacing[0] * fixedSize[0] / 2.0;
  centerFixed[1] = fixedOrigin[1] + fixedSpacing[1] * fixedSize[1] / 2.0;

  TransformType::InputPointType centerMoving;
  transform->SetCenter( centerFixed );
  cerr << centerFixed << "\n"; 
  transform->SetAngle( 0.0 );
	  
  TransformType::OutputVectorType trans = transform->GetTranslation(); 
  trans[0] =   trans[1]  = 0.0; 
  transform->SetTranslation(trans); 
  TransformType::ParametersType params = transform->GetParameters(); 
  
  // evaluate the metric gradient and save it as color image 
  GradientPointerType gradient = ngfmetric->GetGradient(params);

  // make a three channel image out of the two gradient channels
  Gradient3DType::Pointer conversion = Gradient3DType::New(); 
  Gradient3DType::RegionType gradientregion = gradient->GetLargestPossibleRegion(); 
  conversion->SetRegions(gradientregion.GetSize());
  conversion->Allocate(); 
  
  ImageRegionConstIterator< GradientType > cgi(gradient, gradientregion); 
  ImageRegionIterator< Gradient3DType >    oci(conversion, gradientregion); 
  
  cgi.GoToBegin(); 
  oci.GoToBegin(); 

  Gradient3DPixelType outpixel; 
  GradientPixelType   inpixel; 
  float minp= 0.0; 
  float maxp= 0.0; 
  float maxn= 0.0; 

  while (!cgi.IsAtEnd()) {
	  const GradientPixelType&   inpixel = cgi.Value(); 
	  outpixel[0] = inpixel[0]; 
	  outpixel[1] = inpixel[1];
	  outpixel[2] = log( 1.0 + sqrt (inpixel[0] * inpixel[0] + inpixel[1] * inpixel[1])); 

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
  ifwp->SetFileName(argv[4]); 
  ifwp->Update(); 
  
  return EXIT_SUCCESS;
}

