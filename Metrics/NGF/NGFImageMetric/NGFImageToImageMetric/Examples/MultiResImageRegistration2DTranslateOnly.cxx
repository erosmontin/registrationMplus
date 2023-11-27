/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: MultiResImageRegistration2.cxx,v $
  Language:  C++
  Date:      $Date: 2009-03-17 22:14:46 $
  Version:   $Revision: 1.50 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

//  Software Guide : BeginCommandLineArgs
//    INPUTS:  {BrainT1SliceBorder20.png}
//    INPUTS:  {BrainProtonDensitySliceShifted13x17y.png}
//    OUTPUTS: {MultiResImageRegistration2Output.png}
//    100
//    OUTPUTS: {MultiResImageRegistration2CheckerboardBefore.png}
//    OUTPUTS: {MultiResImageRegistration2CheckerboardAfter.png}
//  Software Guide : EndCommandLineArgs

// Software Guide : BeginLatex
//
//  This example illustrates the use of more complex components of the
//  registration framework. In particular, it introduces the use of the
//  \doxygen{AffineTransform} and the importance of fine-tuning the scale
//  parameters of the optimizer.
//
// \index{itk::ImageRegistrationMethod!AffineTransform}
// \index{itk::ImageRegistrationMethod!Scaling parameter space}
// \index{itk::AffineTransform!Image Registration}
//
// The AffineTransform is a linear transformation that maps lines into
// lines. It can be used to represent translations, rotations, anisotropic
// scaling, shearing or any combination of them. Details about the affine
// transform can be seen in Section~\ref{sec:AffineTransform}.
//
// In order to use the AffineTransform class, the following header
// must be included.
//
// \index{itk::AffineTransform!Header}
//
// Software Guide : EndLatex 


// Software Guide : BeginCodeSnippet
#include <itkTranslationTransform.h>
// Software Guide : EndCodeSnippet

#include <itkMultiResolutionImageRegistrationMethod.h>
#include <itkNormalizedGradientFieldImageToImageMetric.h>
#include <itkMeanSquaresImageToImageMetric.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkRegularStepGradientDescentOptimizer.h>
#include <itkRecursiveMultiResolutionPyramidImageFilter.h>
#include <itkImage.h>

#include <itkLevenbergMarquardtOptimizer.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkResampleImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkIntensityWindowingImageFilter.h>
#include <itkDifferenceImageFilter.h>

//  The following section of code implements an observer
//  that will monitor the evolution of the registration process.
//
#include <itkCommand.h>

using namespace itk; 

  typedef   RegularStepGradientDescentOptimizer     OptimizerType;

class CommandIterationUpdate : public Command 
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  Command             Superclass;
  typedef  SmartPointer<Self>  Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate(): m_CumulativeIterationIndex(0) {};
public:

  typedef   const OptimizerType   *           OptimizerPointer;

  void Execute(Object *caller, const EventObject & event)
    {
      Execute( (const Object *)caller, event);
    }

  void Execute(const Object * object, const EventObject & event)
    {
      OptimizerPointer optimizer = 
        dynamic_cast< OptimizerPointer >( object );
      if( !(IterationEvent().CheckEvent( &event )) )
        {
        return;
        }
      std::cout << optimizer->GetCurrentIteration() << "   ";
      std::cout << optimizer->GetValue() << "   ";
      std::cout << optimizer->GetCurrentPosition() << "  " <<
        m_CumulativeIterationIndex++ << std::endl;
    }
private:
  unsigned int m_CumulativeIterationIndex;
};


const    unsigned int    Dimension = 2;
typedef  float  PixelType;
typedef  unsigned char  OutputPixelType;

typedef Image< OutputPixelType, Dimension > OutputImageType;
typedef Image< PixelType, Dimension >  FixedImageType;
typedef Image< PixelType, Dimension >  MovingImageType;
typedef ShiftScaleImageFilter<MovingImageType, MovingImageType>   ShiftScaleType;
typedef IntensityWindowingImageFilter<MovingImageType, OutputImageType> WindowFilter; 


void writeColorOverlayImage(FixedImageType::Pointer fixed, MovingImageType::Pointer transformed, const char *fname)
{
	typedef Image< RGBPixel<float>, Dimension>  ColorImageType; 
	typedef Image< RGBPixel<unsigned char> > ColorImageOutputType; 
	typedef RescaleIntensityImageFilter<ColorImageType, ColorImageOutputType>   RescalerType;
	

	ColorImageType::Pointer colormix = ColorImageType::New(); 
	FixedImageType::RegionType region = fixed->GetLargestPossibleRegion(); 
	colormix->SetRegions(region);
	colormix->Allocate(); 

	ImageRegionIterator< ColorImageType > ci(colormix, region); 
	ImageRegionConstIterator< FixedImageType > fi(fixed, region); 
	ImageRegionConstIterator< MovingImageType > mi(transformed, region); 

	float maxcol[3] = {0,0,0}; 
	float mincol[3] = {0,0,0}; 
	while (!ci.IsAtEnd()) {
		ColorImageType::PixelType p; 
		p[0] = fi.Value(); 
		p[1] = mi.Value(); 
		p[2] = 0.5 * (fi.Value() + mi.Value());
		
		for(int i = 0; i < 3; ++i) {
			if (maxcol[i] < p[i]) 
				maxcol[i] = p[i]; 
			if (mincol[i] > p[i]) 
				mincol[i] = p[i]; 
		}
		ci.Set(p); 
		++fi; 
		++mi; 
		++ci; 
	}

	float scale[3]; 
	for(int i = 0; i < 3; ++i)
		scale[i] = 255.0/ (maxcol[i] - mincol[i]); 
		

	ColorImageOutputType::Pointer colorout = ColorImageOutputType::New(); 
	colorout->SetRegions(region);
	colorout->Allocate(); 

	ImageRegionConstIterator< ColorImageType > colin(colormix, region); 
	ImageRegionIterator< ColorImageOutputType > colout(colorout, region); 

	
	while (!colin.IsAtEnd()) {
		const ColorImageType::PixelType& p = colin.Value();  
		ColorImageOutputType::PixelType o; 
		for(int i = 0; i < 3; ++i)
			o[i] = (p[i] - mincol[i]) * scale[i]; 
		colout.Set(o); 
		++colin; 
		++colout; 
	}

	typedef ImageFileWriter< ColorImageOutputType >  WriterType;
	WriterType::Pointer      writer =  WriterType::New();
	
	writer->SetInput( colorout );
	writer->SetFileName( fname );
	writer->Update();

}

void writeDifferenceImage(FixedImageType::Pointer fixed, MovingImageType::Pointer transformed, const char *fname)
{
	typedef Image< float, Dimension>  DiffImageType; 
	typedef Image< unsigned char,  Dimension> DiffImageOutputType; 
	typedef IntensityWindowingImageFilter<DiffImageType, DiffImageOutputType> WindowFilter; 

	DiffImageType::Pointer diff = DiffImageType::New(); 
	FixedImageType::RegionType region = fixed->GetLargestPossibleRegion(); 
	diff->SetRegions(region);
	diff->Allocate(); 

	ImageRegionIterator< DiffImageType > ci(diff, region); 
	ImageRegionConstIterator< FixedImageType > fi(fixed, region); 
	ImageRegionConstIterator< MovingImageType > mi(transformed, region); 

	float sum = 0; 
	float maxval = 0; 
	float minval = 0; 
	int n = 0; 
	while (!ci.IsAtEnd()) {
		float p = mi.Value() - fi.Value(); 
		sum += p; 
		if (maxval < p) 
			maxval = p; 
		if (minval > p) 
			minval = p; 
		ci.Set(p); 
		++n; 
		++fi; 
		++mi; 
		++ci; 
	}


	const float avg = sum/n; 
	const float d1 = maxval - avg;
	const float d2 = avg - minval;
	const float d = d1 > d2 ? d1 : d2; 
	const float scale = 255.0/ (2.0 * d); 
	std::cerr << "min = " << minval << ", max = " << maxval; 
	std::cerr << "scale = " << scale << ", avg = " << avg << ", d = " << d << "\n"; 
	
	ci.GoToBegin(); 
	while (!ci.IsAtEnd()) {
		ci.Set((ci.Value() + d) * scale); 
		++ci; 
	}
	
	WindowFilter::Pointer windower = WindowFilter::New(); 
	windower->SetOutputMaximum(255.0); 
	windower->SetOutputMinimum(0.0);
	windower->SetWindowMinimum(0.0);
	windower->SetWindowMaximum(255.0); 
	windower->SetInput(diff); 
		

	typedef ImageFileWriter< DiffImageOutputType >  WriterType;
	WriterType::Pointer      writer =  WriterType::New();
	
	writer->SetInput( windower->GetOutput() );
	writer->SetFileName( fname );
	writer->Update();

}

//  The following section of code implements a Command observer
//  that will control the modification of optimizer parameters
//  at every change of resolution level.
//
template <typename TRegistration>
class RegistrationInterfaceCommand : public Command 
{
public:
  typedef  RegistrationInterfaceCommand   Self;
  typedef  Command                   Superclass;
  typedef  SmartPointer<Self>        Pointer;
  itkNewMacro( Self );
protected:
  RegistrationInterfaceCommand() {};
public:
  typedef   TRegistration                              RegistrationType;
  typedef   RegistrationType *                         RegistrationPointer;
  typedef   OptimizerType *                            OptimizerPointer;
  void Execute(Object * object, const EventObject & event)
  {
    if( !(IterationEvent().CheckEvent( &event )) )
      {
      return;
      }
    RegistrationPointer registration =
                        dynamic_cast<RegistrationPointer>( object );
    OptimizerPointer optimizer = dynamic_cast< OptimizerPointer >( 
                       registration->GetOptimizer() );

    std::cout << "-------------------------------------" << std::endl;
    std::cout << "MultiResolution Level : "
              << registration->GetCurrentLevel()  << std::endl;
    std::cout << std::endl;

    if ( registration->GetCurrentLevel() == 0 )
      {
      optimizer->SetMaximumStepLength( 0.1 );
      optimizer->SetMinimumStepLength( 0.01 );
      }
    else
      {
      optimizer->SetMaximumStepLength( optimizer->GetMaximumStepLength() / 0.5 );
      optimizer->SetMinimumStepLength( optimizer->GetMinimumStepLength() / 0.5 );
      }
  }
  void Execute(const Object * , const EventObject & )
    { return; }
};



int main( int argc, char *argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedImageFile  movingImageFile ";
    std::cerr << " outputImagefile [metric]" << std::endl;
    return EXIT_FAILURE;
    }
  
  typedef float     InternalPixelType;
  typedef Image< InternalPixelType, Dimension > InternalImageType;
  typedef TranslationTransform< double, 2 > TransformType;
  typedef LinearInterpolateImageFunction< InternalImageType, double> InterpolatorType;
  typedef NormalizedGradientFieldImageToImageMetric< InternalImageType, InternalImageType >  NGFMetricType;

  typedef NGFMetricType::MovingNGFType MovingNGFType; 
  typedef NGFMetricType::FixedNGFType  FixedNGFType; 
  typedef OptimizerType::ScalesType       OptimizerScalesType;

  typedef MultiResolutionImageRegistrationMethod< 
                                    InternalImageType, 
                                    InternalImageType    > RegistrationType;
  typedef RecursiveMultiResolutionPyramidImageFilter<
                                    InternalImageType,
                                    InternalImageType  >    FixedImagePyramidType;
  typedef RecursiveMultiResolutionPyramidImageFilter<
                                    InternalImageType,
                                    InternalImageType  >   MovingImagePyramidType;

  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();

  registration->SetOptimizer(     optimizer     );
  registration->SetInterpolator(  interpolator  );

  TransformType::Pointer   transform  = TransformType::New();
  registration->SetTransform( transform );

  
  FixedImagePyramidType::Pointer fixedImagePyramid = FixedImagePyramidType::New();
  MovingImagePyramidType::Pointer movingImagePyramid = MovingImagePyramidType::New();

  


  typedef ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef ImageFileReader< MovingImageType > MovingImageReaderType;

  FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

  fixedImageReader->SetFileName(  argv[1] );
  movingImageReader->SetFileName( argv[2] );


  if( argc > 4 ) {
	  if (!strcmp("scalar",argv[4])) {
		  NGFMetricType::Pointer metric = NGFMetricType::New();
		  metric->SetEvaluator(new NGFScalarKernel<MovingNGFType, FixedNGFType>()); 
		  registration->SetMetric( metric  );
	  }else if (!strcmp("cross",argv[4])) {
		  NGFMetricType::Pointer metric = NGFMetricType::New();
		  metric->SetEvaluator(new NGFCrossKernel<MovingNGFType, FixedNGFType>()); 
		  registration->SetMetric( metric  );
	  }else if (!strcmp("scdelta",argv[4])) {
		  NGFMetricType::Pointer metric = NGFMetricType::New();
		  metric->SetEvaluator(new NGFScaledDeltaKernel<MovingNGFType, FixedNGFType>()); 
		  registration->SetMetric( metric  );
	  }else if (!strcmp("ssd",argv[4])) {
		  registration->SetMetric(  MeanSquaresImageToImageMetric<InternalImageType, InternalImageType>::New()); 
	  }else  {
		  std::cerr << "Unsupported kernel = " << argv[4] << "\n"; 
		  return -1; 
	  }
  }

  typedef CastImageFilter< FixedImageType, InternalImageType > FixedCastFilterType;
  typedef CastImageFilter< MovingImageType, InternalImageType > MovingCastFilterType;
  FixedCastFilterType::Pointer fixedCaster   = FixedCastFilterType::New();
  MovingCastFilterType::Pointer movingCaster = MovingCastFilterType::New();

  fixedCaster->SetInput(  fixedImageReader->GetOutput() );
  movingCaster->SetInput( movingImageReader->GetOutput() );

  registration->SetFixedImage(    fixedCaster->GetOutput()    );
  registration->SetMovingImage(   movingCaster->GetOutput()   );

 
  fixedCaster->Update();
  registration->SetInitialTransformParameters( transform->GetParameters() );
  registration->SetFixedImageRegion( fixedCaster->GetOutput()->GetBufferedRegion() );

  MovingImageType::Pointer movingImage = movingImageReader->GetOutput();
  const FixedImageType::RegionType  movingRegion  = movingImage->GetLargestPossibleRegion();
  const FixedImageType::SizeType    movingSize    = movingRegion.GetSize();
  FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();

  typedef OptimizerType::ScalesType       OptimizerScalesType;
  OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );
  const double translationScale = 1.0;

  optimizerScales[0] = translationScale;
  optimizerScales[1] = translationScale;

  optimizer->SetScales( optimizerScales );
  optimizer->SetNumberOfIterations(  1000  );
  optimizer->SetRelaxationFactor( 0.8 );


  // Create the Command observer and register it with the optimizer.
  //
  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( IterationEvent(), observer );


  // Create the Command interface observer and register it with the optimizer.
  //
  typedef RegistrationInterfaceCommand<RegistrationType> CommandType;
  CommandType::Pointer command = CommandType::New();
  registration->AddObserver( IterationEvent(), command );
  registration->SetNumberOfLevels( 3 );

  try 
    { 
    registration->StartRegistration(); 
    } 
  catch( ExceptionObject & err ) 
    { 
    std::cout << "ExceptionObject caught !" << std::endl; 
    std::cout << err << std::endl; 
    return EXIT_FAILURE;
    } 

  std::cout << "Optimizer Stopping Condition = " 
            << optimizer->GetStopCondition() << std::endl;

  typedef RegistrationType::ParametersType ParametersType;
  ParametersType finalParameters = registration->GetLastTransformParameters();
  
  double TranslationAlongX = finalParameters[4];
  double TranslationAlongY = finalParameters[5];
  
  unsigned int numberOfIterations = optimizer->GetCurrentIteration();
  
  double bestValue = optimizer->GetValue();


  // Print out results
  //
  std::cout << "Result = " << std::endl;
  std::cout << " Translation X = " << TranslationAlongX  << std::endl;
  std::cout << " Translation Y = " << TranslationAlongY  << std::endl;
  std::cout << " Iterations    = " << numberOfIterations << std::endl;
  std::cout << " Metric value  = " << bestValue          << std::endl;

  typedef ResampleImageFilter< 
                            MovingImageType, 
                            FixedImageType >    ResampleFilterType;

  TransformType::Pointer finalTransform = TransformType::New();

  finalTransform->SetParameters( finalParameters );
  finalTransform->SetFixedParameters( transform->GetFixedParameters() );

  ResampleFilterType::Pointer resample = ResampleFilterType::New();

  resample->SetTransform( finalTransform );
  resample->SetInput( movingImageReader->GetOutput() );

  
  PixelType backgroundGrayLevel = 100;

  resample->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
  resample->SetOutputOrigin(  fixedImage->GetOrigin() );
  resample->SetOutputSpacing( fixedImage->GetSpacing() );
  resample->SetOutputDirection( fixedImage->GetDirection() );
  resample->SetDefaultPixelValue( backgroundGrayLevel );


  typedef  unsigned char  OutputPixelType;
  typedef Image< OutputPixelType, Dimension > OutputImageType;
  typedef CastImageFilter< 
                        FixedImageType,
                        OutputImageType > CastFilterType;
  typedef ImageFileWriter< OutputImageType >  WriterType;

  WriterType::Pointer      writer =  WriterType::New();
  CastFilterType::Pointer  caster =  CastFilterType::New();

  writer->SetFileName( argv[3] );

  caster->SetInput( resample->GetOutput() );
  writer->SetInput( caster->GetOutput()   );
  writer->Update();

  if( argc > 5 ) {
	  try {	  
		  writeColorOverlayImage(fixedImageReader->GetOutput(),
				       resample->GetOutput(), argv[5]); 
	  }
	  catch( ExceptionObject & excp ){ 
		  std::cerr << "ExceptionObject while writing the resampled image !" << std::endl; 
		  std::cerr << excp << std::endl; 
		  return EXIT_FAILURE;
	  } 
  }

  if( argc > 6 ) {
	  try {	  
		  writeColorOverlayImage(fixedImageReader->GetOutput(),
				       movingImageReader->GetOutput(), argv[6]); 
	  }
	  catch( ExceptionObject & excp ){ 
		  std::cerr << "ExceptionObject while writing the resampled image !" << std::endl; 
		  std::cerr << excp << std::endl; 
		  return EXIT_FAILURE;
	  } 
  }

  return EXIT_SUCCESS;
}

