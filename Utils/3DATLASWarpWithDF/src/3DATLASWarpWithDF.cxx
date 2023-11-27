#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVectorResampleImageFilter.h"
#include "itkImage.h"
#include "itkVector.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkWarpImageFilter.h"


//image type
const     unsigned int   Dimension = 3;
typedef   unsigned short  PixelType;
typedef   itk::Image< PixelType, Dimension > ImageType;


int main(int argc, char * argv[])
{
int NT=2;
int FIXEDSPACE=3;
int no=-1;

if( argc < 2 )
  {
    std::cerr << "Usage: "
              << std::endl
              << argv[0]
              << "\n01 inputImageFile \n02 deformation file \n03 fixxedspace (-1 if not needed) \n04outputImageFile \n05numberofthread"
              << std::endl;

    return EXIT_FAILURE;
  }

if (argc>5)
NT=atoi(argv[5]);


  typedef   double          VectorComponentType;

  typedef   itk::Vector< VectorComponentType, Dimension >    VectorType;
  typedef   itk::Image< VectorType,  Dimension >   DeformationFieldType;

typedef   itk::ImageFileReader<ImageType> ImageReaderType;
typedef   itk::ImageFileReader<DeformationFieldType> DeformationReaderType;
typedef   itk::ImageFileWriter<ImageType> ImageWriterType;

ImageReaderType::Pointer imReader =    ImageReaderType::New();
DeformationReaderType::Pointer defReader =    DeformationReaderType::New();
ImageWriterType::Pointer imWriter =    ImageWriterType::New();
ImageReaderType::Pointer space = ImageReaderType::New();


imReader->SetFileName(argv[1]);
defReader->SetFileName(argv[2]);
imWriter->SetFileName(argv[4]);
if (argc>FIXEDSPACE)
 	  if(atof(argv[FIXEDSPACE])!=no){ space->SetFileName(argv[FIXEDSPACE]);space->Update();}


imReader->Update();
defReader->Update();



  typedef itk::WarpImageFilter< ImageType,ImageType,DeformationFieldType  >  WarpImageFilterType;
  WarpImageFilterType::Pointer warpImageFilter = WarpImageFilterType::New();
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, VectorComponentType >  InterpolatorType;

  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  warpImageFilter->SetInterpolator( interpolator );
 warpImageFilter->SetNumberOfThreads( NT);
  warpImageFilter->SetDisplacementField( defReader->GetOutput() );
  warpImageFilter->SetInput( imReader->GetOutput() );
 


 if (argc>FIXEDSPACE)
	  {
	  if(atof(argv[FIXEDSPACE])!=no){warpImageFilter->SetOutputParametersFromImage(space->GetOutput());
	  }else {warpImageFilter->SetOutputParametersFromImage(imReader->GetOutput());}
	  }  
  
  warpImageFilter->Update();



 // Write the output

  imWriter->SetInput (  warpImageFilter->GetOutput() );
  imWriter->Update();

  return EXIT_SUCCESS;
}

