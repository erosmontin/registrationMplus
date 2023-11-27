#include "itkImageFileReader.h"
#include "itkTranslationTransform.h"
#include "itkImageMomentsCalculator.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"

#include "itkImageRegistrationMethod.h"
#include "itkTransformFileWriter.h"
#include "itkCenteredTransformInitializer.h"
#include "itkVersorRigid3DTransform.h"
#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"
int
main(int argc, char * argv[])
{
	  if (argc < 2)
  {
    std::cout << "Usage: " << argv[0] << " fixedImage movingImage (g or m) outptfile "<< std::endl;
    return EXIT_FAILURE;
  }
    itk::TimeProbesCollectorBase   chronometer;
  itk::MemoryProbesCollectorBase memorymeter;
  using ImageType = itk::Image<double, 3>;
  using ImageReaderType = itk::ImageFileReader<ImageType>;
  ImageReaderType::Pointer     fixedImageReader = ImageReaderType::New();
  ImageReaderType::Pointer     movingImageReader = ImageReaderType::New();
    
  using ImageWriterType = itk::ImageFileWriter<ImageType>;
  ImageWriterType::Pointer     w = ImageWriterType::New();
  
  using TransformType = itk::TranslationTransform<double,3>;


  using MomentCalulatortype = itk::ImageMomentsCalculator<ImageType>;
  MomentCalulatortype::Pointer initializer1 = MomentCalulatortype::New();

  MomentCalulatortype::Pointer initializer2 = MomentCalulatortype::New();
  
  
  fixedImageReader->SetFileName(argv[1]);
  fixedImageReader->Update();
  
  ImageType::Pointer  fixedImage  = ImageType::New();
  fixedImage = fixedImageReader->GetOutput();
  
  initializer1->SetImage(fixedImage);
  initializer1->Compute();
  

  movingImageReader->SetFileName(argv[2]);
  movingImageReader->Update();
  
  ImageType::Pointer  movingImage  = ImageType::New();
  movingImage = movingImageReader->GetOutput();
  
  initializer2->SetImage(movingImage);
  initializer2->Compute();
  
  std::cout<< initializer1->GetCenterOfGravity() << initializer2->GetCenterOfGravity();

  using RigidTransformType = itk::VersorRigid3DTransform<double>;
  using TransformInitializerType =
    itk::CenteredTransformInitializer<RigidTransformType,
                                      ImageType,
                                      ImageType>;
    auto initializer = TransformInitializerType::New();

  auto rigidTransform = RigidTransformType::New();
 
  initializer->SetTransform(rigidTransform);
  initializer->SetFixedImage(fixedImageReader->GetOutput());
  initializer->SetMovingImage(movingImageReader->GetOutput());

  initializer->MomentsOn();
 
  if (argv[3]=="g"){
  std::cout<<"Geometry on!"<<std::endl;
  initializer->GeometryOn();
    }else{
      std::cout<<"Moment on!"<<std::endl;
    } 
  std::cout << "Starting Rigid Transform Initialization " << std::endl;
 
  memorymeter.Start("Rigid Initialization");
  chronometer.Start("Rigid Initialization");
 
  initializer->InitializeTransform();
 
  chronometer.Stop("Rigid Initialization");
  memorymeter.Stop("Rigid Initialization");
 
  std::cout << "Rigid Transform Initialization completed" << std::endl;
  std::cout << std::endl;

ImageType::RegionType fixedRegion = fixedImage->GetBufferedRegion();

  

 

 auto transform = TransformType::New();

RigidTransformType::ParametersType p= rigidTransform->GetParameters();
        std::cout<<rigidTransform->GetFixedParameters()<<std::endl;

        
 TransformType::OutputVectorType vector;
        vector[0] = p[3];
        vector[1] = p[4];
        vector[2] = p[5];

      transform->Translate(vector);
      
      
       using ResampleImageFilterType = itk::ResampleImageFilter<ImageType, ImageType>;
  auto resampleFilter = ResampleImageFilterType::New();
  resampleFilter->SetTransform(transform);
  resampleFilter->SetInput(movingImage);
 
  resampleFilter->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
  resampleFilter->SetOutputOrigin(fixedImage->GetOrigin());
  resampleFilter->SetOutputSpacing(fixedImage->GetSpacing());
  resampleFilter->SetOutputDirection(fixedImage->GetDirection());
 
  resampleFilter->SetDefaultPixelValue(0);
  
  
  try
  {
    w->SetInput(resampleFilter->GetOutput());
    w->SetFileName(argv[4]) ;
    w->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "Error: " << err << std::endl;
    return EXIT_FAILURE;
  }

   if (argc > 5)
  {
    std::cout << "Writing transform parameter file ...";
    using TransformWriterType = itk::TransformFileWriter;
    auto transformWriter = TransformWriterType::New();
    transformWriter->AddTransform(resampleFilter->GetTransform());
    transformWriter->SetFileName(argv[5]);
    transformWriter->Update();
    std::cout << " Done!" << std::endl;
  }
      
  
  
  
  
}
