// #multiple optimizers
#include "itkImageRegistrationMethodv4.h"
#include "itkMeanSquaresImageToImageMetricv4.h"
#include "itkVersorRigid3DTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkRegularStepGradientDescentOptimizerv4.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkConjugateGradientOptimizer.h"
#include "itkCommand.h"
class CommandIterationUpdate : public itk::Command
{
public:
  using Self = CommandIterationUpdate;
  using Superclass = itk::Command;
  using Pointer = itk::SmartPointer<Self>;
  itkNewMacro(Self);
 
protected:
  CommandIterationUpdate() = default;
 
public:
  using OptimizerType = itk::RegularStepGradientDescentOptimizerv4<double>;
  using OptimizerPointer = const OptimizerType *;
  void
  Execute(itk::Object * caller, const itk::EventObject & event) override
  {
    Execute((const itk::Object *)caller, event);
  }
  void
  Execute(const itk::Object * object, const itk::EventObject & event) override
  {
    auto optimizer = static_cast<OptimizerPointer>(object);
    if (!itk::IterationEvent().CheckEvent(&event))
    {
      return;
    }
    std::cout << optimizer->GetCurrentIteration() << "   ";
    std::cout << optimizer->GetValue() << "   ";
    std::cout << optimizer->GetCurrentPosition() << std::endl;
  }
};
 
int
main(int argc, char * argv[])
{
  if (argc < 4)
  {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedImageFile  movingImageFile ";
    std::cerr << " outputImagefile  optimizer (rsgd,amoeba) " << std::endl;
    return EXIT_FAILURE;
  }
  constexpr unsigned int Dimension = 3;
  using PixelType = float;
  using FixedImageType = itk::Image<PixelType, Dimension>;
  using MovingImageType = itk::Image<PixelType, Dimension>;
 
  using TransformType = itk::VersorRigid3DTransform<double>;
 
using OptimizerType1 = itk::RegularStepGradientDescentOptimizerv4<double>;
using OptimizerType2 = itk::ConjugateGradientOptimizer ;

// Declare a pointer to the base optimizer class
itk::ObjectToObjectOptimizerBase::Pointer optimizer;
// Get the optimizer type from the user
std::string optimizerType = argv[4];

// Create the optimizer based on the user's choice
if (optimizerType == "rsgd") {
    optimizer = OptimizerType1::New();
    auto specificOptimizer = dynamic_cast<OptimizerType1*>(optimizer.GetPointer());
    specificOptimizer->SetLearningRate(0.2);
      specificOptimizer->SetMinimumStepLength(0.001);
    specificOptimizer->SetReturnBestParametersAndValue(true);
} else if (optimizerType == "conj") {
    optimizer = OptimizerType2::New();
    auto specificOptimizer = dynamic_cast<OptimizerType2*>(optimizer.GetPointer());

} else {
    std::cerr << "Unknown optimizer type: " << optimizerType << std::endl;
    return EXIT_FAILURE;
}

  using MetricType =
    itk::MeanSquaresImageToImageMetricv4<FixedImageType, MovingImageType>;
  using RegistrationType = itk::
    ImageRegistrationMethodv4<FixedImageType, MovingImageType, TransformType>;
 
  auto metric = MetricType::New();
  auto registration = RegistrationType::New();
 
  registration->SetMetric(metric);
  registration->SetOptimizer(optimizer);
 
  auto initialTransform = TransformType::New();
 
  using FixedImageReaderType = itk::ImageFileReader<FixedImageType>;
  using MovingImageReaderType = itk::ImageFileReader<MovingImageType>;
  auto fixedImageReader = FixedImageReaderType::New();
  auto movingImageReader = MovingImageReaderType::New();
 
  fixedImageReader->SetFileName(argv[1]);
  movingImageReader->SetFileName(argv[2]);
 
  registration->SetFixedImage(fixedImageReader->GetOutput());
  registration->SetMovingImage(movingImageReader->GetOutput());
 
 
  using TransformInitializerType =
    itk::CenteredTransformInitializer<TransformType,
                                      FixedImageType,
                                      MovingImageType>;
  auto initializer = TransformInitializerType::New();
  initializer->SetTransform(initialTransform);
  initializer->SetFixedImage(fixedImageReader->GetOutput());
  initializer->SetMovingImage(movingImageReader->GetOutput());
  initializer->MomentsOn();
  initializer->InitializeTransform();
  // Software Guide : EndCodeSnippet
 
 
  //  Software Guide : BeginLatex
  //
  //  The rotation part of the transform is initialized using a
  //  \doxygen{Versor} which is simply a unit quaternion.  The
  //  \code{VersorType} can be obtained from the transform traits. The versor
  //  itself defines the type of the vector used to indicate the rotation
  //  axis. This trait can be extracted as \code{VectorType}. The following
  //  lines create a versor object and initialize its parameters by passing a
  //  rotation axis and an angle.
  //
  //  Software Guide : EndLatex
 
  // Software Guide : BeginCodeSnippet
  using VersorType = TransformType::VersorType;
  using VectorType = VersorType::VectorType;
  VersorType rotation;
  VectorType axis;
  axis[0] = 0.0;
  axis[1] = 0.0;
  axis[2] = 1.0;
  constexpr double angle = 0;
  rotation.Set(axis, angle);
  initialTransform->SetRotation(rotation);

  registration->SetInitialTransform(initialTransform);
 
  using OptimizerScalesType = OptimizerType1::ScalesType;
  OptimizerScalesType optimizerScales(
    initialTransform->GetNumberOfParameters());
  const double translationScale = 1.0 / 1000.0;
  optimizerScales[0] = 1.0;
  optimizerScales[1] = 1.0;
  optimizerScales[2] = 1.0;
  optimizerScales[3] = translationScale;
  optimizerScales[4] = translationScale;
  optimizerScales[5] = translationScale;
  optimizer->SetScales(optimizerScales);
  optimizer->SetNumberOfIterations(200);

  auto observer = CommandIterationUpdate::New();
  optimizer->AddObserver(itk::IterationEvent(), observer);
 
  constexpr unsigned int numberOfLevels = 1;
 
  RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
  shrinkFactorsPerLevel.SetSize(1);
  shrinkFactorsPerLevel[0] = 1;
 
  RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
  smoothingSigmasPerLevel.SetSize(1);
  smoothingSigmasPerLevel[0] = 0;
 
  registration->SetNumberOfLevels(numberOfLevels);
  registration->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel);
  registration->SetShrinkFactorsPerLevel(shrinkFactorsPerLevel);
 
  try
  {
    registration->Update();
    std::cout << "Optimizer stop condition: "
              << registration->GetOptimizer()->GetStopConditionDescription()
              << std::endl;
  }
  catch (const itk::ExceptionObject& err) {
    
    return EXIT_FAILURE;
}
 
  const TransformType::ParametersType finalParameters =
    registration->GetOutput()->Get()->GetParameters();
 
  const double       versorX = finalParameters[0];
  const double       versorY = finalParameters[1];
  const double       versorZ = finalParameters[2];
  const double       finalTranslationX = finalParameters[3];
  const double       finalTranslationY = finalParameters[4];
  const double       finalTranslationZ = finalParameters[5];
  const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
  const double       bestValue = optimizer->GetValue();
 
  // Print out results
  //
  std::cout << std::endl << std::endl;
  std::cout << "Result = " << std::endl;
  std::cout << " versor X      = " << versorX << std::endl;
  std::cout << " versor Y      = " << versorY << std::endl;
  std::cout << " versor Z      = " << versorZ << std::endl;
  std::cout << " Translation X = " << finalTranslationX << std::endl;
  std::cout << " Translation Y = " << finalTranslationY << std::endl;
  std::cout << " Translation Z = " << finalTranslationZ << std::endl;
  std::cout << " Iterations    = " << numberOfIterations << std::endl;
  std::cout << " Metric value  = " << bestValue << std::endl;
 
  auto finalTransform = TransformType::New();
 
  finalTransform->SetFixedParameters(
    registration->GetOutput()->Get()->GetFixedParameters());
  finalTransform->SetParameters(finalParameters);
 
  // Software Guide : BeginCodeSnippet
  TransformType::MatrixType matrix = finalTransform->GetMatrix();
  TransformType::OffsetType offset = finalTransform->GetOffset();
  std::cout << "Matrix = " << std::endl << matrix << std::endl;
  std::cout << "Offset = " << std::endl << offset << std::endl;
  using ResampleFilterType =
    itk::ResampleImageFilter<MovingImageType, FixedImageType>;
 
  auto resampler = ResampleFilterType::New();
 
  resampler->SetTransform(finalTransform);
  resampler->SetInput(movingImageReader->GetOutput());
 
  FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();
 
  resampler->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
  resampler->SetOutputOrigin(fixedImage->GetOrigin());
  resampler->SetOutputSpacing(fixedImage->GetSpacing());
  resampler->SetOutputDirection(fixedImage->GetDirection());
  resampler->SetDefaultPixelValue(100);
 
  using OutputPixelType = unsigned char;
  using OutputImageType = itk::Image<OutputPixelType, Dimension>;
  using CastFilterType =
    itk::CastImageFilter<FixedImageType, OutputImageType>;
  using WriterType = itk::ImageFileWriter<OutputImageType>;
 
  auto writer = WriterType::New();
  auto caster = CastFilterType::New();
 
  writer->SetFileName(argv[3]);
 
  caster->SetInput(resampler->GetOutput());
  writer->SetInput(caster->GetOutput());
  writer->Update();
 
 return EXIT_SUCCESS;
}