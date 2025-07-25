#include <itkImage.h>
#include <itkCommand.h>
#include <itkLBFGSBOptimizer.h>
#include <chrono>    // << add this
#include <iomanip>                       // << for std::setprecision

#include <iostream>

class LBFGSBOptimizeCommandIterationUpdate : public itk::Command
{
public:
    using Clock = std::chrono::steady_clock;    // << add this alias
    typedef  LBFGSBOptimizeCommandIterationUpdate   Self;
    typedef  itk::Command             Superclass;
    typedef itk::SmartPointer<Self>   Pointer;
    itkNewMacro( Self );
protected:
    LBFGSBOptimizeCommandIterationUpdate()
        : m_StartTime( Clock::now() )              // << initialize start
    {};
private:
    Clock::time_point m_StartTime;               // << store start
public:
    typedef itk::LBFGSBOptimizer    OptimizerType;
    typedef   const OptimizerType * OptimizerPointer;
    void Execute(itk::Object *caller, const itk::EventObject & event)
    {
        Execute( (const itk::Object *)caller, event);
    }
    void Execute(const itk::Object * object, const itk::EventObject & event) override
    {
        OptimizerPointer optimizer = static_cast< OptimizerPointer >( object );
        if( !(itk::IterationEvent().CheckEvent( &event )) )
        {
            return;
        }

        // compute elapsed
        auto now     = Clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - m_StartTime);
        double secs  = elapsed.count() / 1000.0;

        // print with elapsed time
        std::cout << "\r[" << std::fixed << std::setprecision(2)
                  << secs << std::defaultfloat << " s] Iter: "
                  << optimizer->GetCurrentIteration()
                  << "  Metric: " << optimizer->GetCachedValue()
                  << "  Inf-norm proj-grad: "
                  << optimizer->GetInfinityNormOfProjectedGradient()
                  << std::flush;
    }
};

#include <itkRegularStepGradientDescentOptimizer.h>
class RegularStepGradientDescentOptimizerCommandIterationUpdate : public itk::Command
{
public:
    using Clock = std::chrono::steady_clock;    // << add this alias
    typedef  RegularStepGradientDescentOptimizerCommandIterationUpdate Self;
    typedef  itk::Command                                             Superclass;
    typedef itk::SmartPointer<Self>                                   Pointer;
    itkNewMacro( Self );

    itkSetMacro(ShowGradient, bool);
    itkGetMacro(ShowGradient, bool);

protected:
    RegularStepGradientDescentOptimizerCommandIterationUpdate()
      : m_ShowGradient(false),
        m_StartTime( Clock::now() )               // << initialize start clock
    {};

private:
    bool                        m_ShowGradient;
    Clock::time_point           m_StartTime;       // << store start time

public:
    typedef itk::RegularStepGradientDescentOptimizer    OptimizerType;
    typedef const OptimizerType *                       OptimizerPointer;

    void Execute(itk::Object *caller, const itk::EventObject & event) override
    {
        this->Execute( static_cast<const itk::Object*>(caller), event );
    }

    void Execute(const itk::Object * object, const itk::EventObject & event) override
    {
        OptimizerPointer optimizer = static_cast< OptimizerPointer >( object );
        if( ! itk::IterationEvent().CheckEvent( &event ) )
        {
            return;
        }

        // compute elapsed
        auto now     = Clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - m_StartTime);
        double secs  = elapsed.count() / 1000.0;

        long iteration = optimizer->GetCurrentIteration();

        // overwrite the same line
        std::cout << "\r[" << std::fixed << std::setprecision(2)
                  << secs << std::defaultfloat << " s] Iter: " << iteration
                  << "  Value: " << optimizer->GetValue();

        if (m_ShowGradient)
        {
            std::cout << "  Grad: " << optimizer->GetGradient();
        }

        std::cout << "  GradTol: " << optimizer->GetGradientMagnitudeTolerance()
                  << "  Step: "    << optimizer->GetCurrentStepLength()
                  << std::flush;
    }
};

#include "itkParticleSwarmOptimizer.h"

class ParticleSwarmOptimizeCommandIterationUpdate : public itk::Command
{
public:
    typedef  ParticleSwarmOptimizeCommandIterationUpdate   Self;
    typedef  itk::Command             Superclass;
    typedef itk::SmartPointer<Self>   Pointer;
    itkNewMacro( Self );
protected:
    ParticleSwarmOptimizeCommandIterationUpdate() {};
public:
    typedef itk::ParticleSwarmOptimizer    OptimizerType;
    typedef   const OptimizerType * OptimizerPointer;
    void Execute(itk::Object *caller, const itk::EventObject & event)
    {
        Execute( (const itk::Object *)caller, event);
    }
    void Execute(const itk::Object * object, const itk::EventObject & event)
    {
        OptimizerPointer optimizer = static_cast< OptimizerPointer >( object );
        if( !(itk::IterationEvent().CheckEvent( &event )) )
        {
            return;
        }
    if(optimizer == nullptr)
    {
        std::cout << "Optimizer is null!" << std::endl;
    }
    else
    {
        try
        {
            std::cout << optimizer->GetCurrentPosition() << "   ";
            std::cout << optimizer->GetValue() << std::endl;
        }
        catch(const std::exception& e)
        {
            std::cout << "Caught exception: " << e.what() << std::endl;
        }
    }
    }
};



template <typename TRegistration>
class RegistrationInterfaceCommand : public itk::Command
{
public:
  using Self = RegistrationInterfaceCommand;
  using Superclass = itk::Command;
  using Pointer = itk::SmartPointer<Self>;
  itkNewMacro(Self);
 
protected:
  RegistrationInterfaceCommand() = default;
 
public:
  using RegistrationType = TRegistration;
  using RegistrationPointer = RegistrationType *;
  using OptimizerType = itk::RegularStepGradientDescentOptimizer;
  using OptimizerPointer = OptimizerType *;
  void
  Execute(itk::Object * object, const itk::EventObject & event) override
  {
    if (!(itk::IterationEvent().CheckEvent(&event)))
    {
      return;
    }
    auto registration = static_cast<RegistrationPointer>(object);
    auto optimizer =
      static_cast<OptimizerPointer>(registration->GetModifiableOptimizer());
 
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "MultiResolution Level : " << registration->GetCurrentLevel()
              << std::endl;
    std::cout << std::endl;
 
    if (registration->GetCurrentLevel() == 0)
    {
        optimizer->SetMaximumStepLength(0.1);
        optimizer->SetMinimumStepLength(0.01);
        }
        
    else
    {
      optimizer->SetMaximumStepLength(optimizer->GetMaximumStepLength() /
                                      4.0);
      optimizer->SetMinimumStepLength(optimizer->GetMinimumStepLength() /
                                      10.0);
    }
        if (optimizer->GetCurrentIteration() == 0)
    {
        std::cout << "-------------------------------------" << std::endl;
        std::cout << "MultiResolution Level : " << registration->GetCurrentLevel()
                  << std::endl;

        // Print the resolution at the current level
        auto fixedImagePyramid = registration->GetFixedImagePyramid();
        auto movingImagePyramid = registration->GetMovingImagePyramid();
        std::cout << "Fixed Image Resolution : " 
                  << fixedImagePyramid->GetOutput(registration->GetCurrentLevel())->GetSpacing() 
                  << std::endl;
        std::cout << "Moving Image Resolution : " 
                  << movingImagePyramid->GetOutput(registration->GetCurrentLevel())->GetSpacing() 
                  << std::endl;

        std::cout << std::endl;
    }
};

  void
  Execute(const itk::Object *, const itk::EventObject &) override
  {
    return;
  }
};


#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBSplineTransform.h"
#include "itkResampleImageFilter.h"
#include "itkTransformToDisplacementFieldFilter.h"

template <typename TImage, typename TTransform>
typename TImage::Pointer CreateBSplineFromTransform(TTransform* transform, TImage* image, typename TImage::SizeType gridSize)
{
    using ImageType = TImage;
    using TransformType = TTransform;
    using DisplacementFieldType = itk::Image<itk::Vector<float, ImageType::ImageDimension>, ImageType::ImageDimension>;
    using TransformToDisplacementFieldFilterType = itk::TransformToDisplacementFieldFilter<DisplacementFieldType, double>;

    // Create a deformation field from the transform
    typename TransformToDisplacementFieldFilterType::Pointer displacementFieldFilter = TransformToDisplacementFieldFilterType::New();
    displacementFieldFilter->SetTransform(transform);
    displacementFieldFilter->SetSize(image->GetLargestPossibleRegion().GetSize());
    displacementFieldFilter->SetOutputStartIndex(image->GetLargestPossibleRegion().GetIndex());
    displacementFieldFilter->SetOutputSpacing(image->GetSpacing());
    displacementFieldFilter->SetOutputDirection(image->GetDirection());
    displacementFieldFilter->SetOutputOrigin(image->GetOrigin());
    displacementFieldFilter->Update();

    typename DisplacementFieldType::Pointer displacementField = displacementFieldFilter->GetOutput();

    // Initialize a B-spline transform with the deformation field
    typename TransformType::Pointer bsplineTransform = TransformType::New();
    bsplineTransform->SetTransformDomainOrigin(image->GetOrigin());
    bsplineTransform->SetTransformDomainDirection(image->GetDirection());
    bsplineTransform->SetTransformDomainPhysicalDimensions(image->GetSpacing());
    bsplineTransform->SetTransformDomainMeshSize(gridSize);

    // Use the deformation field to set the weights of the B-spline transform's coefficients
    typename TransformType::ParametersType parameters(bsplineTransform->GetNumberOfParameters());
    for (unsigned int i = 0; i < parameters.GetSize(); ++i)
    {
        parameters[i] = displacementField->GetPixel(i);
    }
    bsplineTransform->SetParameters(parameters);

    return bsplineTransform;
}


#include "itkTransformToDeformationFieldSource.h"



template<typename TTransformType,typename TMovingImageType,typename TDeformationFieldType>
typename TDeformationFieldType::Pointer TransformToDeformationField(typename TTransformType::Pointer transform, 
                                                                             typename TMovingImageType::Pointer movingImage)
{
    typedef itk::TransformToDeformationFieldSource<TDeformationFieldType, typename TDeformationFieldType::PixelType::ValueType> TransformToDeformationFieldSourceType;
    typename TransformToDeformationFieldSourceType::Pointer td = TransformToDeformationFieldSourceType::New();
    td->SetOutputParametersFromImage(movingImage);
    td->SetTransform(transform);
    td->Update();

    return td->GetOutput();
}

#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"

template<typename TTransformType>
void WriteTransform(const std::string& fileName, typename TTransformType::Pointer transform)
{
    #if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 5) || ITK_VERSION_MAJOR > 4
        typedef itk::TransformFileWriterTemplate<typename TTransformType::ParametersValueType> WriterType;
    #else
        typedef itk::TransformFileWriter WriterType;
    #endif

    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput(transform);
    writer->SetFileName(fileName);
    writer->Update();
}

#include "itkTransformFileReader.h"

itk::TransformBase::Pointer ReadTransformGeneric(const std::string& fileName)
{
    #if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 5) || ITK_VERSION_MAJOR > 4
        typedef itk::TransformFileReaderTemplate<double> ReaderType;
    #else
        typedef itk::TransformFileReader ReaderType;
    #endif

    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(fileName);
    reader->Update();

    return reader->GetTransformList()->front();
}

#include "itkTransformFileReader.h"

template<typename TTransformType>
typename TTransformType::Pointer ReadTransform(const std::string& fileName)
{
    #if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 5) || ITK_VERSION_MAJOR > 4
        typedef itk::TransformFileReaderTemplate<typename TTransformType::ParametersValueType> ReaderType;
    #else
        typedef itk::TransformFileReader ReaderType;
    #endif

    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(fileName);
    reader->Update();

    return dynamic_cast<TTransformType*>(reader->GetTransformList()->front().GetPointer());
}