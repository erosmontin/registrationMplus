#include <itkImage.h>
#include <itkCommand.h>
#include <itkLBFGSBOptimizer.h>

#include <iostream>

class LBFGSBOptimizeCommandIterationUpdate : public itk::Command
{
public:
    typedef  LBFGSBOptimizeCommandIterationUpdate   Self;
    typedef  itk::Command             Superclass;
    typedef itk::SmartPointer<Self>   Pointer;
    itkNewMacro( Self );
protected:
    LBFGSBOptimizeCommandIterationUpdate() {};
public:
    typedef itk::LBFGSBOptimizer    OptimizerType;
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
        std::cout << optimizer->GetCurrentIteration() << "   ";
        std::cout << optimizer->GetCachedValue() << "   ";
        std::cout << optimizer->GetInfinityNormOfProjectedGradient() << std::endl;
    }
};

#include <itkRegularStepGradientDescentOptimizer.h>
class RegularStepGradientDescentOptimizeCommandIterationUpdate : public itk::Command
{
public:
    typedef  RegularStepGradientDescentOptimizeCommandIterationUpdate   Self;
    typedef  itk::Command             Superclass;
    typedef itk::SmartPointer<Self>   Pointer;
    itkNewMacro( Self );
protected:
    RegularStepGradientDescentOptimizeCommandIterationUpdate() {};
public:
    typedef itk::RegularStepGradientDescentOptimizer    OptimizerType;
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
        std::cout << optimizer->GetCurrentIteration() << "   ";
        std::cout << optimizer->GetValue() << "   ";
    std::cout << optimizer->GetGradientMagnitudeTolerance() << std::endl;
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
