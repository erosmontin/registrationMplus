#include "itkCommand.h"

template <typename TRegistration>
class RegistrationInterfaceCommand : public itk::Command
{
public:
	typedef  RegistrationInterfaceCommand   Self;
	typedef  itk::Command                   Superclass;
	typedef  itk::SmartPointer<Self>        Pointer;
	itkNewMacro( Self );



protected:
	RegistrationInterfaceCommand() {};

public:
	typedef   TRegistration                              RegistrationType;
	typedef   RegistrationType *                         RegistrationPointer;
	typedef   itk::RegularStepGradientDescentOptimizer   OptimizerType;
	typedef   OptimizerType *                            OptimizerPointer;
	

	void Execute(itk::Object * object, const itk::EventObject & event)
	{
		if( !(itk::IterationEvent().CheckEvent( &event )) )
		{
			return;
		}
		RegistrationPointer registration =
				dynamic_cast<RegistrationPointer>( object );

		OptimizerPointer optimizer = dynamic_cast< OptimizerPointer >(
				registration->GetOptimizer() );


		if ( registration->GetCurrentLevel() == 0 )
		{
			optimizer->SetMaximumStepLength( optimizer->GetMaximumStepLength());
			optimizer->SetMinimumStepLength( optimizer->GetMinimumStepLength());
		}
		else
		{

			optimizer->SetMaximumStepLength( optimizer->GetMaximumStepLength() / maxsteplengthdivisor );
			optimizer->SetMinimumStepLength( optimizer->GetMinimumStepLength() / minimumsteplengthdivisor );
			optimizer->SetNumberOfIterations((int ) optimizer->GetNumberOfIterations() / nidivisor );


		}
		std::cout << "-------------------------------------" << std::endl;
		std::cout << "MultiResolution Level   : "<< registration->GetCurrentLevel()  << "\n";
		std::cout << "Min step length         : "<< optimizer->GetMinimumStepLength()  << "\n";
		std::cout << "Max step length         : "<< optimizer->GetMaximumStepLength()  << "\n";
		std::cout << "Niteration for the step : "<<  optimizer->GetNumberOfIterations() << "\n";
		std::cout<< "Why                      : "<< registration->GetOptimizer()->GetStopConditionDescription()<<"\n";

		std::cout << "-------------------------------------" << std::endl;

	}

	void Execute(const itk::Object * , const itk::EventObject & )
	{ return; }


};



//  The following section of code implements an observer
//  that will monitor the evolution of the registration process.
//
class CommandIterationUpdate : public itk::Command
{
public:
	typedef  CommandIterationUpdate   Self;
	typedef  itk::Command             Superclass;
	typedef  itk::SmartPointer<Self>  Pointer;
	itkNewMacro( Self );

protected:
	CommandIterationUpdate() {};

public:
	typedef   itk::RegularStepGradientDescentOptimizer OptimizerType;
	typedef   const OptimizerType *                    OptimizerPointer;

	void Execute(itk::Object *caller, const itk::EventObject & event)
	{
		Execute( (const itk::Object *)caller, event);
	}

	void Execute(const itk::Object * object, const itk::EventObject & event)
	{
		OptimizerPointer optimizer =
				dynamic_cast< OptimizerPointer >( object );
		if( !(itk::IterationEvent().CheckEvent( &event )) )
		{
			return;
		}
		std::cout << optimizer->GetCurrentIteration() << "   ";
		std::cout << optimizer->GetValue() << "   ";
		std::cout << optimizer->GetCurrentPosition() << std::endl;
	}
};
