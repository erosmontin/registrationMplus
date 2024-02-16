#ifndef __itkMANGF2_hxx
#define __itkMANGF2_hxx

#include "itkMANGF2.h"
#include "Code/itkNGFMetricKernel.h"

namespace itk
{
/**
 * Constructor
 */
template <class TFixedImage, class TMovingImage>
MANGF2<TFixedImage, TMovingImage>
::MANGF2()
 {
	m_MA=MattesType::New();
	m_NGF=NGFType::New();
	m_INTERNALL_interpolator=LFType::New();
	m_lambda=1.0;
	m_lambdaDerivative=m_lambda;
	m_NGFnumberOfSamples=20000;
	m_MAnumberOfSamples=20000;
	m_BinNumbers=50;
	m_NumberOfThreads=1;
	m_Evaluator=0;
	m_FixedEta=1;
	m_MovingEta=1;



 }

template <class TFixedImage, class TMovingImage>
MANGF2<TFixedImage, TMovingImage>
::~MANGF2(){}

/**
 * Initialize
 */


template <class TFixedImage, class TMovingImage> 
void
MANGF2<TFixedImage,TMovingImage>
::Initialize(void) throw ( ExceptionObject )
 {

	Superclass::Initialize();


	m_MA->SetFixedImage(this->m_FixedImage);
	m_MA->SetMovingImage(this->m_MovingImage);
	m_MA->SetInterpolator(this->m_Interpolator);
	m_MA->SetTransform(this->m_Transform);
	m_MA->SetFixedImageRegion(this->m_FixedImage->GetBufferedRegion() );
	m_MA->SetNumberOfHistogramBins(this->m_BinNumbers);
	m_MA->SetNumberOfSpatialSamples(this->m_MAnumberOfSamples);
	m_MA->SetNumberOfThreads(this->m_NumberOfThreads);
    m_MA->ReinitializeSeed();
	m_MA->Initialize();




	m_NGF->SetFixedImage(this->m_FixedImage);
	m_NGF->SetMovingImage(this->m_MovingImage);
	m_NGF->SetInterpolator(this->m_INTERNALL_interpolator);
	m_NGF->SetTransform(this->m_Transform);
	m_NGF->SetNumberOfSpatialSamples(this->m_NGFnumberOfSamples);
	m_NGF->SetNumberOfThreads(this->m_NumberOfThreads);
	RegionType inputRegion = this->m_FixedImage->GetLargestPossibleRegion();
	m_NGF->SetFixedNoise(this->m_FixedEta);//this NGF noise can be considered the eta parameter
	m_NGF->SetMovingNoise(this->m_MovingEta);
	m_NGF->SetFixedImageRegion(inputRegion);


switch(m_Evaluator)
{
case (0):
//
// Implementation of scalar product based evaluator 
//
m_NGF->SetEvaluator(new NGFScalarKernel<MovingNGFType,FixedNGFType>());
break;
//
// Implementation of cross product based evaluator 
//
case (1):
m_NGF->SetEvaluator(new NGFCrossKernel<MovingNGFType,FixedNGFType>());
break;
//
// Implementation of scaled difference based evaluator 
//
case (2):
m_NGF->SetEvaluator(new NGFScaledDeltaKernel<MovingNGFType,FixedNGFType>());
break;
//
// Implementation of the delta cost evalator 
// Considers gradients only similar if they point in the same direction 
// 
case (3):
m_NGF->SetEvaluator(new NGFDeltaKernel<MovingNGFType,FixedNGFType>());
break;
//
// Implementation of the squared delta cost evalator 
// 
case (4):
m_NGF->SetEvaluator(new NGFDelta2Kernel<MovingNGFType,FixedNGFType>());
break;
default:
m_NGF->SetEvaluator(new NGFScalarKernel<MovingNGFType,FixedNGFType>());
}

	
	m_NGF->Initialize();





 }


template <class TFixedImage, class TMovingImage>
typename MANGF2<TFixedImage, TMovingImage>::MeasureType
MANGF2<TFixedImage, TMovingImage>
::GetValue(const ParametersType & parameters) const
 {
	double a,b;

	a=this->GetMAValue(parameters);
	b=this->GetNGFValue(parameters);


	return a * this->m_lambda * b ;

 }

template <class TFixedImage, class TMovingImage>
typename MANGF2<TFixedImage, TMovingImage>::MeasureType
MANGF2<TFixedImage, TMovingImage>
::GetNGFValue(const ParametersType & parameters) const
 {
	return static_cast<MeasureType> ( m_NGF->GetValue(parameters) );
 }

template <class TFixedImage, class TMovingImage>
typename MANGF2<TFixedImage, TMovingImage>::MeasureType
MANGF2<TFixedImage, TMovingImage>
::GetMAValue(const ParametersType & parameters) const
 {
	return static_cast<MeasureType> ( m_MA->GetValue(parameters) );

 }


template <class TFixedImage, class TMovingImage>
void
MANGF2<TFixedImage, TMovingImage>
::GetDerivative(const ParametersType & parameters,DerivativeType & derivative) const
 {
	DerivativeType a;
	this->m_MA->GetDerivative(parameters,a);

	DerivativeType b;
	this->m_NGF->GetDerivative(parameters,b);


	derivative=a;
	for(long unsigned int p=0;p<derivative.GetSize();p++)
	{derivative[p]= (a[p]* this->GetNGFValue(parameters))+ this->m_lambdaDerivative*b[p]*this->GetMAValue(parameters);}
		

	


 }


template <class TFixedImage, class TMovingImage>
void
MANGF2<TFixedImage, TMovingImage>
::GetMADerivative(const ParametersType & parameters,DerivativeType & derivative) const
 {
	
	m_MA->GetDerivative(parameters,derivative);	


 }


template <class TFixedImage, class TMovingImage>
void
MANGF2<TFixedImage, TMovingImage>
::GetNGFDerivative(const ParametersType & parameters,DerivativeType & derivative) const
 {
	
	m_NGF->GetDerivative(parameters,derivative);	


 }


template <class TFixedImage, class TMovingImage>
void
MANGF2<TFixedImage, TMovingImage>
::GetValueAndDerivative(const ParametersType & parameters, MeasureType & Value, DerivativeType & Derivative) const
 {
	Value=this->GetValue(parameters);
	this->GetDerivative(parameters,Derivative);
 }

template <class TImageType, class TMovingImage>
void
MANGF2<TImageType, TMovingImage>::
PrintSelf( std::ostream& os, Indent indent ) const
{
	Superclass::PrintSelf(os,indent);

}


} // end namespace itk

#endif
