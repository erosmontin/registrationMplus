/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: itkNormalizedGradientFieldImageToImageMetric.txx,v $
Language:  C++
Date:      $Date: 2008-07-03 22:26:16 $
Version:   $Revision: 1.20 $

Copyright (c) 2002 Insight Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkNormalizedGradientFieldImageToImageMetric_txx
#define __itkNormalizedGradientFieldImageToImageMetric_txx

#include "itkNormalizedGradientFieldImageToImageMetric.h"
#include "itkNumericTraits.h"
#include "itkImageRegionConstIteratorWithIndex.h"
//#include "itkVectorBSplineInterpolateImageFunction.h"
#include <itkScaleTransform.h>
#include <iostream>
#include <iomanip>
#include <cstdio>

NS_BEGIN(itk)

/**
 * Constructor
 */
template <class TFixedImage, class TMovingImage> 
NormalizedGradientFieldImageToImageMetric<TFixedImage,TMovingImage>
::NormalizedGradientFieldImageToImageMetric()
{
}


/**
 * Initialize
 */
template <class TFixedImage, class TMovingImage> 
void
NormalizedGradientFieldImageToImageMetric<TFixedImage,TMovingImage>
::Initialize(void) throw ( ExceptionObject )
{

	Superclass::Initialize();
	
	m_MovingNGFEvaluator =  ImageToNGFFilter<TransformedMovingImageType>::New(); 
	
	m_TransformMovingImageFilter = TransformMovingImageFilterType::New();
	
	m_TransformMovingImageFilter->SetTransform(    this->m_Transform );
	m_TransformMovingImageFilter->SetInterpolator( this->m_Interpolator );
	m_TransformMovingImageFilter->SetInput( this->m_MovingImage );
	
	m_TransformMovingImageFilter->SetDefaultPixelValue( 0 );
	
	m_TransformMovingImageFilter->SetSize( this->m_FixedImage->GetLargestPossibleRegion().GetSize() );
	m_TransformMovingImageFilter->SetOutputOrigin( this->m_FixedImage->GetOrigin() );
	m_TransformMovingImageFilter->SetOutputSpacing( this->m_FixedImage->GetSpacing() );
	m_TransformMovingImageFilter->SetOutputDirection( this->m_FixedImage->GetDirection() );
	
	m_MovingNGFEvaluator->SetInput(m_TransformMovingImageFilter->GetOutput()); 
	
	typename  ImageToNGFFilter<FixedImageType>::Pointer ngf_eval = ImageToNGFFilter<FixedImageType>::New(); 
	ngf_eval->SetInput(this->m_FixedImage); 
	ngf_eval->Update(); 
	m_FixedNGF = ngf_eval->GetOutput(); 
	
	m_MovingNGFEvaluator->Update(); 
	m_MovingNGF = m_MovingNGFEvaluator->GetOutput(); 

	for (int i=0; i<TMovingImage::ImageDimension; i++) {
		m_GradientOperators[i].SetDirection( i );
		m_GradientOperators[i].CreateDirectional();
		
		m_GradientFilters[i] = GradientFilter::New();
		
		m_GradientFilters[i]->OverrideBoundaryCondition( &m_MovedBoundCond );
		m_GradientFilters[i]->SetOperator( m_GradientOperators[i] );
		
		m_GradientFilters[i]->SetInput( m_MovingNGF );
		
		m_GradientFilters[i]->UpdateLargestPossibleRegion();
		m_GradientComponent[i] = m_GradientFilters[i]->GetOutput(); 
	}

	if (!m_Evaluator.get())
		m_Evaluator.reset(new NGFScaledDeltaKernel<MovingNGFType,FixedNGFType>); 
}

template <class TFixedImage, class TMovingImage> 
void 
NormalizedGradientFieldImageToImageMetric<TFixedImage,TMovingImage>
::Initialize(const FixedImageType *fixed, 
	     const MovingImageType *moving, 
	     TransformType   *transform, 
	     InterpolatorType *interp, 
	     const RegionType&      fixedRegion)
{
	this->SetFixedImage(fixed); 
	this->SetMovingImage(moving); 
	this->SetTransform(transform); 
	this->SetInterpolator(interp); 
	this->SetFixedImageRegion(fixedRegion); 
	Initialize(); 
}

template <class FI, class MI> 
void NormalizedGradientFieldImageToImageMetric<FI,MI>::SetEvaluator(EvaluatorKernelType *evaluator)
{
	m_Evaluator.reset(evaluator);
}

/**
 * Get the value of the similarity measure
 */
/** Get the derivatives of the match measure. */
template <class FI, class MI> 
void NormalizedGradientFieldImageToImageMetric<FI,MI>::GetDerivative(const TransformParametersType & parameters,
								      DerivativeType  & derivative ) const
{
	typename MovingNGFType::Pointer gradient = GetGradient(parameters); 
	const unsigned int ParametersDimension = this->GetNumberOfParameters();
	derivative = DerivativeType( ParametersDimension );
	
	ImageRegionConstIterator<MovingNGFType>         iti(m_MovingNGF, this->GetFixedImageRegion()); 
	ImageRegionConstIteratorWithIndex<FixedNGFType> ifi(m_FixedNGF, this->GetFixedImageRegion());
	ImageRegionConstIterator<MovingNGFType>         igrad(gradient, this->GetFixedImageRegion());
	
	
	ifi.GoToBegin();
	iti.GoToBegin();
	igrad.GoToBegin();
	
	derivative.Fill( NumericTraits<ITK_TYPENAME DerivativeType::ValueType>::Zero );		
	/*
	  When using a B-Spline transformation, 82% of the registration process is spent in this 
	  loop (according to valgrind). 
	  For this case it should be optimized like in the MattesMutualInformation-Img2ImgMetric, 
	  because most of the time we are adding zeros. 
	*/
	while (! ifi.IsAtEnd()) {
		typename FixedImageType::IndexType index = ifi.GetIndex();
		typename FixedImageType::PointType inputPoint;
		this->m_FixedImage->TransformIndexToPhysicalPoint( index, inputPoint );
		const TransformJacobianType & jacobian = this->m_Transform->GetJacobian( inputPoint ); 
		
		for(unsigned int par=0; par < ParametersDimension; par++)	{
			RealType sum = RealType();
			for(unsigned int dim = 0; dim <  FI::ImageDimension; dim++) {
				sum += jacobian( dim, par ) * igrad.Value()[dim];
			}
			derivative[par] += sum;
		}
		++ifi; 
		++iti; 
		++igrad; 
	}
}

template <class FI, class MI> 
typename NormalizedGradientFieldImageToImageMetric<FI,MI>::MeasureType
NormalizedGradientFieldImageToImageMetric<FI,MI>::GetValue( const TransformParametersType & parameters ) const
{
	this->m_Transform->SetParameters(parameters); 
	m_MovingNGFEvaluator->Update(); 
	return DoGetValue(); 
}

template <class FI, class MI> 
typename NormalizedGradientFieldImageToImageMetric<FI,MI>::MeasureType
NormalizedGradientFieldImageToImageMetric<FI,MI>::DoGetValue(  ) const
{
	ImageRegionConstIterator<MovingNGFType> iti(m_MovingNGF, this->GetFixedImageRegion()); 
	ImageRegionConstIterator<FixedNGFType> ifi(m_FixedNGF, this->GetFixedImageRegion());
	
	MeasureType value = MeasureType();
	
	while (!iti.IsAtEnd())  {
		value += m_Evaluator->Value(iti.Value(), ifi.Value()); 
		++iti; 
		++ifi; 
	}
	return value / this->GetFixedImageRegion().GetNumberOfPixels(); 
}

template <class FI, class MI> 
void	
NormalizedGradientFieldImageToImageMetric<FI,MI>::GetValueAndDerivative( const TransformParametersType & parameters, 
									  MeasureType& Value, DerivativeType& derivative ) const
{
	// this is certainly not the most optimal version, because 
	// some evaluations are run twice,
	// but GetDerivative is already messy enough, and duplicating it's code 
	// here is really not a nice idea. 
	GetDerivative(parameters, derivative); 
	Value = DoGetValue(); 
}

template <class FI, class MI> 
typename NormalizedGradientFieldImageToImageMetric<FI,MI>::MovingNGFType::Pointer
NormalizedGradientFieldImageToImageMetric<FI,MI>::GetGradient(const TransformParametersType & parameters) const
{
	// transform moving image 
	this->m_Transform->SetParameters(parameters); 
	m_MovingNGFEvaluator->Update(); 

	ImageRegionConstIterator<MovingNGFType> iti(m_MovingNGF, this->GetFixedImageRegion()); 
	ImageRegionConstIterator<FixedNGFType> ifi(m_FixedNGF, this->GetFixedImageRegion());

	typename MovingNGFType::Pointer gradient = MovingNGFType::New(); 
	gradient->SetRegions(this->GetFixedImageRegion().GetSize()); 
	gradient->Allocate(); 
	
	ImageRegionIterator<MovingNGFType> iout( gradient, this->GetFixedImageRegion());
	
	ImageRegionConstIterator<MovingNGFType> igrad[MI::ImageDimension]; 
	
	// evaluate the gradient NGF 
	for (size_t i = 0; i < MI::ImageDimension; ++i) {
		m_GradientFilters[i]->Update(); 
		igrad[i] = ImageRegionConstIterator<MovingNGFType>(m_GradientComponent[i], 
								   this->GetFixedImageRegion()); 
	}

	while (!iti.IsAtEnd())  {
		
		m_Evaluator->Gradient(iout.Value(),iti.Value(),igrad,ifi.Value()); 
		
		// evaluating the gradinet manually shows that for some 
		// reason the sign is wrong - maybe the DerivativeNeighborhoodOperator 
		// returns the values with an unexpected sign?
		
		iout.Value() *= -1; 

		++iti; 
		++ifi; 
		++iout; 
		for (size_t i = 0; i < MI::ImageDimension; ++i) {
			++igrad[i]; 
		}
	}
	return gradient; 
}

NS_END(itk)


#endif
