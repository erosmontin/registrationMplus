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
{m_FixedNoise=1;
	m_MovingNoise=1;
}


/**
 * Initialize
 */
template <class TFixedImage, class TMovingImage> 
void
NormalizedGradientFieldImageToImageMetric<TFixedImage,TMovingImage>
::Initialize(void)
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
	m_MovingNGFEvaluator->SetNoise( m_MovingNoise );
	typename  ImageToNGFFilter<FixedImageType>::Pointer ngf_eval = ImageToNGFFilter<FixedImageType>::New(); 
	ngf_eval->SetInput(this->m_FixedImage); 
	ngf_eval->SetNoise( m_FixedNoise );
	ngf_eval->Update(); 
	m_FixedNGF = ngf_eval->GetOutput(); 
	
	m_MovingNGFEvaluator->Update(); 
	m_MovingNGF = m_MovingNGFEvaluator->GetOutput(); 

	for (unsigned int i=0; i<TMovingImage::ImageDimension; i++) {
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

    m_CachedParameters = this->m_Transform->GetParameters();
    m_CachedGradient = this->GetGradient(m_CachedParameters);
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
void NormalizedGradientFieldImageToImageMetric<FI,MI>::GetDerivative(
    const TransformParametersType & parameters,
    DerivativeType  & derivative ) const
{
    if (!m_CachedGradient || parameters != m_CachedParameters) {
        m_CachedGradient = this->GetGradient(parameters);
        m_CachedParameters = parameters;
    }
    typename MovingNGFType::Pointer gradient = m_CachedGradient;
    const unsigned int ParametersDimension = this->GetNumberOfParameters();
    derivative = DerivativeType( ParametersDimension );
    derivative.Fill( NumericTraits<ITK_TYPENAME DerivativeType::ValueType>::Zero );

    // Get region and size
    const RegionType& region = this->GetFixedImageRegion();
    const size_t numPixels = region.GetNumberOfPixels();

    // Use raw buffer pointers for cache-friendly access
    const auto* gradientBuffer = gradient->GetBufferPointer();

    // Precompute all indices for the region (cache-friendly, sequential access)
    std::vector<typename FixedImageType::IndexType> indices;
    indices.reserve(numPixels); // Reserve memory to avoid reallocations
    itk::ImageRegionConstIteratorWithIndex<FixedNGFType> ifi(m_FixedNGF, region);
    for (size_t idx = 0; idx < numPixels; ++idx, ++ifi) {
        indices.push_back(ifi.GetIndex());
    }

    // Main loop: cache Jacobian and gradient vector per pixel
    int numThreads = 1;
    #ifdef _OPENMP
    numThreads = omp_get_max_threads();
    #endif
    std::vector<DerivativeType> threadDerivatives(numThreads, DerivativeType(ParametersDimension));
    for (int t = 0; t < numThreads; ++t)
        threadDerivatives[t].Fill(NumericTraits<typename DerivativeType::ValueType>::Zero);

    #pragma omp parallel
    {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif

        #pragma omp for
        for (size_t idx = 0; idx < numPixels; ++idx) {
            const auto& index = indices[idx];
            typename FixedImageType::PointType inputPoint;
            this->m_FixedImage->TransformIndexToPhysicalPoint(index, inputPoint);
            const TransformJacobianType& jacobian = this->m_Transform->GetJacobian(inputPoint);
            const auto& gradVec = gradientBuffer[idx];

            for (unsigned int par = 0; par < ParametersDimension; ++par) {
                bool allZero = true;
                for (unsigned int dim = 0; dim < FI::ImageDimension; ++dim) {
                    if (jacobian(dim, par) != 0.0) {
                        allZero = false;
                        break;
                    }
                }
                if (allZero) continue;

                RealType sum = NumericTraits<RealType>::ZeroValue();
                for (unsigned int dim = 0; dim < FI::ImageDimension; ++dim) {
                    auto jacobianValue = jacobian(dim, par);
                    auto gradientValue = gradVec[dim];
                    if (jacobianValue != 0.0 && gradientValue != 0.0) {
                        sum += jacobianValue * gradientValue;
                    }
                }
                threadDerivatives[tid][par] += sum;
            }
        }
    }

    // Reduce
    for (int t = 0; t < numThreads; ++t)
        derivative += threadDerivatives[t];
}

template <class FI, class MI> 
typename NormalizedGradientFieldImageToImageMetric<FI,MI>::MeasureType
NormalizedGradientFieldImageToImageMetric<FI,MI>::GetValue( const TransformParametersType & parameters ) const
{
    // Only update if parameters have changed
    if (m_CachedValueParameters != parameters) {
        this->m_Transform->SetParameters(parameters); 
        m_MovingNGFEvaluator->Update(); 
        m_CachedValueParameters = parameters;
    }
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
		
		
		// Eros Montin update 09/23/25
		//The sign flip is needed because the ITK derivative operator (or your kernel) returns the gradient with a sign opposite to what your metric expects. This is a common issue when using finite difference operators.
// If possible, fix the sign at the source (operator or kernel). If not, keep the flip and document it.
		
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
