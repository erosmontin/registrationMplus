/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkNormalizedGradientFieldsImageToImageMetric.h,v $
  Language:  C++
  Date:      $Date: 2008-06-29 01:24:10 $
  Version:   $Revision: 1.9 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkNormalizedGradientFieldsImageToImageMetric_h
#define __itkNormalizedGradientFieldsImageToImageMetric_h

#include <itkImageToImageMetric.h>
#include <itkVectorResampleImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkSobelOperator.h>
#include <itkDerivativeOperator.h>
#include <itkVectorNeighborhoodOperatorImageFilter.h>
#include <itkZeroFluxNeumannBoundaryCondition.h>

#include "itkMyDefines.h"
#include "itkImageToNGFFilter.h" 
#include "itkNGFMetricKernel.h" 

NS_BEGIN(itk)

/** 
 * \class NormalizedGradientFieldsImageToImageMetric
 * \brief Computes similarity between two objects to be registered base on 
 * Normalized Gradient Fields.  
 *
 * E. Haber and J. Modersitzki, "Beyond Mutual Information: A Simple and Robust
 * Alternative", in Bildverarbeitung für die Medizin 2005, eds. Meinzer, Handels, 
 * Horsch and Tolxdorff, 2005, pp. 350-354
 * and 
 * G. Wollny, M. J. Ledesma-Carbayo, P. Kellman, and A. Santos, "Exploiting 
 * Periodicity in Motion Compensation of Free-Breathing Myocardial Perfusion MRI", 
 * submitted to IEEE Trans. Med. Imag.   
 *
 * This Class is templated over the type of the images to be compared 
 * The code requires that both images create the same type of vector. 
 * \todo TheMetric doesn't support masking yet
 *
 * \ingroup RegistrationMetrics
 *
 */

template < class TFixedImage, class TMovingImage >
class ITK_EXPORT NormalizedGradientFieldImageToImageMetric :
public ImageToImageMetric< TFixedImage, TMovingImage>
{
public:

  /** Standard class typedefs. */
  typedef NormalizedGradientFieldImageToImageMetric           Self;
  typedef ImageToImageMetric<TFixedImage, TMovingImage > Superclass;

  typedef SmartPointer<Self>         Pointer;
  typedef SmartPointer<const Self>   ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(NormalizedGradientFieldImageToImageMetric, ImageToImageMetric);


  /** Types transferred from the base class */
// Work around a Visual Studio .NET bug
#if defined(_MSC_VER) && (_MSC_VER == 1300)
  typedef double RealType;
#else
  typedef typename Superclass::RealType                 RealType;
#endif
  typedef typename Superclass::TransformType            TransformType;
  typedef typename Superclass::TransformPointer         TransformPointer;
  typedef typename Superclass::TransformParametersType  TransformParametersType;
  typedef typename Superclass::TransformJacobianType    TransformJacobianType;

  typedef typename Superclass::MeasureType              MeasureType;
  typedef typename Superclass::DerivativeType           DerivativeType;
  typedef typename Superclass::FixedImageType           FixedImageType;
  typedef typename Superclass::MovingImageType          MovingImageType;
  typedef typename Superclass::FixedImageConstPointer   FixedImageConstPointer;
  typedef typename Superclass::MovingImageConstPointer  MovingImageConstPointer;
  typedef typename Superclass::InterpolatorType         InterpolatorType;
  typedef typename FixedImageType::RegionType           RegionType; 


  typedef typename TFixedImage::PixelType               FixedImagePixelType;
  typedef typename TMovingImage::PixelType              MovedImagePixelType;
  typedef typename TMovingImage::RegionType::SizeType   SizeType;

  itkStaticConstMacro(FixedImageDimension, unsigned int, TFixedImage::ImageDimension);

  /** Types for transforming the moving image */
  typedef Image< FixedImagePixelType,
	  itkGetStaticConstMacro( FixedImageDimension ) >
	  TransformedMovingImageType;

  typedef typename ImageToNGFFilter<TransformedMovingImageType>::OutputImageType  MovingNGFType; 
  typedef typename ImageToNGFFilter<FixedImageType>::OutputImageType              FixedNGFType; 

  typedef NGFKernel<MovingNGFType,FixedNGFType> EvaluatorKernelType; 

  /** Get the derivatives of the match measure. */
  void GetDerivative( const TransformParametersType& parameters, DerivativeType& derivative ) const;

  /**  Get the value for single valued optimizers. */
  MeasureType GetValue( const TransformParametersType& parameters ) const;

  /**  Get value and derivatives for multiple valued optimizers. */
  void GetValueAndDerivative( const TransformParametersType & parameters,
                              MeasureType& Value, DerivativeType& derivative ) const;

  /** Initialize the Metric by making sure that all the components
   *  are present and plugged together correctly     */
  virtual void Initialize(void);

  /** Initialize the Metric by setting all the required components 
      \param fixed 
      \param moving
      \param transform
      \param interp
      \param fixedRegion
   */
  virtual void Initialize(const FixedImageType *fixed, 
			  const MovingImageType *moving, 
			  TransformType   *transform, 
			  InterpolatorType   *interp, 
			  const RegionType&         fixedRegion); 
  
  /**
    The the gradint of the cost function, public for testing purpouses
  */
  typename MovingNGFType::Pointer GetGradient(const TransformParametersType & parameters) const; 

  /**
     Sets the actual kernel evaluator to something else then the standard \a NGFScaledDeltaKernel
     \param evaluator
   */

  void SetEvaluator(EvaluatorKernelType *evaluator); 

	itkGetMacro( FixedNoise, double);
	itkSetMacro( FixedNoise, double);

	itkGetMacro( MovingNoise, double);
	itkSetMacro( MovingNoise, double);


protected:
  NormalizedGradientFieldImageToImageMetric();
  virtual ~NormalizedGradientFieldImageToImageMetric() {};

private:
  
  NormalizedGradientFieldImageToImageMetric(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  // evaluate the cost value 
  MeasureType DoGetValue() const; 

  typedef ResampleImageFilter< MovingImageType, TransformedMovingImageType >
	  TransformMovingImageFilterType;

  typename TransformMovingImageFilterType::Pointer m_TransformMovingImageFilter;
  
  typename ImageToNGFFilter<TransformedMovingImageType>::Pointer m_MovingNGFEvaluator; 
  typename MovingNGFType::Pointer m_MovingNGF; 
  typename FixedNGFType::Pointer m_FixedNGF; 

  typename EvaluatorKernelType::Pointer m_Evaluator; 
  typedef VectorNeighborhoodOperatorImageFilter<MovingNGFType, MovingNGFType > GradientFilter;

  ZeroFluxNeumannBoundaryCondition<MovingNGFType>            m_MovedBoundCond;
  DerivativeOperator<float, TMovingImage::ImageDimension>    m_GradientOperators[TMovingImage::ImageDimension];
  
  typename GradientFilter::Pointer            m_GradientFilters[TMovingImage::ImageDimension]; 
  typename GradientFilter::OutputImagePointer m_GradientComponent[TMovingImage::ImageDimension];

	double m_FixedNoise;
	double m_MovingNoise;
  
  mutable typename MovingNGFType::Pointer m_CachedGradient;
  mutable TransformParametersType m_CachedParameters;
  mutable TransformParametersType m_CachedValueParameters;

  
};

NS_END(itk)

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkNormalizedGradientFieldImageToImageMetric.txx"
#endif

#endif
