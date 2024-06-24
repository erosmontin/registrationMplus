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

#ifndef itkNGFMetricKernel_h
#define itkNGFMetricKernel_h

#include <itkCovariantVector.h>
#include <itkImage.h>
#include <itkImageRegionConstIterator.h>
#include "itkMyDefines.h"

#ifdef _MSC_VER
#if _MSC_VER < 1310 
#warning You are using Visual C++ < 7.1, consider updating to a compiler that supports partial template specialisation
#endif
#endif

NS_BEGIN(itk) 


/**
   This iterator trait is added to be able to mock an iterator 
   in the test code. 
 */

template <typename T, unsigned int dim> 
struct  NGFIteratorTrait {
	typedef void ConstIterator; 
}; 

/**
   Template specialization for the "normal" use case of the iterator. 
 */
template <typename T, unsigned int dim> 
struct  NGFIteratorTrait<Image<T, dim>, dim > {
	typedef ImageRegionConstIterator<Image<T, dim> > ConstIterator; 
}; 


/**
   Abstract base class of the NGF kernel evaluators. 
 */
template <class MovingGF, class FixedGF> 
class NGFKernel {
public: 
	typedef typename MovingGF::PixelType MV; 
	typedef typename FixedGF::PixelType  FV; 
	typedef typename MovingGF::PixelType::RealValueType Real; 
	typedef typename NGFIteratorTrait<FixedGF, FixedGF::ImageDimension>::ConstIterator GradientIteratorType; // changed from MovingGF to FixedGF to allow for different image types EROS.	
	typedef std::unique_ptr<NGFKernel> Pointer;
	// typedef std::auto_ptr<NGFKernel> Pointer; 
	
	/**
	   Evalute the cost function value based on the given NGF vectors.
	   \param moving Gradient field value of the moving image 
	   \param fixed  Gradient field value of the fixed image 
	   \remark: Abstract function that needs to be specified in derived classes. 
	*/
	virtual typename NGFKernel<MovingGF,FixedGF>::Real
		Value(const MV& moving, const FV& fixed) const = 0;  


	/**
	   Evalute the cost function value and its gradient based on the given NGF vectors.
	   \retval result Cost function Gradient at the requested position 
	   \param moving Gradient field value of the moving image 
	   \param gradmoving Gradient of the moving image at the requested position 
	   \param fixed  Gradient field value of the fixed image 
	   \remark: Abstract function that needs to be specified in derived classes. 
	 */
	virtual typename NGFKernel<MovingGF,FixedGF>::Real
		Gradient(MV& result,
			 const MV& moving,
			 const GradientIteratorType gradmoving[],
			 const FV& fixed) const = 0;
};


/**
   This NGF evaluator kernel implements the scalar product based version of the 
   NGF metric as proposed in Haber et.al. 
*/
template <class MovingGF, class FixedGF> 
class NGFScalarKernel: public NGFKernel<MovingGF,FixedGF> {
public: 
	typedef NGFKernel<MovingGF,FixedGF> Superclass; 
	typedef typename Superclass::MV MV; 
	typedef typename Superclass::FV FV; 
	typedef typename Superclass::Real Real; 
	typedef typename Superclass::GradientIteratorType GradientIteratorType; 
	
	virtual typename NGFKernel<MovingGF,FixedGF>::Real
		Value(const MV& moving, const FV& fixed)const;  
	virtual typename NGFKernel<MovingGF,FixedGF>::Real
		Gradient(MV& result, 
			   const MV& moving,
			   const GradientIteratorType gradmoving[], 
			   const FV& fixed)const ; 
}; 

/**
   This NGF evaluator kernel implements the cross product based version of the 
   NGF metric as proposed in Haber et.al. 
   It is only defined for images of two and three dimensions. 
*/
template <class MovingGF, class FixedGF> 
class NGFCrossKernel: public NGFKernel<MovingGF,FixedGF> {
public: 
	typedef NGFKernel<MovingGF,FixedGF> Superclass; 
	typedef typename Superclass::MV MV; 
	typedef typename Superclass::FV FV; 
	typedef typename Superclass::Real Real; 
	typedef typename Superclass::GradientIteratorType GradientIteratorType; 

	virtual typename NGFKernel<MovingGF,FixedGF>::Real
		Value(const MV& moving, const FV& fixed) const;  
	virtual typename NGFKernel<MovingGF,FixedGF>::Real
		Gradient(MV& result, 
			   const MV& moving,
			   const GradientIteratorType gradmoving[], 
			   const FV& fixed)const ; 
}; 


/**
   This NGF evaluator kernel implements the scaled scalar product based version of the 
   NGF metric as proposed in Wollny et al. (submitted)
*/
template <class MovingGF, class FixedGF> 
class NGFScaledDeltaKernel: public NGFKernel<MovingGF,FixedGF> {
public: 
	typedef NGFKernel<MovingGF,FixedGF> Superclass; 
	typedef typename Superclass::MV MV; 
	typedef typename Superclass::FV FV; 
	typedef typename Superclass::Real Real; 
	typedef typename Superclass::GradientIteratorType GradientIteratorType; 

	virtual typename NGFKernel<MovingGF,FixedGF>::Real
		Value(const MV& moving, const FV& fixed)const;  
	virtual typename NGFKernel<MovingGF,FixedGF>::Real
		Gradient(MV& result, 
			 const MV& moving,
			 const GradientIteratorType gradmoving[], 
			 const FV& fixed)const ; 
}; 

/**
   This NGF evaluator kernel implements another version of the metric 
   Only to be used for testing purpouses. 
*/

template <class MovingGF, class FixedGF> 
class NGFDeltaKernel: public NGFKernel<MovingGF,FixedGF> {
public: 
	typedef NGFKernel<MovingGF,FixedGF> Superclass; 
	typedef typename Superclass::MV MV; 
	typedef typename Superclass::FV FV; 
	typedef typename Superclass::Real Real; 
	typedef typename Superclass::GradientIteratorType GradientIteratorType; 

	virtual typename NGFKernel<MovingGF,FixedGF>::Real
		Value(const MV& moving, const FV& fixed)const;  
	virtual typename NGFKernel<MovingGF,FixedGF>::Real
		Gradient(MV& result, 
			 const MV& moving,
			 const GradientIteratorType gradmoving[], 
			 const FV& fixed)const ; 
}; 

/**
   This NGF evaluator kernel implements another version of the metric 
   Only to be used for testing purpouses. 
*/
template <class MovingGF, class FixedGF> 
class NGFDelta2Kernel: public NGFKernel<MovingGF,FixedGF> {
public: 
	typedef NGFKernel<MovingGF,FixedGF> Superclass; 
	typedef typename Superclass::MV MV; 
	typedef typename Superclass::FV FV; 
	typedef typename Superclass::Real Real; 
	typedef typename Superclass::GradientIteratorType GradientIteratorType; 

	virtual typename NGFKernel<MovingGF,FixedGF>::Real
		Value(const MV& moving, const FV& fixed)const;  
	virtual typename NGFKernel<MovingGF,FixedGF>::Real
		Gradient(MV& result, 
			 const MV& moving,
			 const GradientIteratorType gradmoving[], 
			 const FV& fixed)const ; 
}; 


NS_END(itk) 

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkNGFMetricKernel.txx"
#endif

#endif
