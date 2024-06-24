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


#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API

#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "itkNormalizedGradientFieldImageToImageMetric.h"
#include "itkBSplineDeformableTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkCenteredAffineTransform.h"

using namespace std; 
using namespace itk; 

template <class MovingGF, class FixedGF> 
class MockNGFKernel: public NGFKernel<MovingGF,FixedGF> {
public: 
	typedef NGFKernel<MovingGF,FixedGF> Superclass; 
	typedef typename Superclass::MV MV; 
	typedef typename Superclass::FV FV; 
	typedef typename Superclass::Real Real; 
	typedef typename Superclass::GradientIteratorType GradientIteratorType; 
	
	virtual typename NGFKernel<MovingGF,FixedGF>::Real
	Value(const MV& moving, const FV& fixed)const {
		return 1.0; 
	}
	
	virtual typename NGFKernel<MovingGF,FixedGF>::Real
		Gradient(MV& result, 
			 const MV& moving,
			 const GradientIteratorType gradmoving[], 
			 const FV& fixed)const 
	{
		for (unsigned int i = 0; i < MovingGF::ImageDimension; ++i) 
			result[i] = i + 1.0;
		return 1.0;
	}
}; 


template <int dim>
struct ImageCreator {
	typedef Image<float, dim> CImage; 
	typedef typename CImage::Pointer PImage; 
	
	static PImage create_image(typename CImage::SizeType size, const float *data); 
}; 

template <int dim, class Transform>
struct NGFMetricFixture: public ImageCreator<dim> {

	typedef typename ImageCreator<dim>::CImage CImage;
	typedef typename ImageCreator<dim>::PImage PImage; 
	typedef  ImageToNGFFilter<CImage> CImageToNGF; 
	typedef typename CImageToNGF::Pointer PImageToNGF;
	typedef NormalizedGradientFieldImageToImageMetric< CImage, CImage> CNGFMetric; 
	typedef typename CNGFMetric::Pointer PNGFMetric; 
	typedef Transform CTransform; 
	typedef typename CTransform::Pointer PTransform; 
	typedef NearestNeighborInterpolateImageFunction< CImage,double > CInterpolator; 
	typedef typename CInterpolator::Pointer PInterpolator; 

	PImage src; 
	PImage ref; 
	PNGFMetric metric;
	PTransform    transform; 
	PInterpolator interpolator; 

	void Initialize(typename CImage::SizeType size,  const float *src_data, const float *ref_data); 

	virtual void InitializeTransform() = 0; 
}; 

template <int dim>
struct CenteredAffineNGFMetricFixture: public NGFMetricFixture<dim, CenteredAffineTransform<double, dim> > {

	typedef NGFMetricFixture<dim, CenteredAffineTransform<double, dim> > Superclass; 
	typedef typename Superclass::CImage CImage;
	typedef typename Superclass::PImage PImage; 
	typedef ImageToNGFFilter<CImage> CImageToNGF; 
	typedef typename CImageToNGF::Pointer PImageToNGF;
	typedef typename Superclass::CNGFMetric CNGFMetric; 
	typedef typename CNGFMetric::Pointer PNGFMetric; 
	typedef typename Superclass::CTransform CTransform; 
	typedef typename CTransform::Pointer PTransform; 
	typedef typename Superclass::CInterpolator CInterpolator; 
	typedef typename CInterpolator::Pointer PInterpolator; 

	virtual void InitializeTransform(); 
}; 


// Test evaluation of the metric and the gradient 
BOOST_FIXTURE_TEST_CASE( test_metric_ngf_2d, CenteredAffineNGFMetricFixture<2>) 
{
	CImage::SizeType size = {{5, 6}}; 

	const float ref_data[5 * 6] =  {
		0, 0, 0, 0, 0, 
		0, 0, 0, 0, 0, 
		0, 0, 4, 0, 0, 
		0, 4, 0, 0, 0, 
		0, 0, 0, 4, 0, 
		0, 0, 0, 0, 0
	}; 	

	const float src_data[5 * 6] =  {
		0, 0, 0, 0, 0, 
		0, 0, 4, 0, 0, 
		0, 4, 0, 0, 0, 
		0, 0, 0, 4, 0, 
		0, 0, 0, 0, 0, 
		0, 0, 0, 0, 0
	}; 	

	Initialize(size,  src_data, ref_data); 

	CTransform::ParametersType params = transform->GetParameters(); 
	
	BOOST_CHECK_CLOSE(metric->GetValue(params), 1.0, 0.1); 

	CNGFMetric::MovingNGFType::Pointer g = metric->GetGradient(params); 
	CImage::RegionType r = g->GetLargestPossibleRegion(); 

	ImageRegionConstIterator<CNGFMetric::MovingNGFType> it(g, r); 

	while (! it.IsAtEnd()) {
		BOOST_CHECK_CLOSE( it.Value()[0], -1.0f, 0.1f); 
		BOOST_CHECK_CLOSE( it.Value()[1], -2.0f, 0.1f); 
		++it; 
	}
}


/*
  Implementations of the fixture. 
 */

template <int dim>
typename ImageCreator<dim>::PImage 
ImageCreator<dim>::create_image(typename CImage::SizeType size, const float *data)
{
	PImage image = CImage::New(); 
	
	image->SetRegions(size); 
	image->Allocate(); 

	typename CImage::RegionType region = image->GetLargestPossibleRegion (); 
	
	ImageRegionIterator<CImage> it(image, region); 

	while (!it.IsAtEnd()) {
		it.Set(*data++); 
		++it; 
	}
	return image;
}


template <int dim, class Transform>
void NGFMetricFixture<dim, Transform >::Initialize(typename CImage::SizeType size,  
				       const float *src_data, const float *ref_data)
{
	src = create_image(size, src_data); 
	ref = create_image(size, ref_data); 
	interpolator = CInterpolator::New(); 
	transform = CTransform::New(); 

	metric = CNGFMetric::New(); 
	metric->SetEvaluator(new MockNGFKernel<typename CNGFMetric::MovingNGFType, 
			     typename CNGFMetric::FixedNGFType>()); 
	InitializeTransform(); 
	metric->Initialize(ref, src, transform, interpolator, ref->GetLargestPossibleRegion()); 

}

template <int dim>
void CenteredAffineNGFMetricFixture<dim>::InitializeTransform()
{
	
}


