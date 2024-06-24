/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

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
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <itkImageRegionConstIterator.h>
#include <itkCenteredAffineTransform.h>
#include <itkNearestNeighborInterpolateImageFunction.h>

#include "itkNormalizedGradientFieldImageToImageMetric.h"


using namespace itk; 

template <int dim>
struct ImageCreator {
	typedef Image<float, dim> CImage; 
	typedef typename CImage::Pointer PImage; 
	
	static PImage create_image(typename CImage::SizeType size, const float *data); 
}; 

template <int dim>
struct Image2NGFFixture : public ImageCreator<dim> {
	
	typedef typename ImageCreator<dim>::CImage CImage;
	typedef typename ImageCreator<dim>::PImage PImage; 
	typedef  ImageToNGFFilter<CImage> CImageToNGF; 
	typedef typename CImageToNGF::Pointer PImageToNGF;

	Image2NGFFixture(); 
	void Initialize(typename CImage::SizeType size,  const float *data); 

	
	PImage src; 
	PImageToNGF filter; 


}; 

template <int dim>
struct NGFMetricFixture: public ImageCreator<dim> {

	typedef typename ImageCreator<dim>::CImage CImage;
	typedef typename ImageCreator<dim>::PImage PImage; 
	typedef  ImageToNGFFilter<CImage> CImageToNGF; 
	typedef typename CImageToNGF::Pointer PImageToNGF;
	typedef NormalizedGradientFieldImageToImageMetric< CImage, CImage> CNGFMetric; 
	typedef typename CNGFMetric::Pointer PNGFMetric; 
	typedef CenteredAffineTransform<double, dim> CTransform; 
	typedef typename CTransform::Pointer PTransform; 
	typedef NearestNeighborInterpolateImageFunction< CImage,double > CInterpolator; 
	typedef typename CInterpolator::Pointer PInterpolator; 

	PImage src; 
	PImage ref; 
	PNGFMetric metric;
	PTransform    transform; 
	PInterpolator interpolator; 

	void Initialize(typename CImage::SizeType size,  const float *src_data, const float *ref_data); 
}; 


// Test the evaluation of a normalized gradient field from an image
BOOST_FIXTURE_TEST_CASE( test_image_to_ngf_2d, Image2NGFFixture<2>) 
{
	Image2NGFFixture<2>::CImage::SizeType size = {{5, 6}}; 
	
	const float src_data[5 * 6] =  {
		0, 0, 0, 0, 0, 
		0, 0, 0, 0, 0, 
		0, 0, 4, 0, 0, 
		0, 4, 0, 0, 0, 
		0, 0, 0, 4, 0, 
		0, 0, 0, 0, 0
	}; 


	
	const float test_data_x[5 * 6] =  {
		0,       0,       0,       0,        0, 
		0,       0,       0,       0,        0, 
		0,       0.66941,   0,    -0.90108,  0,
		0.90108, 0,      -0.66941, 0,        0, 
		0,       0,       0.90108, 0,       -0.90108,           
		0,       0,       0,       0,        0
	}; 

	const float test_data_y[5 * 6] =  {
		0, 0,        0,       0,   0,
		0, 0,        0.90108, 0,   0,
		0, 0.66941,  0,       0,   0,
		0, 0,       -0.66941, 0.90108, 0,
		0, -0.90108, 0,       0,        0, 
		0, 0,        0,      -0.90108, 0, 
		
	}; 
	

	Initialize(size, src_data); 
	Image2NGFFixture<2>::CImageToNGF::OutputPointer result = filter->GetOutput(); 
	filter->Update(); 

	// test size
	BOOST_CHECK_EQUAL(result->GetImageDimension(), 2); 
	
	Image2NGFFixture<2>::CImageToNGF::OutputImageType::RegionType r = result->GetLargestPossibleRegion(); 
	Image2NGFFixture<2>::CImageToNGF::OutputImageType::RegionType::SizeType rsize = r.GetSize(); 
	
	BOOST_CHECK_EQUAL(rsize[0], 5); 
	BOOST_CHECK_EQUAL(rsize[1], 6); 
	
	// test data
	ImageRegionConstIterator<Image2NGFFixture<2>::CImageToNGF::OutputImageType> it(result, r); 
	const float *xtest = test_data_x; 
	const float *ytest = test_data_y; 

	while (! it.IsAtEnd()) {
		BOOST_CHECK_CLOSE( it.Value()[0], *xtest, 0.1f); 
		BOOST_CHECK_CLOSE( it.Value()[1], *ytest, 0.1f); 
		++xtest; 
		++ytest; 
		++it; 
	}
}

// Test evaluation of the metric and the gradient 
BOOST_FIXTURE_TEST_CASE( test_metric_ngf_2d, NGFMetricFixture<2>) 
{
	Image2NGFFixture<2>::CImage::SizeType size = {{5, 6}}; 

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



	const float test_grad_x[5 * 6] =  {
		0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
		0.00000, 0.00000, 0.48976, 0.00000, 0.00000,
		0.00000, 1.34380, 0.00000, 0.48976, 0.00000,
		0.00000, 0.00000, 0.90443, 0.00000, 0.00000,
		0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
		0.00000, 0.00000, 0.00000, 0.00000, 0.00000

	}; 

	const float test_grad_y[5 * 6] =  {
		0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
		0.00000, 0.00000, 1.14901, 0.00000, 0.00000,
		0.00000, 1.34380, 0.00000, 0.65925, 0.00000,
		0.65925, 0.00000, 1.34380, 1.31851, 0.00000,
		0.00000, 0.65925, 0.65925, 0.00000, 0.65925,
		0.00000, 0.00000, 0.00000, 0.65925, 0.00000
	}; 

	Initialize(size,  src_data, ref_data); 

	CTransform::ParametersType params = transform->GetParameters(); 
	
	BOOST_CHECK_CLOSE(metric->GetValue(params), 0.27752, 0.1); 

	CNGFMetric::MovingNGFType::Pointer g = metric->GetGradient(params); 
	Image2NGFFixture<2>::CImageToNGF::OutputImageType::RegionType r = g->GetLargestPossibleRegion(); 

	ImageRegionConstIterator<Image2NGFFixture<2>::CImageToNGF::OutputImageType> it(g, r); 
	const float *xtest = test_grad_x; 
	const float *ytest = test_grad_y; 

	while (! it.IsAtEnd()) {
		BOOST_CHECK_CLOSE( it.Value()[0], *xtest, 0.1f); 
		BOOST_CHECK_CLOSE( it.Value()[1], *ytest, 0.1f); 
		++xtest; 
		++ytest; 
		++it; 
	}
}

// Test evaluation of the derivative of an affine transformation  
BOOST_FIXTURE_TEST_CASE( test_derivative_ngf_2d, NGFMetricFixture<2>) 
{
	Image2NGFFixture<2>::CImage::SizeType size = {{5, 6}}; 

	const float src_data[5 * 6] =  {
		0, 0, 0, 0, 0, 
		0, 0, 0, 0, 0, 
		0, 0, 4, 0, 0, 
		0, 4, 0, 0, 0, 
		0, 0, 0, 4, 0, 
		0, 0, 0, 0, 0
	}; 	

	const float test_der[8] = {
		5.60146, 6.870018, 18.85528, 26.3216, 0, 0, 3.22775, 9.11066
	}; 

	Initialize(size,  src_data, src_data); 
	CTransform::ParametersType params = transform->GetParameters(); 

	CNGFMetric::DerivativeType derivative; 
	CNGFMetric::MeasureType    measure; 

	metric->GetValueAndDerivative(params, measure, derivative); 

	for (size_t i = 0; i < derivative.GetSize(); ++i) 
		BOOST_CHECK_EQUAL(derivative[i], 0.0); 
	
	params[7] = 1; 
	
	// now evaluate the derivative w.r.t. the transform 
		
	metric->GetValueAndDerivative(params, measure, derivative); 

	BOOST_CHECK_CLOSE(measure, 0.27752, 0.1);

	for (size_t i = 0; i < derivative.GetSize(); ++i) 
		BOOST_CHECK_CLOSE(derivative[i]+1.0, test_der[i] + 1.0, 0.1); 


}

/*
  Implementations of the fixture. 
 */

template <int dim>
Image2NGFFixture<dim>::Image2NGFFixture()
{
}


template <int dim>
void Image2NGFFixture<dim>::Initialize(typename CImage::SizeType size, const float *data)
{
	src = create_image(size, data); 
	filter = CImageToNGF::New(); 
	filter->SetInput(src); 
}

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


template <int dim>
void NGFMetricFixture<dim>::Initialize(typename CImage::SizeType size,  
				       const float *src_data, const float *ref_data)
{
	src = create_image(size, src_data); 
	ref = create_image(size, ref_data); 
	interpolator = CInterpolator::New(); 
	transform = CTransform::New(); 

	metric = CNGFMetric::New(); 
	// this kernel was used to evaluate the test data, better would be 
	// to use a mock kernel
	metric->SetEvaluator(new NGFDeltaKernel<typename CNGFMetric::MovingNGFType, 
			     typename CNGFMetric::FixedNGFType>()); 
	metric->Initialize(ref, src, transform, interpolator, ref->GetLargestPossibleRegion()); 

}

