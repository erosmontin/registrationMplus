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
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "itkNGFMetricKernel.h"

NS_BEGIN(itk)

template <unsigned int dim> 
struct MockConstIterator {
	MockConstIterator(){}; 
	MockConstIterator(const CovariantVector<float,dim>& v); 
	
	const CovariantVector<float,dim>& Value() const; 
private: 
	CovariantVector<float,dim> m_Value; 
}; 

template <typename T, unsigned int dim> 
struct MockGF {

	typedef CovariantVector<float,dim> PixelType; 
	typedef CovariantVector<float,dim> VectorType; 
	//typedef float RealValueType; 
	typedef MockConstIterator<dim> ConstIterator; 
	typedef ConstIterator GradientIteratorType; 
	
	const static unsigned int ImageDimension = dim; 
};

template <typename T, unsigned int dim> 
struct  NGFIteratorTrait<MockGF<T, dim>, dim >{
	typedef typename MockGF<T, dim>::ConstIterator ConstIterator; 
}; 


template <unsigned int dim> 
struct KernelDataFixture {


	typedef MockGF<float, dim> GradientField; 
	typedef typename GradientField::PixelType::RealValueType RealValueType; 
	typedef typename GradientField::PixelType PixelType;  
	typedef typename GradientField::VectorType VectorType; 
	typedef typename GradientField::GradientIteratorType GradientIteratorType; 

	void do_check(const RealValueType test_value, 
		      const VectorType& test_gradient, 
		      const GradientIteratorType igrad[], 
		      const PixelType& moving, 
		      const PixelType& fixed,
		      const NGFKernel<GradientField,GradientField> & kernel); 
	
}; 

struct Kernel2DDataFixture: public KernelDataFixture<2> {
	typedef KernelDataFixture<2>::VectorType VectorType; 
	typedef KernelDataFixture<2>::GradientIteratorType GradientIteratorType; 
	typedef MockGF<float, 2> MockGF2D; 

	Kernel2DDataFixture(); 
	
	void check(const RealValueType& test_value, 
		   const VectorType& test_gradient, 
		   const NGFKernel<GradientField,GradientField> & kernel); 
private: 
	VectorType moving; 
	VectorType fixed; 
	VectorType gx; 
	VectorType gy; 
	GradientIteratorType grad[2]; 
}; 

struct Kernel3DDataFixture: public KernelDataFixture<3> {
	typedef KernelDataFixture<3>::VectorType VectorType; 
	typedef KernelDataFixture<3>::GradientIteratorType GradientIteratorType; 
	typedef MockGF<float, 3> MockGF3D; 

	Kernel3DDataFixture(); 
	
	void check(const RealValueType& test_value, 
		   const VectorType& test_gradient, 
		   const NGFKernel<GradientField,GradientField> & kernel); 
	VectorType moving; 
	VectorType fixed; 
	VectorType gx; 
	VectorType gy; 
	VectorType gz; 
	GradientIteratorType grad[3]; 
}; 

/*
	moving <0.4, 0.8>
	fixed  <0.1, 0.7>
	delta  <0.3, 0.1>
	sum    <0.5, 1.5>

	gx = <1.0, 2.0>
	gy = <3.0, -4.0>
*/

BOOST_FIXTURE_TEST_CASE( test_scalar_ngf_kernel_2d, Kernel2DDataFixture )
{
	Kernel2DDataFixture::VectorType test_grad; 
	Kernel2DDataFixture::RealValueType test_value; 

	test_value = -0.36;  // 0.4 * 0.1 + 0.8 *0.7 = 0.04 + 0.56 = 0.6 -> ^2 0.36
	test_grad[0] = - 1.2 * 1.5; // - 2.0 * 0.6 * ( 0.1 * 1.0 + 2.0 * 0.7) = -1.2 * 1.5 
	test_grad[1] = + 1.2 * 2.5; // - 2.0 * 0.6 * ( 0.1 * 3.0 - 4.0 * 0.7) = -1.2 * (-2.5)
	
	check(test_value, test_grad, NGFScalarKernel<MockGF2D,MockGF2D>()); 
}

BOOST_FIXTURE_TEST_CASE( test_delta2_ngf_kernel_2d, Kernel2DDataFixture )
{
	Kernel2DDataFixture::VectorType test_grad; 
	Kernel2DDataFixture::RealValueType test_value; 

	test_value =    0.09;  // 0.16 + 0.64 - 0.01 - 0.49 = 0.3 (^2) = 0.09
	test_grad[0] =  2.4; // 0.3 * (0.4 * 1.0 + 0.8 * 2.0) = 0.3 (0.4 + 1.6) =   .6 ->(*4)  2.4
	test_grad[1] = -2.4;  // 0.3 * (0.4 * 3.0 - 0.8 * 4.0) = 0.3 (1.2 - 3.2) =  -.6 ->(*4) -2.4
	
	check(test_value, test_grad, NGFDelta2Kernel<MockGF2D,MockGF2D >()); 
}

BOOST_FIXTURE_TEST_CASE( test_delta_ngf_kernel_2d, Kernel2DDataFixture )
{
	Kernel2DDataFixture::VectorType test_grad; 
	Kernel2DDataFixture::RealValueType test_value; 

	test_value = 0.01;  // 0.1 *0.3 + 0.7 *0.1 = 0.1 ->^2 = 0.01
	test_grad[0] =  .3; // 0.1 * (1.0 * 0.1 + 2.0 * 0.7) = 0.1(0.1 + 1.4) =   .15 ->(*2)  0.3
	test_grad[1] = -.5; // 0.1 * (3.0 * 0.1 - 4.0 * 0.7) = 0.1(0.3 - 2.8) =  -.25 ->(*2) -0.5
	
	check(test_value, test_grad, NGFDeltaKernel<MockGF2D,MockGF2D >()); 
}


BOOST_FIXTURE_TEST_CASE( test_scaled_delta_ngf_kernel_2d, Kernel2DDataFixture )
{
	Kernel2DDataFixture::VectorType test_grad; 
	Kernel2DDataFixture::RealValueType test_value; 

	// see "octave deltadot.oct"
	test_value = 0.0023950;  
	test_grad[0] =  .098488; 
	test_grad[1] = -.229804; 
	
	check(test_value, test_grad, NGFScaledDeltaKernel<MockGF2D,MockGF2D>()); 
}


BOOST_FIXTURE_TEST_CASE( test_cross_ngf_kernel_2d, Kernel2DDataFixture )
{
	Kernel2DDataFixture::VectorType test_grad; 
	Kernel2DDataFixture::RealValueType test_value; 

	test_value = 0.04;  // (0.4 * 0.7 - (0.8 * 0.1)  = 0.28 + 0.08 = 0.2 -> ^2 0.04
	test_grad[0] =  .2; // 0.2 * (1.0 * 0.7 - 2.0 * 0.1) = 0.2(0.7 - 0.2) =  .1 (*2)
	test_grad[1] = 1.0; // 0.2 * (3.0 * 0.7 + 4.0 * 0.1) = 0.2(2.1 + 0.4) =  0.5 (*2)
	
	check(test_value, test_grad, NGFCrossKernel<MockGF2D,MockGF2D >()); 
}

/*
  moving = <0.4, 0.8, 0.2>
  fixed  = <0.1, 0.7, 2.0>
  delta  = <0.3, 0.1, -1.8>
  gx     = <1.0, -2.0, 0.5>
  gy     = <3.0, -4.0, -2.0>
  gz     = <1.0, -2.0, -5.0>
*/

BOOST_FIXTURE_TEST_CASE( test_delta_kernel_3d, Kernel3DDataFixture )
{
	Kernel3DDataFixture::VectorType test_grad; 
	Kernel3DDataFixture::RealValueType test_value; 

	test_value = 3.5 * 3.5;  //  - (0.03 + 0.07 - 3.6 = -3.5)^2 = -1
	test_grad[0] = -7 * (0.1 - 1.4 + 1.0);
	test_grad[1] = -7 * (0.3 - 2.8 - 4.0); 
	test_grad[2] = -7 * (0.1 - 1.4 - 10.0); 
	
	check(test_value, test_grad, NGFDeltaKernel<MockGF3D,MockGF3D >()); 
}

BOOST_FIXTURE_TEST_CASE( test_scalar_ngf_kernel_3d, Kernel3DDataFixture )
{
	Kernel3DDataFixture::VectorType test_grad; 
	Kernel3DDataFixture::RealValueType test_value; 

	test_value = -1;  //  - (0.04 + 0.56 + .4)^2 = -1
	test_grad[0] = -2. * (0.1 - 1.4 + 1.0);
	test_grad[1] = -2. * (0.3 - 2.8 - 4.0); 
	test_grad[2] = -2. * (0.1 - 1.4 - 10.0); 
	
	check(test_value, test_grad, NGFScalarKernel<MockGF3D,MockGF3D>()); 
}

BOOST_FIXTURE_TEST_CASE( test_cross_ngf_kernel_3d, Kernel3DDataFixture )
{
	Kernel3DDataFixture::VectorType test_grad; 
	Kernel3DDataFixture::RealValueType test_value; 
	Kernel3DDataFixture::VectorType help; 
	Kernel3DDataFixture::VectorType ghelp0; 
	Kernel3DDataFixture::VectorType ghelp1; 
	Kernel3DDataFixture::VectorType ghelp2; 
	
	help[0] = 2.0 * 0.8 - 0.2 * 0.7; 
	help[1] = 0.2 * 0.1 - 2.0 * 0.4; 
	help[2] = 0.4 * 0.7 - 0.1 * 0.8; 
	
	ghelp0[0] = gx[1] * fixed[2] - gx[2] * fixed[1]; 
	ghelp0[1] = gx[2] * fixed[0] - gx[0] * fixed[2]; 
	ghelp0[2] = gx[0] * fixed[1] - gx[1] * fixed[0]; 


	ghelp1[0] = gy[1] * fixed[2] - gy[2] * fixed[1]; 
	ghelp1[1] = gy[2] * fixed[0] - gy[0] * fixed[2]; 
	ghelp1[2] = gy[0] * fixed[1] - gy[1] * fixed[0]; 

	ghelp2[0] = gz[1] * fixed[2] - gz[2] * fixed[1]; 
	ghelp2[1] = gz[2] * fixed[0] - gz[0] * fixed[2]; 
	ghelp2[2] = gz[0] * fixed[1] - gz[1] * fixed[0]; 

	test_value = help * help; 
	test_grad[0] = 2.0 * (help * ghelp0); 
	test_grad[1] = 2.0 * (help * ghelp1); 
	test_grad[2] = 2.0 * (help * ghelp2); 
	
	check(test_value, test_grad, NGFCrossKernel<MockGF3D,MockGF3D>()); 
}





template <unsigned int dim> 
void KernelDataFixture<dim>::do_check(const RealValueType test_value, 
				      const VectorType& test_gradient, 
				      const GradientIteratorType igrad[], 
				      const PixelType& moving, 
				      const PixelType& fixed,
				      const NGFKernel<MockGF<float, dim>, MockGF<float, dim> >& 
				      kernel)
{
	BOOST_CHECK_CLOSE(kernel.Value(moving, fixed), test_value, 0.1f); 

	VectorType test_result; 
	BOOST_CHECK_CLOSE(kernel.Gradient(test_result, moving, igrad, fixed), test_value, 0.1); 
	for (unsigned int i = 0; i < dim; ++i) 
		BOOST_CHECK_CLOSE(test_result[i], test_gradient[i], 0.1);
}

Kernel2DDataFixture::Kernel2DDataFixture()
{
	moving[0] = 0.4; 
	moving[1] = 0.8; 
	
	fixed[0] = 0.1; 
	fixed[1] = 0.7; 
	
	gx[0] = 1.0; 
	gx[1] = 2.0; 

	gy[0] = 3.0; 
	gy[1] = -4.0; 

	grad[0] = GradientIteratorType(gx); 
	grad[1] = GradientIteratorType(gy); 
}

void Kernel2DDataFixture::check(const RealValueType& test_value, 
				const VectorType& test_gradient, 
				const NGFKernel<GradientField,GradientField> & kernel)
{
	do_check(test_value, test_gradient, grad, moving, fixed, kernel); 
}


Kernel3DDataFixture::Kernel3DDataFixture()
{
	moving[0] = 0.4; 
	moving[1] = 0.8; 
	moving[2] = 0.2; 
	
	fixed[0] = 0.1; 
	fixed[1] = 0.7;
	fixed[2] = 2.0;
 
	
	gx[0] = 1.0; 
	gx[1] =-2.0; 
	gx[2] = 0.5; 

	gy[0] =  3.0; 
	gy[1] = -4.0; 
	gy[2] = -2.0; 

	gz[0] =  1.0; 
	gz[1] = -2.0; 
	gz[2] = -5.0; 

	grad[0] = GradientIteratorType(gx); 
	grad[1] = GradientIteratorType(gy); 
	grad[2] = GradientIteratorType(gz); 
}

void Kernel3DDataFixture::check(const RealValueType& test_value, 
				const VectorType& test_gradient, 
				const NGFKernel<GradientField,GradientField> & kernel)
{
	do_check(test_value, test_gradient, grad, moving, fixed, kernel); 
}

template <unsigned int dim> 
MockConstIterator<dim>::MockConstIterator(const CovariantVector<float,dim>& v):
	m_Value(v)
{
}

template <unsigned int dim> 
const CovariantVector<float,dim>& MockConstIterator<dim>::Value() const
{
	return m_Value; 
}
	
NS_END(itk)
