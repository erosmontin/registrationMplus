/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision:  $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef itkNGFMetricKernel_txx
#define itkNGFMetricKernel_txx

#include <iostream>
#ifdef HAVE_BOOST
#include <boost/static_assert.hpp>
#endif

NS_BEGIN(itk)

//
// Implementation of scalar product based evaluator 
//
template <class MovingGF, class FixedGF> 
inline typename NGFKernel<MovingGF,FixedGF>::Real 
NGFScalarKernel<MovingGF,FixedGF>::Value(const MV& moving, const FV& fixed)const 
{
	Real v = moving * fixed; 
	return - v * v; 
}

template <class MovingGF, class FixedGF> 
inline typename NGFKernel<MovingGF,FixedGF>::Real 
NGFScalarKernel<MovingGF,FixedGF>::Gradient(MV& result, 
				   const MV& moving,
				   const GradientIteratorType gradmoving[], 
				   const FV& fixed)const 
{
	Real v = moving * fixed; 
	#pragma omp simd

	for (unsigned long i = 0; i < result.Size(); ++i) 
		result[i] = - 2.0 * v * (gradmoving[i].Value() * fixed); 
	
	return -v * v; 
}

#if defined(_MSC_VER) && _MSC_VER < 1310 
template <typename M, typename F, unsigned int dim>
struct __DispatchCrossproduct {
	typedef M ReturnType; 
	static void Apply(ReturnType &result,  const M &a,  const F &b)
	{
		/*
		  Well, you insist on using an outdated compiler, so you get a switch 
		  statement in the inner loop. Don't blame me if your program is slow.
		  In addition, in the 2D case you get two more useless FLOPS: 
		  a multiplication 0.0*0.0 and the addition of this zero.
		*/
		switch (dim) {
		case 2:
			result[0] = a[0] * b[1] - a[1] *b[0]; 
			result[1] = 0.0; 
			break; 
		case 3: 
			result[0] = a[1] * b[2] - a[2] *b[1];
			result[1] = a[2] * b[0] - a[0] *b[2]; 
			result[2] = a[0] * b[1] - a[1] *b[0]; 
			break; 
		default:
			itkGenericExceptionMacro("Cross product only supports dimensions 2 and 3 but " << 
					 dim << "was requested"); 
		}

	}
}; 
#else
//
// Implementation of cross product based evaluator 
//
template <typename M, typename F, unsigned int dim>
struct __DispatchCrossproduct {
	typedef float ReturnType; 
	static void Apply(ReturnType &,  const M &,  const F &)
	{
#ifdef HAVE_BOOST
		// if this is instanciated, then the dimension is not supported
		BOOST_STATIC_ASSERT(dim == 2 || dim == 3);
#else
		itkGenericExceptionMacro("Cross product only supports dimensions 2 and 3 but " << 
					 dim << "was requested"); 
#endif
	}
}; 

template <typename M, typename F>
struct __DispatchCrossproduct<M,F,3> {
	typedef M ReturnType; 
	static void Apply(ReturnType& result,
			  const M &a,
			  const F &b)
	{
		result[0] = a[1] * b[2] - a[2] *b[1];
		result[1] = a[2] * b[0] - a[0] *b[2]; 
		result[2] = a[0] * b[1] - a[1] *b[0]; 
	}
}; 

template <typename M, typename F>
struct __DispatchCrossproduct<M,F,2> {
	typedef typename M::RealValueType ReturnType; 
	static void Apply(ReturnType& result,
			  const M &a,
			  const F &b)
	{
		result = static_cast<ReturnType>(a[0] * b[1] - a[1] *b[0]); 
	}
}; 
#endif

template <class MovingGF, class FixedGF> 
inline typename NGFKernel<MovingGF,FixedGF>::Real 
NGFCrossKernel<MovingGF,FixedGF>::Value(const MV& moving, const FV& fixed)const 
{
	typedef __DispatchCrossproduct<MV, FV, MovingGF::ImageDimension> DispatchCrossproduct; 
	typename DispatchCrossproduct::ReturnType x; 
	DispatchCrossproduct::Apply(x, moving, fixed);
	return x*x; 
}

template <class MovingGF, class FixedGF> 
inline typename NGFKernel<MovingGF,FixedGF>::Real 
NGFCrossKernel<MovingGF,FixedGF>::Gradient(MV& result, 
					   const MV& moving,
					   const GradientIteratorType gradmoving[], 
					   const FV& fixed)const 
{
	typedef __DispatchCrossproduct<MV, FV, MovingGF::ImageDimension> DispatchCrossproduct; 
	
	typename DispatchCrossproduct::ReturnType x; 
	typename DispatchCrossproduct::ReturnType help; 

	DispatchCrossproduct::Apply(x, moving, fixed);
	#pragma omp simd
	for (unsigned long i = 0; i < result.Size(); ++i) {
		DispatchCrossproduct::Apply(help, gradmoving[i].Value(), fixed); 
		result[i] = 2.0 * (help * x); 
	}
	return x * x; 
}


//
// Implementation of scaled difference based evaluator 
//


template <class MovingGF, class FixedGF> 
inline typename NGFKernel<MovingGF,FixedGF>::Real 
NGFScaledDeltaKernel<MovingGF,FixedGF>::Value(const MV& moving, const FV& fixed)const 
{
	const Real dot_mf =  moving * fixed;
	const Real dot_mm =  moving * moving;
	const Real dot_ff =  fixed * fixed;
	const Real f = dot_mm * dot_ff;
	const Real cxy = f != 0.0 ? dot_mf / sqrt(f) : 0.0; 
	const Real helper0 = dot_ff - dot_mf * cxy; 

	return 0.5 * (helper0 * helper0); 
}

template <class MovingGF, class FixedGF> 
inline typename NGFKernel<MovingGF,FixedGF>::Real 
NGFScaledDeltaKernel<MovingGF,FixedGF>::Gradient(MV& result, 
					const MV& moving,
					const GradientIteratorType g[], 
					const FV& fixed)const 
{
	const Real dot_mf =  moving * fixed;
	const Real dot_mm =  moving * moving;
	const Real dot_ff =  fixed * fixed;
	const Real f = dot_mm * dot_ff;
	const Real cxy = f != 0.0 ? dot_mf / sqrt(f) : 0.0; 
	const Real dotrsbymm = dot_mm != 0.0f ? dot_mf / dot_mm : 0.0f; 
	const Real helper0 = dot_ff - dot_mf * cxy; 

	const Real helper = - helper0 * cxy; 
	#pragma omp simd
	for (unsigned long i = 0; i < result.Size(); ++i)
		result[i] = helper * ((2.0 * (fixed * g[i].Value())) - dotrsbymm * (moving * g[i].Value())); 

	return 0.5 * helper0 * helper0; 
}

//
// Implementation of the delta cost evalator 
// Considers gradients only similar if they point in the same direction 
// 

template <class MovingGF, class FixedGF> 
inline typename NGFKernel<MovingGF,FixedGF>::Real 
NGFDeltaKernel<MovingGF,FixedGF>::Value(const MV& moving, const FV& fixed)const 
{
	Real v = (moving - fixed) * fixed; 
	return v * v; 
}

template <class MovingGF, class FixedGF> 
inline typename NGFKernel<MovingGF,FixedGF>::Real 
NGFDeltaKernel<MovingGF,FixedGF>::Gradient(MV& result, 
					      const MV& moving,
					      const GradientIteratorType gradmoving[], 
					      const FV& fixed)const 
{
	Real v = (moving - fixed) * fixed; 
	#pragma omp simd
	for (unsigned long i = 0; i < result.Size(); ++i) 
		result[i] = 2.0 * v * (gradmoving[i].Value() * fixed); 
	return v * v; 
}


//
// Implementation of the squared delta cost evalator 
// 

template <class MovingGF, class FixedGF> 
inline typename NGFKernel<MovingGF,FixedGF>::Real 
NGFDelta2Kernel<MovingGF,FixedGF>::Value(const MV& moving, const FV& fixed)const 
{
	Real v = (moving - fixed) * (moving + fixed); 
	return v * v; 
}

template <class MovingGF, class FixedGF> 
inline typename NGFKernel<MovingGF,FixedGF>::Real 
NGFDelta2Kernel<MovingGF,FixedGF>::Gradient(MV& result, 
					      const MV& moving,
					      const GradientIteratorType gradmoving[], 
					      const FV& fixed)const 
{
	Real v = (moving - fixed) * (moving + fixed); 
	#pragma omp simd
	for (unsigned long i = 0; i < result.Size(); ++i) 
		result[i] = 4.0 * v * (gradmoving[i].Value() * moving); 
	return v * v; 
}

NS_END(itk)

#endif
