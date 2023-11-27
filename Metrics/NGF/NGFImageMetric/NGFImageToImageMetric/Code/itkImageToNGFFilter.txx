/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkImageToNGFFilter.txx $
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "itkGetImageNoiseFunction.h"

NS_BEGIN(itk)

template <class TInputImage>
ImageToNGFFilter<TInputImage>::ImageToNGFFilter():
m_Noise(-1.0)
{
}

template <class TInputImage>
void ImageToNGFFilter<TInputImage>::SetNoise(double noise)
{
	m_Noise = noise; 
}

template <class TInputImage>
void ImageToNGFFilter<TInputImage>::Update()
{
	Superclass::Update();
	typename OutputImageType::Pointer o = this->GetOutput(); 
	typename TInputImage::RegionType region = o->GetLargestPossibleRegion (); 

	// automatic initialization doesn't seem to be ITK style 
	// OTOH, never let a man do a machines job 
	if (m_Noise <= 0.0)
		SetNoise(GetImageNoise<TInputImage>(this->GetInput()));

	// Evaluate "jump" value
	ImageRegionIterator<GradientType> ig(o, region); 
	double eta = 0.0;
	while (!ig.IsAtEnd()) {
		eta += ig.Value().GetNorm();
		++ig;
	}
	eta *= m_Noise / region.GetNumberOfPixels();

	const double eta2 = eta * eta;
	
	// normalize gradient 
	ig.GoToBegin();
	while (!ig.IsAtEnd()) {
		double n = ig.Value().GetSquaredNorm() + eta2;
		if (n > 0) {
			ig.Set(ig.Value() / sqrt(n));
		}
		++ig;
	}
}

NS_END(itk)
