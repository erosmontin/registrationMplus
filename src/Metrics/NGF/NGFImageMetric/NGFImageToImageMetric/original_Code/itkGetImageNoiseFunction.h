/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkImageToNGFFilter.h $
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef itkGetImageNoise_h
#define itkGetImageNoise_h

#include <stdexcept>
#include <vector>
#include <itkNoiseImageFilter.h>
#include "itkMyDefines.h"

NS_BEGIN(itk)

/**
   This function takes an input image and evaluates its noise as the 
   median of the regional noise evaluated by using NoiseImageFilter
   with its default parameters. 
   \param image input image 
   \returns the estimated noise
*/

template <class TInputImage>
typename TInputImage::PixelType ITK_EXPORT GetImageNoise(const TInputImage *image) 
{
	typedef typename TInputImage::PixelType PixelType; 
	typedef NoiseImageFilter<TInputImage ,TInputImage> NoiseImageFilterType; 
	typedef typename NoiseImageFilterType::Pointer     NoiseImageFilterPointer; 

	typename TInputImage::RegionType region = image->GetLargestPossibleRegion (); 
	
	NoiseImageFilterPointer noise_filter = NoiseImageFilterType::New(); 
	noise_filter->SetInput(image); 
	typename TInputImage::Pointer noise_image = noise_filter->GetOutput(); 
	noise_image->Update(); 
	
	ImageRegionConstIterator<TInputImage> it(noise_image, region); 
	std::vector<PixelType> help(region.GetNumberOfPixels()); 
	typename std::vector<PixelType>::iterator hi = help.begin(); 
	
	it.GoToBegin(); 
	while (!it.IsAtEnd()) {
		*hi  = it.Value(); 
		++it; 
		++hi; 
	}
	
	if (help.empty()) 
		throw std::invalid_argument("GetImageNoise: Input image is empty"); 
	
	typename std::vector<PixelType>::iterator noise = help.begin() + help.size() / 2; 
	std::nth_element(help.begin(), noise , help.end()); 
	return  *noise; 
}

NS_END(itk)
#endif
