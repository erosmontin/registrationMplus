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

#ifndef itkImageToNGFFilter_h
#define itkImageToNGFFilter_h

#include <itkGradientImageFilter.h>
#include "itkMyDefines.h"

NS_BEGIN(itk)

/**
   This filter evaluated the Normalized Gradient Field of an image according to 
   E. Haber and J. Modersitzki, "Beyond Mutual Information: A Simple and Robust
   Alternative", in Bildverarbeitung für die Medizin 2005, eds. Meinzer, Handels, 
   Horsch and Tolxdorff, 2005, pp. 350-354
   The input of the filter is a scalar valued image, and the output an vector valued
   image containing the normalized gradient. 
   It is possible to pass a noise value to the filter before executing Update(), 
   however, if no value is given, then the filter will estimate the noise 
   by using the GetImageNoise() function.  
   \todo This code is most likely not in line with the ITK coding standard
   \todo Somehow the GenerateData function should be used instead of Update()
*/

template <class TInputImage>
class ITK_EXPORT ImageToNGFFilter : public GradientImageFilter<TInputImage> {
public: 
	  typedef ImageToNGFFilter                                Self;
	  typedef SmartPointer<Self>                              Pointer;
	  typedef SmartPointer<const Self>                        ConstPointer;
	  itkNewMacro(Self);  
	  itkTypeMacro(ImageToNGFFilter, GradientImageFilter);
	  
	  typedef GradientImageFilter< TInputImage>               Superclass;
	  typedef typename TInputImage::InternalPixelType         PixelType;
	  typedef typename Superclass::OutputImageType            OutputImageType;
	  typedef typename Superclass::OutputImageType::Pointer   OutputPointer;
	  typedef typename Superclass::InputImageType             InputImageType;

	  ImageToNGFFilter(); 

	  virtual void Update(); 
	  void SetNoise(double noise); 
private:  
	  typedef typename Superclass::InputImageRegionType RegionType;
	  typedef typename GradientImageFilter<TInputImage>::OutputImageType GradientType; 
	  double m_Noise; 
}; 

NS_END(itk)

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkImageToNGFFilter.txx"
#endif

#endif
