/*=========================================================================
NOise is the ETA parameter of haber
=========================================================================*/

#include "itkGetImageNoiseFunction.h"
#include <omp.h>
NS_BEGIN(itk)

template <class TInputImage>
ImageToNGFFilter<TInputImage>::ImageToNGFFilter()//:m_Noise(-1.0)
{

}

template <class TInputImage>
void ImageToNGFFilter<TInputImage>::aSetNoise(double noise)
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
		aSetNoise(GetImageNoise<TInputImage>(this->GetInput()));

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
