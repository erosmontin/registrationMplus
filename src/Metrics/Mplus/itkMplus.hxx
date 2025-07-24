#ifndef __itkMplus_hxx
#define __itkMplus_hxx

#include "itkMplus.h"
#include "../NGF/NGFImageMetric/NGFImageToImageMetric/Code/itkNGFMetricKernel.h"
#include "itkMinimumMaximumImageCalculator.h"

// <<< add for auto-eta
#include "itkGradientMagnitudeImageFilter.h"
#include "itkImageRegionIterator.h"
#include <vector>
#include <algorithm>
// >>>

namespace itk
{
	/**
	 * Constructor
	 */
	template <class TFixedImage, class TMovingImage>
	Mplus<TFixedImage, TMovingImage>::Mplus()
	{
		m_MA = MattesType::New();
		m_NGF = NGFType::New();
		m_MSE = MSEType::New();
		m_HMI = HMIType::New();
		m_NMI = NMIType::New();
		m_INTERNALL_interpolator = LFType::New();
		m_Lambda = 0.0;
		m_LambdaDerivative = m_Lambda;
		m_NGFNumberOfSamples = 20000;
		m_MANumberOfSamples = 20000;
		m_MSENumberOfSamples = 20000;
		m_HMINumberOfSamples = 20000;
		m_NMINumberOfSamples = 20000;
		m_BinNumbers = 50;
		m_NumberOfThreads = 1;
		m_Evaluator = 0;
		m_FixedEta = -1;
		m_MovingEta = -1;
		m_Alpha = 1.0;
		m_AlphaDerivative = 1.0;
		m_Nu = 0.0;
		m_NuDerivative = 0.0;
		m_Yota = 0.0;
		m_YotaDerivative = 0.0;
		m_Rho = 0.0;
		m_RhoDerivative = 0.0;
		m_UseCachingOfBSplineWeights = true;
		m_UseExplicitPDFDerivatives = true;
		m_NormalizeDerivatives = false;
		m_NGFSpacing.Fill(4.0);
		m_AutoEstimateEta = false;
		m_RangeDerivatives=0.0;

	}
	template <class TFixedImage, class TMovingImage>
	void
	Mplus<TFixedImage, TMovingImage>::NormalizeDerivative(DerivativeType &derivative) const
	{
		if (m_NormalizeDerivatives)
		{
			this->NormalizeComponents(derivative);
		}
	}

	template <class TFixedImage, class TMovingImage>
	Mplus<TFixedImage, TMovingImage>::~Mplus() {}

	/**
	 * Initialize
	 */

	template <class TFixedImage, class TMovingImage>
	void
	Mplus<TFixedImage, TMovingImage>::Initialize(void)
	{

		Superclass::Initialize();


		if ((this->m_Alpha!=0.0) || (this->m_AlphaDerivative!=0.0))
		{
		m_MA->SetFixedImage(this->m_FixedImage);
		m_MA->SetMovingImage(this->m_MovingImage);
		m_MA->SetInterpolator(this->m_Interpolator);
		m_MA->SetTransform(this->m_Transform);
		m_MA->SetFixedImageRegion(this->m_FixedImage->GetRequestedRegion());
		m_MA->SetNumberOfHistogramBins(this->m_BinNumbers);
		m_MA->SetNumberOfSpatialSamples(this->m_MANumberOfSamples);
		m_MA->SetNumberOfThreads(this->m_NumberOfThreads);
		m_MA->SetUseCachingOfBSplineWeights(this->m_UseCachingOfBSplineWeights);
		m_MA->SetUseExplicitPDFDerivatives(this->m_UseExplicitPDFDerivatives);

		m_MA->ReinitializeSeed();
		m_MA->Initialize();
		}

		if ((this->m_Nu != 0.0) || (this->m_NuDerivative != 0.0))
		{
			m_MSE->SetFixedImage(this->m_FixedImage);
			m_MSE->SetMovingImage(this->m_MovingImage);
			m_MSE->SetInterpolator(this->m_Interpolator);
			m_MSE->SetTransform(this->m_Transform);
			m_MSE->SetFixedImageRegion(this->m_FixedImage->GetRequestedRegion());
			m_MSE->UseAllPixelsOff();
			m_MSE->SetNumberOfThreads(this->m_NumberOfThreads);
			m_MSE->SetUseCachingOfBSplineWeights(this->m_UseCachingOfBSplineWeights);
			m_MSE->SetNumberOfSpatialSamples(this->m_MSENumberOfSamples);
			m_MSE->ReinitializeSeed();
			m_MSE->Initialize();
		}
		if ((this->m_Rho != 0.0) || (this->m_RhoDerivative != 0.0))
		{
			m_HMI->SetFixedImage(this->m_FixedImage);
			m_HMI->SetMovingImage(this->m_MovingImage);
			m_HMI->SetInterpolator(this->m_Interpolator);
			m_HMI->SetTransform(this->m_Transform);
			m_HMI->SetFixedImageRegion(this->m_FixedImage->GetRequestedRegion());

			{
				// fill a 2‑D histogram size (fixed×moving) with BinNumbers in each dimension
				typename HMIType::HistogramSizeType histogramSize;
				histogramSize.Fill(this->m_BinNumbers);
				m_HMI->SetHistogramSize(histogramSize); // :contentReference[oaicite:0]{index=0}
			}
			m_HMI->SetNumberOfSpatialSamples(this->m_HMINumberOfSamples);
			m_HMI->SetNumberOfThreads(this->m_NumberOfThreads);
			m_HMI->UseAllPixelsOff();
			m_HMI->ReinitializeSeed();
			m_HMI->Initialize();
		}
		if (this->m_Yota != 0.0 || this->m_YotaDerivative != 0.0)
		{
			m_NMI->SetTransform(this->GetTransform());
			m_NMI->SetInterpolator(this->m_Interpolator);
			m_NMI->SetFixedImage(this->m_FixedImage);
			m_NMI->SetMovingImage(this->m_MovingImage);
			m_NMI->SetFixedImageRegion(this->m_FixedImage->GetRequestedRegion());

			// histogram bins (v3 API)
			{
				typename NMIType::HistogramSizeType histSize;
				histSize.Fill(this->m_BinNumbers);
				m_NMI->SetHistogramSize(histSize);
			}

			m_NMI->SetNumberOfSpatialSamples(this->m_NMINumberOfSamples);
			m_NMI->SetNumberOfThreads(this->m_NumberOfThreads);
			m_NMI->UseAllPixelsOff();
			m_NMI->ReinitializeSeed();
			m_NMI->Initialize();

		}
		// add a resampling filter to the NGF metric
		if ((m_Lambda != 0) || (m_LambdaDerivative != 0))
		{

			if (this->m_AutoEstimateEta)
			{
				std::cout << "Auto-estimating η for NGF metric..." << std::endl;
				using GradFilterType = itk::GradientMagnitudeImageFilter<FixedImageType, FixedImageType>;
				typename GradFilterType::Pointer gradFilter = GradFilterType::New();
				gradFilter->SetInput(this->m_FixedImage);
				gradFilter->Update();
	
				std::vector<double> mags;
				mags.reserve(100000);
				itk::ImageRegionConstIterator<FixedImageType> it(
					gradFilter->GetOutput(),
					gradFilter->GetOutput()->GetLargestPossibleRegion());
				size_t count = 0;
				for (; !it.IsAtEnd(); ++it)
				{
					mags.push_back(it.Get());
					if (++count >= 100000)
						break;
				}
				std::sort(mags.begin(), mags.end());
				const double percentile = 0.10;
				size_t idx = static_cast<size_t>(percentile * mags.size());
				double eta = mags[idx];
				constexpr double kMinEta = 1e-8;
				if (eta < kMinEta)
					eta = kMinEta;
				this->SetFixedEta(eta);
				this->SetMovingEta(eta);
				std::cout << "Estimated η: " << eta << std::endl;
			}

			// Desired NGF spacing per dimension
			const typename TFixedImage::SpacingType newSpacing = this->m_NGFSpacing;

			// Get original size and spacing
			typename TFixedImage::SizeType originalSize = this->m_FixedImage->GetLargestPossibleRegion().GetSize();
			typename TFixedImage::SpacingType originalSpacing = this->m_FixedImage->GetSpacing();

			// Compute downsampled size so that
			// downsampledSize[d] * newSpacing[d] >= originalSize[d] * originalSpacing[d]
			typename TFixedImage::SizeType downsampledSize;
			for (unsigned int d = 0; d < FixedImageType::ImageDimension; ++d)
			{
				const double extent = static_cast<double>(originalSize[d]) * originalSpacing[d];
				downsampledSize[d] = static_cast<typename SizeType::SizeValueType>(
					std::ceil(extent / newSpacing[d]));
			}

			// Fixed image resampler
			using FixedResampleFilterType = itk::ResampleImageFilter<FixedImageType, FixedImageType>;
			auto fixedResampler = FixedResampleFilterType::New();
			fixedResampler->SetInput(this->m_FixedImage);
			fixedResampler->SetSize(downsampledSize);
			fixedResampler->SetOutputSpacing(newSpacing);
			fixedResampler->SetOutputOrigin(this->m_FixedImage->GetOrigin());
			fixedResampler->SetOutputDirection(this->m_FixedImage->GetDirection());
			fixedResampler->Update();
			auto downsampledFixed = fixedResampler->GetOutput();

			// Moving image resampler
			using MovingResampleFilterType = itk::ResampleImageFilter<MovingImageType, MovingImageType>;
			auto movingResampler = MovingResampleFilterType::New();
			movingResampler->SetInput(this->m_MovingImage);
			movingResampler->SetSize(downsampledSize);
			movingResampler->SetOutputSpacing(newSpacing);
			movingResampler->SetOutputOrigin(this->m_MovingImage->GetOrigin());
			movingResampler->SetOutputDirection(this->m_MovingImage->GetDirection());
			movingResampler->Update();
			auto downsampledMoving = movingResampler->GetOutput();

			// Set downsampled images on NGF
			m_NGF->SetFixedImage(downsampledFixed);
			m_NGF->SetMovingImage(downsampledMoving);

			m_NGF->SetInterpolator(this->m_INTERNALL_interpolator);
			m_NGF->SetTransform(this->m_Transform);
			m_NGF->SetNumberOfSpatialSamples(this->m_NGFNumberOfSamples);
			m_NGF->SetNumberOfThreads(this->m_NumberOfThreads);
			m_NGF->SetUseCachingOfBSplineWeights(this->m_UseCachingOfBSplineWeights);
			m_NGF->UseAllPixelsOff();
			// RegionType inputRegion = this->m_FixedImage->GetRequestedRegion();
			RegionType inputRegion = this->m_NGF->GetFixedImage()->GetRequestedRegion();
			m_NGF->SetFixedNoise(this->m_FixedEta); // this NGF noise can be considered the eta parameter
			m_NGF->SetMovingNoise(this->m_MovingEta);
			m_NGF->SetFixedImageRegion(inputRegion);

			switch (m_Evaluator)
			{
			case (0):
				//
				// Implementation of scalar product based evaluator
				//
				m_NGF->SetEvaluator(new NGFScalarKernel<MovingNGFType, FixedNGFType>());
				break;
			//
			// Implementation of cross product based evaluator
			//
			case (1):
				m_NGF->SetEvaluator(new NGFCrossKernel<MovingNGFType, FixedNGFType>());
				break;
			//
			// Implementation of scaled difference based evaluator
			//
			case (2):
				m_NGF->SetEvaluator(new NGFScaledDeltaKernel<MovingNGFType, FixedNGFType>());
				break;
			//
			// Implementation of the delta cost evalator
			// Considers gradients only similar if they point in the same direction
			//
			case (3):
				m_NGF->SetEvaluator(new NGFDeltaKernel<MovingNGFType, FixedNGFType>());
				break;
			//
			// Implementation of the squared delta cost evalator
			//
			case (4):
				m_NGF->SetEvaluator(new NGFDelta2Kernel<MovingNGFType, FixedNGFType>());
				break;
			default:
				m_NGF->SetEvaluator(new NGFScalarKernel<MovingNGFType, FixedNGFType>());
			}

			m_NGF->ReinitializeSeed();
			m_NGF->Initialize();
		}
	}
	template<class TFixedImage, class TMovingImage>
	double
	Mplus<TFixedImage,TMovingImage>
	::ComputeDerivativeMean(const DerivativeType & der) const
	{
	  const auto N = der.Size();
	  if (N == 0) return 0.0;
	  double sum = 0.0;

	  for (unsigned i = 0; i < N; ++i)
	  {
		sum += der[i];
	  }
	  return sum / static_cast<double>(N);
	}
	
	template<class TFixedImage, class TMovingImage>
	double
	Mplus<TFixedImage,TMovingImage>
	::ComputeDerivativeStdDev(const DerivativeType & der) const
	{
	  const auto N = der.Size();
	  if (N == 0) return 0.0;
	  const double mean = this->ComputeDerivativeMean(der);
	  double sumSq = 0.0;
	  for (unsigned i = 0; i < N; ++i)
	  {
		const double diff = der[i] - mean;
		sumSq += diff * diff;
	  }
	  // population standard deviation:
	  return std::sqrt( sumSq / static_cast<double>(N) );
	}

	template <class TFixedImage, class TMovingImage>
	typename Mplus<TFixedImage, TMovingImage>::MeasureType
	Mplus<TFixedImage, TMovingImage>::GetValue(const ParametersType &parameters) const
	{
		double a, b, c, d, e;

		a = 0.0;
		if (this->m_Alpha != 0.0)
			a = this->GetMAValue(parameters);
		b = 0.0;
		if (this->m_Lambda != 0.0)
			b = this->GetNGFValue(parameters);
		c = 0.0;
		if (this->m_Nu != 0.0)
			c = this->GetMSEValue(parameters);
		d = 0.0;

		if (this->m_Rho != 0.0)
			d = this->GetHMIValue(parameters);

		e = 0.0;
		if (this->m_Yota != 0.0)
			e = this->GetNMIValue(parameters);

		// std::cout << "a: " << a << " b: " << b << " c: " << c << " d: " << d << std::endl;
		return a + b + c + d + e;
	}

	template <class TFixedImage, class TMovingImage>
	typename Mplus<TFixedImage, TMovingImage>::MeasureType
	Mplus<TFixedImage, TMovingImage>::GetNGFValue(const ParametersType &parameters) const
	{
		return static_cast<MeasureType>(m_NGF->GetValue(parameters) * this->m_Lambda);
	}

	template <class TFixedImage, class TMovingImage>
	typename Mplus<TFixedImage, TMovingImage>::MeasureType
	Mplus<TFixedImage, TMovingImage>::GetMAValue(const ParametersType &parameters) const
	{
		return static_cast<MeasureType>(m_MA->GetValue(parameters) * this->m_Alpha);
	}

	template <class TFixedImage, class TMovingImage>
	typename Mplus<TFixedImage, TMovingImage>::MeasureType
	Mplus<TFixedImage, TMovingImage>::GetMSEValue(const ParametersType &parameters) const
	{
		return static_cast<MeasureType>(m_MSE->GetValue(parameters) * this->m_Nu);
	}

	template <class TFixedImage, class TMovingImage>
	typename Mplus<TFixedImage, TMovingImage>::MeasureType
	Mplus<TFixedImage, TMovingImage>::GetHMIValue(const ParametersType &parameters) const
	{
		return static_cast<MeasureType>(m_HMI->GetValue(parameters) * this->m_Yota);
	}
	template <class TFixedImage, class TMovingImage>
	typename Mplus<TFixedImage, TMovingImage>::MeasureType
	Mplus<TFixedImage, TMovingImage>::GetNMIValue(const ParametersType &parameters) const
	{

		return static_cast<MeasureType>(m_NMI->GetValue(parameters) * this->m_Yota);
	}

	template <class TFixedImage, class TMovingImage>
	void
	Mplus<TFixedImage, TMovingImage>::GetDerivative(const ParametersType &parameters, DerivativeType &derivative) const
	{

		
		


		DerivativeType a;
		a = parameters;		
		if (this->m_AlphaDerivative != 0.0)
			this->GetMADerivative(parameters, a);
		else
			a.Fill(0.0);

		DerivativeType b;
		b = parameters;
		if (this->m_LambdaDerivative != 0.0)
			this->GetNGFDerivative(parameters, b);
		else
			b.Fill(0.0);

		DerivativeType c;
		c = parameters;
		if (this->m_NuDerivative != 0.0)
			this->GetMSEDerivative(parameters, c);
		else
			c.Fill(0.0);

		DerivativeType d;
		d = parameters;
		if (this->m_RhoDerivative != 0.0)
			this->GetHMIDerivative(parameters, d);
		else
			d.Fill(0.0);
		DerivativeType e;
		e = parameters;
		if (this->m_YotaDerivative != 0.0)
			this->GetNMIDerivative(parameters, e);
		else
			e.Fill(0.0);

		derivative = a;
#pragma omp parallel for
		for (long unsigned int p = 0; p < derivative.GetSize(); ++p)
		{
			derivative[p] =
				a[p] * this->m_AlphaDerivative + this->m_LambdaDerivative * b[p] + this->m_NuDerivative * c[p] + this->m_RhoDerivative * d[p] + this->m_YotaDerivative * e[p];
			if (this->m_NormalizeDerivatives)
			{
				// normalize the derivative

				derivative[p] *= this->m_LastComponentNorm;
			}
		}
	}

	template<class TFixedImage, class TMovingImage>
	double
	Mplus<TFixedImage,TMovingImage>
	::ComputeDerivativeRange(const DerivativeType & der) const
	{
	  const auto N = der.Size();
	  if(N == 0) return 0.0;
	  double minVal = std::numeric_limits<double>::infinity();
	  double maxVal = -std::numeric_limits<double>::infinity();

	  for(unsigned i = 0; i < N; ++i)
	  {
		const double v = der[i];
		if(v < minVal) minVal = v;
		if(v > maxVal) maxVal = v;
	  }
	  return maxVal - minVal;
	}



	template<class TFixedImage, class TMovingImage>
	double
	Mplus<TFixedImage,TMovingImage>
	::ComputeDerivativeNorm(const DerivativeType & der) const
	{
	  
	double norm = 0.0;
	#pragma omp parallel for reduction(+:norm)
	for (unsigned int i = 0; i < der.size(); ++i) {
		norm += der[i] * der[i];
	}
	return std::sqrt(norm);
}

	template <class TFixedImage, class TMovingImage>
	void
	Mplus<TFixedImage, TMovingImage>::GetMADerivative(const ParametersType &parameters, DerivativeType &derivative) const
	{
		m_MA->GetDerivative(parameters, derivative);
		// normalize the derivative
		// this->m_MeanDerivatives = this->ComputeDerivativeMean(derivative);
		// this->m_STDDerivatives = this->ComputeDerivativeStdDev(derivative);
		this->m_LastComponentNorm = this->ComputeDerivativeNorm(derivative);
		this->NormalizeDerivative(derivative);
		
	}

	template <class TFixedImage, class TMovingImage>
	void
	Mplus<TFixedImage, TMovingImage>::GetHMIDerivative(const ParametersType &parameters, DerivativeType &derivative) const
	{
		m_HMI->GetDerivative(parameters, derivative);
		this->NormalizeDerivative(derivative);
	}

	template <class TFixedImage, class TMovingImage>
	void
	Mplus<TFixedImage, TMovingImage>::GetNMIDerivative(const ParametersType &parameters, DerivativeType &derivative) const
	{
		m_NMI->GetDerivative(parameters, derivative);

		this->NormalizeDerivative(derivative);
	}
	template <class TFixedImage, class TMovingImage>
	void
	Mplus<TFixedImage, TMovingImage>::GetNGFDerivative(const ParametersType &parameters, DerivativeType &derivative) const
	{

		m_NGF->GetDerivative(parameters, derivative);

		this->NormalizeDerivative(derivative);
	}

	template <class TFixedImage, class TMovingImage>
	void
	Mplus<TFixedImage, TMovingImage>::GetMSEDerivative(const ParametersType &parameters, DerivativeType &derivative) const
	{
		m_MSE->GetDerivative(parameters, derivative);

		this->NormalizeDerivative(derivative);
	}

	template <class TFixedImage, class TMovingImage>
	void
	Mplus<TFixedImage, TMovingImage>::GetValueAndDerivative(const ParametersType &parameters, MeasureType &Value, DerivativeType &Derivative) const
	{
		Value = this->GetValue(parameters);
		this->GetDerivative(parameters, Derivative);
	}

	template <class TImageType, class TMovingImage>
	void
	Mplus<TImageType, TMovingImage>::
		PrintSelf(std::ostream &os, Indent indent) const
	{
		Superclass::PrintSelf(os, indent);
	}

} // end namespace itk

#endif
