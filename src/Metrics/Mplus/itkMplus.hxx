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

// <<< label metric support
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include <set>
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
		m_GD = GDType::New();
		m_NMI = NMIType::New();
		m_INTERNALL_interpolator = LFType::New();
		m_Lambda = 0.0;
		m_LambdaDerivative = m_Lambda;
		m_NGFNumberOfSamples = 20000;
		m_MANumberOfSamples = 20000;
		m_MSENumberOfSamples = 20000;
		m_NCNumberOfSamples = 20000;
		m_GDNumberOfSamples = 20000;
		m_NMINumberOfSamples = 20000;
		m_NMIBinNumbers = 64;
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
		m_Sigma = 0.0;
		m_SigmaDerivative = 0.0;
		m_UseCachingOfBSplineWeights = true;
		m_UseExplicitPDFDerivatives = true;
		m_DerivativeMode = 0;
		m_MainMetricIndex = 0;
		m_NGFSpacing.Fill(4.0);
		m_AutoEstimateEta = false;
		m_RangeDerivatives=0.0;
		m_NGFPrecomputeGradient = false;
		m_ComputeOverlap   = true;    // default: compute overlap
		m_OverlapPadding   = 20;
		m_FixedImageThreshold   = 0.0;
		m_UseFixedImageThreshold = false;

		// cached scale factors (mode 2)
		m_ScaleMA = 1.0; m_ScaleNGF = 1.0; m_ScaleMSE = 1.0;
		m_ScaleNC = 1.0; m_ScaleLabel = 1.0;
		m_ScaleGD = 1.0; m_ScaleNMI = 1.0;

		// label metric defaults
		m_LabelKappa           = 0.0;
		m_LabelKappaDerivative = 0.0;
		m_LabelNumberOfSamples = 20000;

		// per-sub-metric cached values
		m_LastValMI    = 0.0;
		m_LastValNGF   = 0.0;
		m_LastValMSE   = 0.0;
		m_LastValNC    = 0.0;
		m_LastValGD    = 0.0;
		m_LastValNMI   = 0.0;
		m_LastValLabel = 0.0;
		m_LastValTotal = 0.0;
	}
	template <class TFixedImage, class TMovingImage>
	void
	Mplus<TFixedImage, TMovingImage>::NormalizeDerivative(DerivativeType &derivative) const
	{
		if (m_DerivativeMode == 1)
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



			// 1) grab the fixed‐image region
	auto fixedRegion = this->m_FixedImage->GetRequestedRegion();
 // 1) grab the fixed-image region
	typename FixedImageType::RegionType overlap;

 if ( this->m_ComputeOverlap )
 {
	// 2) compute all 8 corners of the moving image in physical space
	auto movingLargest  = this->m_MovingImage->GetLargestPossibleRegion();
	std::vector<typename FixedImageType::PointType> physPts;
	physPts.reserve(8);
	for (unsigned corner=0; corner<8; ++corner)
	{
	typename MovingImageType::IndexType idx;
	// build the corner index: each bit of ‘corner’ chooses min/max
	for (unsigned d=0; d<FixedImageType::ImageDimension; ++d)
		idx[d] = ((corner>>d)&1)
				? movingLargest.GetIndex()[d] + static_cast<long>(movingLargest.GetSize()[d]) - 1
				: movingLargest.GetIndex()[d];
	// to physical
	typename MovingImageType::PointType p;
	this->m_MovingImage->TransformIndexToPhysicalPoint(idx, p);
	// through current transform to fixed‐space
	physPts.push_back(this->m_Transform->TransformPoint(p));
	}

	// 3) convert those physical points into fixed‐image indices & find min/max
	typename FixedImageType::IndexType minIdx, maxIdx;
	for (unsigned d=0; d<FixedImageType::ImageDimension; ++d)
	{
	minIdx[d] = std::numeric_limits<long>::max();
	maxIdx[d] = std::numeric_limits<long>::min();
	}
	for (auto &p : physPts)
	{
	typename FixedImageType::IndexType idx;
	this->m_FixedImage->TransformPhysicalPointToIndex(p, idx);
	for (unsigned d=0; d<FixedImageType::ImageDimension; ++d)
	{
		minIdx[d] = std::min(minIdx[d], idx[d]);
		maxIdx[d] = std::max(maxIdx[d], idx[d]);
	}
	}

	// 4) build the overlapping region & crop it against fixedRegion
	typename FixedImageType::SizeType   ovSize;
	for (unsigned d=0; d<FixedImageType::ImageDimension; ++d)
	{
	ovSize[d] = maxIdx[d] - minIdx[d] + 1;
	overlap.SetIndex(d, minIdx[d]);
	overlap.SetSize(d,  ovSize[d]);
	}
	overlap.Crop(fixedRegion);



		const unsigned int pad = this->m_OverlapPadding;

		// compute padded min/max in fixed‐index space
		typename FixedImageType::IndexType paddedMin, paddedMax;
		auto fixedIdx0 = fixedRegion.GetIndex();
		auto fixedSz  = fixedRegion.GetSize();
		for (unsigned d = 0; d < FixedImageType::ImageDimension; ++d)
		{
		// clamp so we don’t go outside the fixed image
		long low  = std::max(fixedIdx0[d],  minIdx[d] - static_cast<long>(pad));
		long high = std::min(fixedIdx0[d] + static_cast<long>(fixedSz[d]) - 1,
							maxIdx[d] + static_cast<long>(pad));
		paddedMin[d] = low;
		paddedMax[d] = high;
		}

		// build the new region
		typename FixedImageType::SizeType paddedSize;
		for (unsigned d = 0; d < FixedImageType::ImageDimension; ++d)
		paddedSize[d] = paddedMax[d] - paddedMin[d] + 1;

		overlap.SetIndex(paddedMin);
		overlap.SetSize(paddedSize);


	 }
	 else
	 {
		 // if we do not compute the overlap, use the fixed image region
		 overlap =fixedRegion;
	 }

		if ((this->m_Alpha!=0.0) || (this->m_AlphaDerivative!=0.0))
		{
		m_MA->SetFixedImage(this->m_FixedImage);
		m_MA->SetMovingImage(this->m_MovingImage);
		m_MA->SetInterpolator(this->m_Interpolator);
		m_MA->SetTransform(this->m_Transform);
		// m_MA->SetFixedImageRegion(this->m_FixedImage->GetRequestedRegion());
		m_MA->SetFixedImageRegion(overlap);
		m_MA->UseAllPixelsOff(); // use all pixels is not implemented yet
		m_MA->SetNumberOfHistogramBins(this->m_BinNumbers);
		m_MA->SetNumberOfSpatialSamples(this->m_MANumberOfSamples);
		m_MA->SetNumberOfThreads(this->m_NumberOfThreads);
		m_MA->SetUseCachingOfBSplineWeights(this->m_UseCachingOfBSplineWeights);
		m_MA->SetUseExplicitPDFDerivatives(this->m_UseExplicitPDFDerivatives);

		if (this->m_UseFixedImageThreshold)
			m_MA->SetFixedImageSamplesIntensityThreshold(this->m_FixedImageThreshold);
		m_MA->ReinitializeSeed();
		m_MA->Initialize();
		}

		if ((this->m_Nu != 0.0) || (this->m_NuDerivative != 0.0))
		{
			m_MSE->SetFixedImage(this->m_FixedImage);
			m_MSE->SetMovingImage(this->m_MovingImage);
			m_MSE->SetInterpolator(this->m_Interpolator);
			m_MSE->SetTransform(this->m_Transform);
			// m_MSE->SetFixedImageRegion(this->m_FixedImage->GetRequestedRegion());
			m_MSE->SetFixedImageRegion(overlap);
			m_MSE->UseAllPixelsOff();
			m_MSE->SetNumberOfThreads(this->m_NumberOfThreads);
			m_MSE->SetUseCachingOfBSplineWeights(this->m_UseCachingOfBSplineWeights);
			m_MSE->SetNumberOfSpatialSamples(this->m_MSENumberOfSamples);
			if (this->m_UseFixedImageThreshold)
				m_MSE->SetFixedImageSamplesIntensityThreshold(this->m_FixedImageThreshold);
			m_MSE->ReinitializeSeed();
			m_MSE->Initialize();
		}

		if (this->m_Yota != 0.0 || this->m_YotaDerivative != 0.0)
		{
			m_NC = NCType::New();
			m_NC->SetFixedImage(this->GetFixedImage());
			m_NC->SetMovingImage(this->GetMovingImage());
			m_NC->SetTransform(this->GetTransform());
			m_NC->SetInterpolator(this->GetInterpolator());
			// m_NC->SetFixedImageRegion(this->GetFixedImage()->GetRequestedRegion());
			m_NC->SetFixedImageRegion(overlap);
			m_NC->UseAllPixelsOff(); // use all pixels is not implemented yet
			m_NC->SetNumberOfSpatialSamples(this->m_NCNumberOfSamples);
			m_NC->SetNumberOfThreads(this->m_NumberOfThreads);
			m_NC->SetUseCachingOfBSplineWeights(this->m_UseCachingOfBSplineWeights);
			if (this->m_UseFixedImageThreshold)
				m_NC->SetFixedImageSamplesIntensityThreshold(this->m_FixedImageThreshold);
			m_NC->Initialize();

		}

		// ── Gradient Difference (GD) ──────────────────────────────────────────
		if (this->m_Rho != 0.0 || this->m_RhoDerivative != 0.0)
		{
			m_GD = GDType::New();
			m_GD->SetFixedImage(this->GetFixedImage());
			m_GD->SetMovingImage(this->GetMovingImage());
			m_GD->SetTransform(this->GetTransform());
			m_GD->SetInterpolator(this->GetInterpolator());
			m_GD->SetFixedImageRegion(overlap);
			m_GD->SetDerivativeDelta(0.001);
			m_GD->Initialize();
		}

		// ── Normalized Mutual Information (NMI) ──────────────────────────────
		if (this->m_Sigma != 0.0 || this->m_SigmaDerivative != 0.0)
		{
			m_NMI = NMIType::New();
			m_NMI->SetFixedImage(this->GetFixedImage());
			m_NMI->SetMovingImage(this->GetMovingImage());
			m_NMI->SetTransform(this->GetTransform());
			m_NMI->SetInterpolator(this->GetInterpolator());
			m_NMI->SetFixedImageRegion(overlap);
			// Histogram size: [bins_fixed, bins_moving]
			typename NMIType::HistogramType::SizeType histSize(2);
			histSize.Fill(static_cast<typename NMIType::HistogramType::SizeType::ValueType>(this->m_NMIBinNumbers));
			m_NMI->SetHistogramSize(histSize);
			m_NMI->Initialize();
		}
		// add a resampling filter to the NGF metric
		if ((m_Lambda != 0) || (m_LambdaDerivative != 0))
		{

			if (this->m_AutoEstimateEta)
			{
				std::cout << "Auto-estimating η for NGF metric..." << std::endl;
				constexpr double kMinEta = 1e-8;
				constexpr double percentile = 0.10;
				constexpr double sampleFraction = 0.05; // sample ~5% of voxels

				// Fixed image eta
				{
					using GradFilterType = itk::GradientMagnitudeImageFilter<FixedImageType, FixedImageType>;
					typename GradFilterType::Pointer gradFilter = GradFilterType::New();
					gradFilter->SetInput(this->m_FixedImage);
					gradFilter->Update();

					const unsigned long totalPixF = this->m_FixedImage->GetLargestPossibleRegion().GetNumberOfPixels();
					const size_t maxSamplesF = std::max(static_cast<size_t>(10000),
					                                    static_cast<size_t>(totalPixF * sampleFraction));
					std::vector<double> mags;
					mags.reserve(maxSamplesF);
					itk::ImageRegionConstIterator<FixedImageType> it(
						gradFilter->GetOutput(),
						gradFilter->GetOutput()->GetLargestPossibleRegion());
					for (; !it.IsAtEnd() && mags.size() < maxSamplesF; ++it)
						mags.push_back(it.Get());
					std::sort(mags.begin(), mags.end());
					double etaF = mags[static_cast<size_t>(percentile * mags.size())];
					if (etaF < kMinEta) etaF = kMinEta;
					this->SetFixedEta(etaF);
					std::cout << "  Fixed η: " << etaF << std::endl;
				}

				// Moving image eta (computed independently)
				{
					using GradFilterType = itk::GradientMagnitudeImageFilter<MovingImageType, MovingImageType>;
					typename GradFilterType::Pointer gradFilter = GradFilterType::New();
					gradFilter->SetInput(this->m_MovingImage);
					gradFilter->Update();

					const unsigned long totalPixM = this->m_MovingImage->GetLargestPossibleRegion().GetNumberOfPixels();
					const size_t maxSamplesM = std::max(static_cast<size_t>(10000),
					                                    static_cast<size_t>(totalPixM * sampleFraction));
					std::vector<double> mags;
					mags.reserve(maxSamplesM);
					itk::ImageRegionConstIterator<MovingImageType> it(
						gradFilter->GetOutput(),
						gradFilter->GetOutput()->GetLargestPossibleRegion());
					for (; !it.IsAtEnd() && mags.size() < maxSamplesM; ++it)
						mags.push_back(it.Get());
					std::sort(mags.begin(), mags.end());
					double etaM = mags[static_cast<size_t>(percentile * mags.size())];
					if (etaM < kMinEta) etaM = kMinEta;
					this->SetMovingEta(etaM);
					std::cout << "  Moving η: " << etaM << std::endl;
				}
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

			// Rescale downsampled images to [0,1] for consistent η across scanners
			using RescaleFixedType = itk::RescaleIntensityImageFilter<FixedImageType, FixedImageType>;
			typename RescaleFixedType::Pointer rfx = RescaleFixedType::New();
			rfx->SetInput(downsampledFixed);
			rfx->SetOutputMinimum(0.0);
			rfx->SetOutputMaximum(1.0);
			rfx->Update();
			typename FixedImageType::Pointer scaledFixed = rfx->GetOutput();
			scaledFixed->DisconnectPipeline();

			using RescaleMovingType = itk::RescaleIntensityImageFilter<MovingImageType, MovingImageType>;
			typename RescaleMovingType::Pointer rmv = RescaleMovingType::New();
			rmv->SetInput(downsampledMoving);
			rmv->SetOutputMinimum(0.0);
			rmv->SetOutputMaximum(1.0);
			rmv->Update();
			typename MovingImageType::Pointer scaledMoving = rmv->GetOutput();
			scaledMoving->DisconnectPipeline();

			// Set downsampled + rescaled images on NGF
			m_NGF->SetFixedImage(scaledFixed);
			m_NGF->SetMovingImage(scaledMoving);

			m_NGF->SetInterpolator(this->m_INTERNALL_interpolator);
			m_NGF->SetTransform(this->m_Transform);
			m_NGF->SetNumberOfSpatialSamples(this->m_NGFNumberOfSamples);
			m_NGF->SetNumberOfThreads(this->m_NumberOfThreads);
			m_NGF->SetUseCachingOfBSplineWeights(this->m_UseCachingOfBSplineWeights);
			m_NGF->UseAllPixelsOff();
			RegionType ngfFullRegion = m_NGF->GetFixedImage()->GetLargestPossibleRegion();
			m_NGF->SetFixedNoise(this->m_FixedEta);
			m_NGF->SetMovingNoise(this->m_MovingEta);
			// Map the overlap region from original image space to downsampled NGF space
			{
				RegionType ngfOverlap;
				for (unsigned d = 0; d < FixedImageType::ImageDimension; ++d)
				{
					const double scale = originalSpacing[d] / newSpacing[d];
					long newIdx  = static_cast<long>(std::floor(overlap.GetIndex()[d] * scale));
					long newSz   = static_cast<long>(std::ceil(overlap.GetSize()[d] * scale));
					if (newIdx < 0) newIdx = 0;
					ngfOverlap.SetIndex(d, newIdx);
					ngfOverlap.SetSize(d, static_cast<unsigned long>(newSz));
				}
				ngfOverlap.Crop(ngfFullRegion);
				m_NGF->SetFixedImageRegion(ngfOverlap);
			}

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
			m_NGF->SetPrecomputeGradient(this->m_NGFPrecomputeGradient);
			m_NGF->Initialize();
		}

		// Initialize label metric (distance maps) if label maps have been set
		this->InitializeLabelMetric();
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
			d = this->GetGDValue(parameters);

		e = 0.0;
		if (this->m_Yota != 0.0)
			e = this->GetNCValue(parameters);

		double f = 0.0;
		if (this->m_LabelKappa != 0.0 && m_FixedLabelMap && m_MovingLabelMap)
			f = this->GetKappaValue(parameters);

		double g = 0.0;
		if (this->m_Sigma != 0.0)
			g = this->GetNMIValue(parameters);

		// Cache per-sub-metric weighted contributions (pre mode-2 scaling)
		this->m_LastValMI    = a;
		this->m_LastValNGF   = b;
		this->m_LastValMSE   = c;
		this->m_LastValGD    = d;
		this->m_LastValNC    = e;
		this->m_LastValLabel = f;
		this->m_LastValNMI   = g;

		if (this->m_DerivativeMode == 2)
		{
			// Mode 2: apply the same per-component scale factors cached by
			// GetDerivative() so that value and gradient stay consistent.
			const double total2 = this->m_ScaleMA    * a
			     + this->m_ScaleNGF   * b
			     + this->m_ScaleMSE   * c
			     + this->m_ScaleGD    * d
			     + this->m_ScaleNC    * e
			     + this->m_ScaleLabel * f
			     + this->m_ScaleNMI   * g;
			this->m_LastValTotal = total2;
			return total2;
		}
		const double total0 = a + b + c + d + e + f + g;
		this->m_LastValTotal = total0;
		return total0;
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

	// template <class TFixedImage, class TMovingImage>
	// typename Mplus<TFixedImage, TMovingImage>::MeasureType
	// Mplus<TFixedImage, TMovingImage>::GetCHValue(const ParametersType &parameters) const
	// {
	// 	return static_cast<MeasureType>(m_CH->GetValue(parameters) * this->m_Rho);
	// }

	template <class TFixedImage, class TMovingImage>
	typename Mplus<TFixedImage, TMovingImage>::MeasureType
	Mplus<TFixedImage, TMovingImage>::GetGDValue(const ParametersType &parameters) const
	{
		return static_cast<MeasureType>(m_GD->GetValue(parameters) * this->m_Rho);
	}

	template <class TFixedImage, class TMovingImage>
	typename Mplus<TFixedImage, TMovingImage>::MeasureType
	Mplus<TFixedImage, TMovingImage>::GetNMIValue(const ParametersType &parameters) const
	{
		return static_cast<MeasureType>(m_NMI->GetValue(parameters) * this->m_Sigma);
	}

	template <class TFixedImage, class TMovingImage>
	typename Mplus<TFixedImage, TMovingImage>::MeasureType
	Mplus<TFixedImage, TMovingImage>::GetNCValue(const ParametersType &parameters) const
	{

		return static_cast<MeasureType>(m_NC->GetValue(parameters) * this->m_Yota);
	}

	template <class TFixedImage, class TMovingImage>
	void
	Mplus<TFixedImage, TMovingImage>::GetDerivative(const ParametersType &parameters, DerivativeType &derivative) const
	{
		// ── Collect raw derivatives and merge ────────────────────────────────
		//
		// Three modes, controlled by m_DerivativeMode:
		//
		// ▸ 0 (default, LBFGS-B-safe) – simple weighted sum using the same
		//   weights as GetValue(), preserving value/derivative consistency
		//   required by quasi-Newton optimizers.
		//
		// ▸ 1 (RSGD only) – normalize each component to unit length so the
		//   user‐supplied weights act as pure *direction ratios*; then rescale
		//   the merged gradient to the weighted average of the original norms.
		//   Not suitable for LBFGS-B (breaks Wolfe line search / Hessian approx).
		//
		// ▸ 2 (LBFGS-B-safe) – "main-metric adaptive scaling".  Pick one metric
		//   as reference, scale all other derivatives so their norms match the
		//   main one, then apply user weights.  GetValue uses the same cached
		//   scale factors for value/derivative consistency.

		constexpr double kMinNorm = 1.0e-12;

		// ---------- per‐component raw derivatives --------------------------------
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
			this->GetGDDerivative(parameters, d);
		else
			d.Fill(0.0);

		DerivativeType e;
		e = parameters;
		if (this->m_YotaDerivative != 0.0)
			this->GetNCDerivative(parameters, e);
		else
			e.Fill(0.0);

		DerivativeType f;
		f = parameters;
		if (this->m_LabelKappaDerivative != 0.0 && m_FixedLabelMap && m_MovingLabelMap)
			this->GetKappaDerivative(parameters, f);
		else
			f.Fill(0.0);

		DerivativeType g;
		g = parameters;
		if (this->m_SigmaDerivative != 0.0)
			this->GetNMIDerivative(parameters, g);
		else
			g.Fill(0.0);

		if (this->m_DerivativeMode == 1)
		{
			// ── Mode 1: normalize → merge → rescale (RSGD only) ─────────────
			double normA = this->ComputeDerivativeNorm(a);
			if (normA > kMinNorm)
				for (unsigned int i = 0; i < a.size(); ++i) a[i] /= normA;

			double normB = this->ComputeDerivativeNorm(b);
			if (normB > kMinNorm)
				for (unsigned int i = 0; i < b.size(); ++i) b[i] /= normB;

			double normC = this->ComputeDerivativeNorm(c);
			if (normC > kMinNorm)
				for (unsigned int i = 0; i < c.size(); ++i) c[i] /= normC;

			double normD = this->ComputeDerivativeNorm(d);
			if (normD > kMinNorm)
				for (unsigned int i = 0; i < d.size(); ++i) d[i] /= normD;

			double normE = this->ComputeDerivativeNorm(e);
			if (normE > kMinNorm)
				for (unsigned int i = 0; i < e.size(); ++i) e[i] /= normE;

			double normF = this->ComputeDerivativeNorm(f);
			if (normF > kMinNorm)
				for (unsigned int i = 0; i < f.size(); ++i) f[i] /= normF;

			double normG = this->ComputeDerivativeNorm(g);
			if (normG > kMinNorm)
				for (unsigned int i = 0; i < g.size(); ++i) g[i] /= normG;

			const double wA = std::abs(this->m_AlphaDerivative);
			const double wB = std::abs(this->m_LambdaDerivative);
			const double wC = std::abs(this->m_NuDerivative);
			const double wD = std::abs(this->m_RhoDerivative);
			const double wE = std::abs(this->m_YotaDerivative);
			const double wF = std::abs(this->m_LabelKappaDerivative);
			const double wG = std::abs(this->m_SigmaDerivative);

			derivative = a;
#pragma omp parallel for
			for (long unsigned int p = 0; p < derivative.GetSize(); ++p)
			{
				derivative[p] =
					  wA * a[p]
					+ wB * b[p]
					+ wC * c[p]
					+ wD * d[p]
					+ wE * e[p]
					+ wF * f[p]
					+ wG * g[p];
			}

			// rescale to weighted average of original norms
			const double wSum = wA + wB + wC + wD + wE + wF + wG;
			if (wSum > kMinNorm)
			{
				const double avgNorm = (wA * normA + wB * normB + wC * normC
				                        + wD * normD + wE * normE + wF * normF
				                        + wG * normG) / wSum;
				const double mergedNorm = this->ComputeDerivativeNorm(derivative);
				if (mergedNorm > kMinNorm)
				{
					const double scale = avgNorm / mergedNorm;
#pragma omp parallel for
					for (long unsigned int p = 0; p < derivative.GetSize(); ++p)
						derivative[p] *= scale;
				}
			}
		}
		else if (this->m_DerivativeMode == 2)
		{
			// ── Mode 2: main-metric adaptive scaling (LBFGS-B safe) ─────────
			// Pick a "main" metric, compute its derivative norm, then scale
			// every other component so its norm matches the main one.  The
			// user-supplied weights then act as pure direction ratios without
			// the user needing to compensate for magnitude differences.
			//
			// Effective derivative:
			//   ∇V = wA·sA·∇V_MI + wB·sB·∇V_NGF + wC·sC·∇V_MSE
			//       + wE·sE·∇V_NC + wF·sF·∇V_Label
			// where sX = mainNorm / normX  (sMain = 1 by definition).
			//
			// The same scale factors are cached and re-used in GetValue()
			// so that  V = wA·sA·V_MI + …  , keeping the value/derivative
			// pair consistent for quasi-Newton line searches.

			const double normA = this->ComputeDerivativeNorm(a);
			const double normB = this->ComputeDerivativeNorm(b);
			const double normC = this->ComputeDerivativeNorm(c);
			const double normD = this->ComputeDerivativeNorm(d);
			const double normE = this->ComputeDerivativeNorm(e);
			const double normF = this->ComputeDerivativeNorm(f);
			const double normG = this->ComputeDerivativeNorm(g);

			// Identify the main metric's norm
			double mainNorm = 0.0;
			switch (this->m_MainMetricIndex)
			{
				case 0:  mainNorm = normA; break;  // MI
				case 1:  mainNorm = normB; break;  // NGF
				case 2:  mainNorm = normC; break;  // MSE
				case 3:  mainNorm = normE; break;  // NC
				case 4:  mainNorm = normF; break;  // Label
				case 5:  mainNorm = normD; break;  // GD
				case 6:  mainNorm = normG; break;  // NMI
				default: mainNorm = normA; break;
			}

			// Compute per-component scale factors
			// If mainNorm is tiny everything is near-zero → leave scales at 1
			if (mainNorm > kMinNorm)
			{
				this->m_ScaleMA    = (normA > kMinNorm) ? mainNorm / normA : 1.0;
				this->m_ScaleNGF   = (normB > kMinNorm) ? mainNorm / normB : 1.0;
				this->m_ScaleMSE   = (normC > kMinNorm) ? mainNorm / normC : 1.0;
				this->m_ScaleGD    = (normD > kMinNorm) ? mainNorm / normD : 1.0;
				this->m_ScaleNC    = (normE > kMinNorm) ? mainNorm / normE : 1.0;
				this->m_ScaleLabel = (normF > kMinNorm) ? mainNorm / normF : 1.0;
				this->m_ScaleNMI   = (normG > kMinNorm) ? mainNorm / normG : 1.0;
			}
			else
			{
				this->m_ScaleMA = this->m_ScaleNGF = this->m_ScaleMSE = 1.0;
				this->m_ScaleGD = this->m_ScaleNC = this->m_ScaleLabel = 1.0;
				this->m_ScaleNMI = 1.0;
			}

			derivative = a;
#pragma omp parallel for
			for (long unsigned int p = 0; p < derivative.GetSize(); ++p)
			{
				derivative[p] =
					  this->m_AlphaDerivative      * this->m_ScaleMA    * a[p]
					+ this->m_LambdaDerivative     * this->m_ScaleNGF   * b[p]
					+ this->m_NuDerivative         * this->m_ScaleMSE   * c[p]
					+ this->m_RhoDerivative        * this->m_ScaleGD    * d[p]
					+ this->m_YotaDerivative       * this->m_ScaleNC    * e[p]
					+ this->m_LabelKappaDerivative * this->m_ScaleLabel * f[p]
					+ this->m_SigmaDerivative      * this->m_ScaleNMI   * g[p];
			}
		}
		else
		{
			// ── Mode 0: consistent weighted sum (LBFGS-B, default) ──────────
			// Uses the *derivative* weights, matching GetValue()'s convention:
			//   V = α·V_MI + λ·V_NGF + ν·V_MSE + ρ·V_GD + γ·V_NC + κ·V_Label + σ·V_NMI
			//   ∇V = α·∇V_MI + λ·∇V_NGF + ν·∇V_MSE + ρ·∇V_GD + γ·∇V_NC + κ·∇V_Label + σ·∇V_NMI
			derivative = a;
#pragma omp parallel for
			for (long unsigned int p = 0; p < derivative.GetSize(); ++p)
			{
				derivative[p] =
					  this->m_AlphaDerivative      * a[p]
					+ this->m_LambdaDerivative     * b[p]
					+ this->m_NuDerivative         * c[p]
					+ this->m_RhoDerivative        * d[p]
					+ this->m_YotaDerivative       * e[p]
					+ this->m_LabelKappaDerivative * f[p]
					+ this->m_SigmaDerivative      * g[p];
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
		// Normalization is now handled centrally in GetDerivative().
	}

	// template <class TFixedImage, class TMovingImage>
	// void
	// Mplus<TFixedImage, TMovingImage>::GetCHDerivative(const ParametersType &parameters, DerivativeType &derivative) const
	// {
	// 	m_CH->GetDerivative(parameters, derivative);
	// 	this->NormalizeDerivative(derivative);
	// }

	template <class TFixedImage, class TMovingImage>
	void
	Mplus<TFixedImage, TMovingImage>::GetGDDerivative(const ParametersType &parameters, DerivativeType &derivative) const
	{
		m_GD->GetDerivative(parameters, derivative);
	}

	template <class TFixedImage, class TMovingImage>
	void
	Mplus<TFixedImage, TMovingImage>::GetNMIDerivative(const ParametersType &parameters, DerivativeType &derivative) const
	{
		m_NMI->GetDerivative(parameters, derivative);
	}

	template <class TFixedImage, class TMovingImage>
	void
	Mplus<TFixedImage, TMovingImage>::GetNCDerivative(const ParametersType &parameters, DerivativeType &derivative) const
	{
		m_NC->GetDerivative(parameters, derivative);
	}
	template <class TFixedImage, class TMovingImage>
	void
	Mplus<TFixedImage, TMovingImage>::GetNGFDerivative(const ParametersType &parameters, DerivativeType &derivative) const
	{
		m_NGF->GetDerivative(parameters, derivative);
	}

	template <class TFixedImage, class TMovingImage>
	void
	Mplus<TFixedImage, TMovingImage>::GetMSEDerivative(const ParametersType &parameters, DerivativeType &derivative) const
	{
		m_MSE->GetDerivative(parameters, derivative);
	}

	template <class TFixedImage, class TMovingImage>
	void
	Mplus<TFixedImage, TMovingImage>::GetValueAndDerivative(const ParametersType &parameters, MeasureType &Value, DerivativeType &Derivative) const
	{
		// Optimized: call each sub-metric's GetValueAndDerivative once instead of
		// calling GetValue and GetDerivative separately (which each iterate samples twice).
		constexpr double kMinNorm = 1.0e-12;
		const unsigned int nParams = parameters.GetSize();

		MeasureType rawValA = 0, rawValB = 0, rawValC = 0, rawValD = 0, rawValE = 0, rawValF = 0, rawValG = 0;
		DerivativeType rawDerA(nParams), rawDerB(nParams), rawDerC(nParams),
		               rawDerD(nParams), rawDerE(nParams), rawDerF(nParams), rawDerG(nParams);
		rawDerA.Fill(0.0); rawDerB.Fill(0.0); rawDerC.Fill(0.0);
		rawDerD.Fill(0.0); rawDerE.Fill(0.0); rawDerF.Fill(0.0); rawDerG.Fill(0.0);

		// MA (Mattes MI)
		if (this->m_Alpha != 0.0 && this->m_AlphaDerivative != 0.0)
			m_MA->GetValueAndDerivative(parameters, rawValA, rawDerA);
		else if (this->m_Alpha != 0.0)
			rawValA = m_MA->GetValue(parameters);
		else if (this->m_AlphaDerivative != 0.0)
			m_MA->GetDerivative(parameters, rawDerA);

		// NGF
		if (this->m_Lambda != 0.0 && this->m_LambdaDerivative != 0.0)
			m_NGF->GetValueAndDerivative(parameters, rawValB, rawDerB);
		else if (this->m_Lambda != 0.0)
			rawValB = m_NGF->GetValue(parameters);
		else if (this->m_LambdaDerivative != 0.0)
			m_NGF->GetDerivative(parameters, rawDerB);

		// MSE
		if (this->m_Nu != 0.0 && this->m_NuDerivative != 0.0)
			m_MSE->GetValueAndDerivative(parameters, rawValC, rawDerC);
		else if (this->m_Nu != 0.0)
			rawValC = m_MSE->GetValue(parameters);
		else if (this->m_NuDerivative != 0.0)
			m_MSE->GetDerivative(parameters, rawDerC);

		// GD (Gradient Difference)
		if (this->m_Rho != 0.0 && this->m_RhoDerivative != 0.0)
			m_GD->GetValueAndDerivative(parameters, rawValD, rawDerD);
		else if (this->m_Rho != 0.0)
			rawValD = m_GD->GetValue(parameters);
		else if (this->m_RhoDerivative != 0.0)
			m_GD->GetDerivative(parameters, rawDerD);

		// NC
		if (this->m_Yota != 0.0 && this->m_YotaDerivative != 0.0)
			m_NC->GetValueAndDerivative(parameters, rawValE, rawDerE);
		else if (this->m_Yota != 0.0)
			rawValE = m_NC->GetValue(parameters);
		else if (this->m_YotaDerivative != 0.0)
			m_NC->GetDerivative(parameters, rawDerE);

		// Kappa (label metric)
		if ((this->m_LabelKappa != 0.0 || this->m_LabelKappaDerivative != 0.0)
		    && m_FixedLabelMap && m_MovingLabelMap) {
			if (this->m_LabelKappa != 0.0 && this->m_LabelKappaDerivative != 0.0)
				this->GetKappaValueAndDerivative(parameters, rawValF, rawDerF);
			else if (this->m_LabelKappa != 0.0)
				rawValF = this->GetKappaValue(parameters);
			else
				this->GetKappaDerivative(parameters, rawDerF);
		}

		// NMI (Normalized Mutual Information)
		if (this->m_Sigma != 0.0 && this->m_SigmaDerivative != 0.0)
			m_NMI->GetValueAndDerivative(parameters, rawValG, rawDerG);
		else if (this->m_Sigma != 0.0)
			rawValG = m_NMI->GetValue(parameters);
		else if (this->m_SigmaDerivative != 0.0)
			m_NMI->GetDerivative(parameters, rawDerG);

		// ── Cache weighted per-sub-metric contributions ───────────────────
		this->m_LastValMI    = this->m_Alpha  * rawValA;
		this->m_LastValNGF   = this->m_Lambda * rawValB;
		this->m_LastValMSE   = this->m_Nu     * rawValC;
		this->m_LastValGD    = this->m_Rho    * rawValD;
		this->m_LastValNC    = this->m_Yota   * rawValE;
		this->m_LastValLabel = this->m_LabelKappa * rawValF;
		this->m_LastValNMI   = this->m_Sigma  * rawValG;
		this->m_LastValTotal = this->m_LastValMI + this->m_LastValNGF
		                     + this->m_LastValMSE + this->m_LastValGD
		                     + this->m_LastValNC  + this->m_LastValLabel
		                     + this->m_LastValNMI;

		// ── Combine derivatives and value ────────────────────────────────
		if (this->m_DerivativeMode == 1)
		{
			// Mode 1: normalize + merge + rescale (RSGD only)
			double normA = this->ComputeDerivativeNorm(rawDerA);
			if (normA > kMinNorm) for (unsigned i = 0; i < nParams; ++i) rawDerA[i] /= normA;
			double normB = this->ComputeDerivativeNorm(rawDerB);
			if (normB > kMinNorm) for (unsigned i = 0; i < nParams; ++i) rawDerB[i] /= normB;
			double normC = this->ComputeDerivativeNorm(rawDerC);
			if (normC > kMinNorm) for (unsigned i = 0; i < nParams; ++i) rawDerC[i] /= normC;
			double normD = this->ComputeDerivativeNorm(rawDerD);
			if (normD > kMinNorm) for (unsigned i = 0; i < nParams; ++i) rawDerD[i] /= normD;
			double normE = this->ComputeDerivativeNorm(rawDerE);
			if (normE > kMinNorm) for (unsigned i = 0; i < nParams; ++i) rawDerE[i] /= normE;
			double normF = this->ComputeDerivativeNorm(rawDerF);
			if (normF > kMinNorm) for (unsigned i = 0; i < nParams; ++i) rawDerF[i] /= normF;
			double normG = this->ComputeDerivativeNorm(rawDerG);
			if (normG > kMinNorm) for (unsigned i = 0; i < nParams; ++i) rawDerG[i] /= normG;

			const double wA = std::abs(this->m_AlphaDerivative);
			const double wB = std::abs(this->m_LambdaDerivative);
			const double wC = std::abs(this->m_NuDerivative);
			const double wD = std::abs(this->m_RhoDerivative);
			const double wE = std::abs(this->m_YotaDerivative);
			const double wF = std::abs(this->m_LabelKappaDerivative);
			const double wG = std::abs(this->m_SigmaDerivative);

			Derivative.SetSize(nParams);
#pragma omp parallel for
			for (long unsigned int p = 0; p < nParams; ++p)
				Derivative[p] = wA*rawDerA[p] + wB*rawDerB[p] + wC*rawDerC[p]
				              + wD*rawDerD[p] + wE*rawDerE[p] + wF*rawDerF[p]
				              + wG*rawDerG[p];

			const double wSum = wA + wB + wC + wD + wE + wF + wG;
			if (wSum > kMinNorm) {
				const double avgNorm = (wA*normA + wB*normB + wC*normC + wD*normD
				                        + wE*normE + wF*normF + wG*normG) / wSum;
				const double mergedNorm = this->ComputeDerivativeNorm(Derivative);
				if (mergedNorm > kMinNorm) {
					const double scale = avgNorm / mergedNorm;
#pragma omp parallel for
					for (long unsigned int p = 0; p < nParams; ++p)
						Derivative[p] *= scale;
				}
			}

			Value = this->m_Alpha * rawValA + this->m_Lambda * rawValB
			      + this->m_Nu * rawValC + this->m_Rho * rawValD
			      + this->m_Yota * rawValE + this->m_LabelKappa * rawValF
			      + this->m_Sigma * rawValG;
		}
		else if (this->m_DerivativeMode == 2)
		{
			// Mode 2: main-metric adaptive scaling
			const double normA = this->ComputeDerivativeNorm(rawDerA);
			const double normB = this->ComputeDerivativeNorm(rawDerB);
			const double normC = this->ComputeDerivativeNorm(rawDerC);
			const double normD = this->ComputeDerivativeNorm(rawDerD);
			const double normE = this->ComputeDerivativeNorm(rawDerE);
			const double normF = this->ComputeDerivativeNorm(rawDerF);
			const double normG = this->ComputeDerivativeNorm(rawDerG);

			double mainNorm = 0.0;
			switch (this->m_MainMetricIndex) {
				case 0: mainNorm = normA; break;
				case 1: mainNorm = normB; break;
				case 2: mainNorm = normC; break;
				case 3: mainNorm = normE; break;
				case 4: mainNorm = normF; break;
				case 5: mainNorm = normD; break;
				case 6: mainNorm = normG; break;
				default: mainNorm = normA; break;
			}
			if (mainNorm > kMinNorm) {
				this->m_ScaleMA    = (normA > kMinNorm) ? mainNorm / normA : 1.0;
				this->m_ScaleNGF   = (normB > kMinNorm) ? mainNorm / normB : 1.0;
				this->m_ScaleMSE   = (normC > kMinNorm) ? mainNorm / normC : 1.0;
				this->m_ScaleGD    = (normD > kMinNorm) ? mainNorm / normD : 1.0;
				this->m_ScaleNC    = (normE > kMinNorm) ? mainNorm / normE : 1.0;
				this->m_ScaleLabel = (normF > kMinNorm) ? mainNorm / normF : 1.0;
				this->m_ScaleNMI   = (normG > kMinNorm) ? mainNorm / normG : 1.0;
			} else {
				this->m_ScaleMA = this->m_ScaleNGF = this->m_ScaleMSE = 1.0;
				this->m_ScaleGD = this->m_ScaleNC = this->m_ScaleLabel = 1.0;
				this->m_ScaleNMI = 1.0;
			}

			Derivative.SetSize(nParams);
#pragma omp parallel for
			for (long unsigned int p = 0; p < nParams; ++p)
				Derivative[p] = this->m_AlphaDerivative  * this->m_ScaleMA    * rawDerA[p]
				              + this->m_LambdaDerivative * this->m_ScaleNGF   * rawDerB[p]
				              + this->m_NuDerivative     * this->m_ScaleMSE   * rawDerC[p]
				              + this->m_RhoDerivative    * this->m_ScaleGD    * rawDerD[p]
				              + this->m_YotaDerivative   * this->m_ScaleNC    * rawDerE[p]
				              + this->m_LabelKappaDerivative * this->m_ScaleLabel * rawDerF[p]
				              + this->m_SigmaDerivative  * this->m_ScaleNMI   * rawDerG[p];

			Value = this->m_ScaleMA  * this->m_Alpha  * rawValA
			      + this->m_ScaleNGF * this->m_Lambda * rawValB
			      + this->m_ScaleMSE * this->m_Nu     * rawValC
			      + this->m_ScaleGD  * this->m_Rho    * rawValD
			      + this->m_ScaleNC  * this->m_Yota   * rawValE
			      + this->m_ScaleLabel * this->m_LabelKappa * rawValF
			      + this->m_ScaleNMI * this->m_Sigma  * rawValG;
		}
		else
		{
			// Mode 0: consistent weighted sum
			Derivative.SetSize(nParams);
#pragma omp parallel for
			for (long unsigned int p = 0; p < nParams; ++p)
				Derivative[p] = this->m_AlphaDerivative  * rawDerA[p]
				              + this->m_LambdaDerivative * rawDerB[p]
				              + this->m_NuDerivative     * rawDerC[p]
				              + this->m_RhoDerivative    * rawDerD[p]
				              + this->m_YotaDerivative   * rawDerE[p]
				              + this->m_LabelKappaDerivative * rawDerF[p]
				              + this->m_SigmaDerivative  * rawDerG[p];

			Value = this->m_Alpha * rawValA + this->m_Lambda * rawValB
			      + this->m_Nu * rawValC + this->m_Rho * rawValD
			      + this->m_Yota * rawValE + this->m_LabelKappa * rawValF
			      + this->m_Sigma * rawValG;
		}
	}

	template <class TImageType, class TMovingImage>
	void
	Mplus<TImageType, TMovingImage>::
		PrintSelf(std::ostream &os, Indent indent) const
	{
		Superclass::PrintSelf(os, indent);
	}

	// ── Label-map / ROI metric implementations ───────────────────────────────

	template <class TFixedImage, class TMovingImage>
	typename TFixedImage::Pointer
	Mplus<TFixedImage, TMovingImage>::BuildBinaryFromLabel(
	    LabelImageConstPointer                      labelMap,
	    const typename TFixedImage::SizeType      & refSize,
	    const typename TFixedImage::SpacingType   & refSpacing,
	    const typename TFixedImage::PointType     & refOrigin,
	    const typename TFixedImage::DirectionType & refDirection,
	    LabelPixelType                              L) const
	{
	    // Cast label map to float
	    using CastType = itk::CastImageFilter<LabelImageType, TFixedImage>;
	    auto caster = CastType::New();
	    caster->SetInput(labelMap);

	    // Resample onto the reference geometry with nearest-neighbor
	    using ResampleType = itk::ResampleImageFilter<TFixedImage, TFixedImage>;
	    using NNType       = itk::NearestNeighborInterpolateImageFunction<TFixedImage, double>;
	    auto resampler = ResampleType::New();
	    resampler->SetInput(caster->GetOutput());
	    resampler->SetInterpolator(NNType::New());
	    resampler->SetSize(refSize);
	    resampler->SetOutputSpacing(refSpacing);
	    resampler->SetOutputOrigin(refOrigin);
	    resampler->SetOutputDirection(refDirection);
	    resampler->SetDefaultPixelValue(0.0f);

	    // Threshold to binary: 1.0 where pixel is within ±0.5 of L
	    using ThreshType = itk::BinaryThresholdImageFilter<TFixedImage, TFixedImage>;
	    auto thresh = ThreshType::New();
	    thresh->SetInput(resampler->GetOutput());
	    thresh->SetLowerThreshold(static_cast<float>(L) - 0.5f);
	    thresh->SetUpperThreshold(static_cast<float>(L) + 0.5f);
	    thresh->SetInsideValue(1.0f);
	    thresh->SetOutsideValue(0.0f);
	    thresh->Update();

	    typename TFixedImage::Pointer out = thresh->GetOutput();
	    out->DisconnectPipeline();
	    return out;
	}

	template <class TFixedImage, class TMovingImage>
	typename TFixedImage::Pointer
	Mplus<TFixedImage, TMovingImage>::ComputeSignedDist(
	    const typename TFixedImage::Pointer & binaryImage) const
	{
	    using UCImage     = itk::Image<unsigned char, TFixedImage::ImageDimension>;
	    using CastToUC    = itk::CastImageFilter<TFixedImage, UCImage>;
	    using DistFilter  = itk::SignedMaurerDistanceMapImageFilter<UCImage, TFixedImage>;

	    auto castUC = CastToUC::New();
	    castUC->SetInput(binaryImage);

	    auto dist = DistFilter::New();
	    dist->SetInput(castUC->GetOutput());
	    dist->SetSquaredDistance(false);
	    dist->SetUseImageSpacing(true);
	    dist->SetInsideIsPositive(false);   // negative inside, positive outside
	    dist->Update();

	    typename TFixedImage::Pointer out = dist->GetOutput();
	    out->DisconnectPipeline();
	    return out;
	}

	template <class TFixedImage, class TMovingImage>
	typename Mplus<TFixedImage, TMovingImage>::GradientImagePointer
	Mplus<TFixedImage, TMovingImage>::ComputeGradient(
	    const typename TFixedImage::Pointer & image) const
	{
	    using GradFilter = itk::GradientRecursiveGaussianImageFilter<TFixedImage, GradientImageType>;
	    auto grad = GradFilter::New();
	    grad->SetInput(image);
	    grad->SetSigma(1.0);   // mild smoothing; distance maps are already smooth
	    grad->Update();

	    GradientImagePointer out = grad->GetOutput();
	    out->DisconnectPipeline();
	    return out;
	}

	template <class TFixedImage, class TMovingImage>
	void
	Mplus<TFixedImage, TMovingImage>::InitializeLabelMetric()
	{
	    m_FixedDistMaps.clear();
	    m_MovingDistMaps.clear();
	    m_MovingDistGradMaps.clear();
	    m_MovingDistInterps.clear();
	    m_MovingDistGradInterps.clear();
	    m_LabelValues.clear();
	    m_LastDice.clear();

	    if (!m_FixedLabelMap || !m_MovingLabelMap) return;
	    if (m_LabelKappa == 0.0 && m_LabelKappaDerivative == 0.0) return;

	    // Collect unique non-zero labels from both maps
	    std::set<LabelPixelType> labelSet;
	    {
	        itk::ImageRegionConstIterator<LabelImageType> it(
	            m_FixedLabelMap, m_FixedLabelMap->GetLargestPossibleRegion());
	        for (; !it.IsAtEnd(); ++it)
	            if (it.Get() != 0) labelSet.insert(it.Get());
	    }
	    {
	        itk::ImageRegionConstIterator<LabelImageType> it(
	            m_MovingLabelMap, m_MovingLabelMap->GetLargestPossibleRegion());
	        for (; !it.IsAtEnd(); ++it)
	            if (it.Get() != 0) labelSet.insert(it.Get());
	    }
	    m_LabelValues.assign(labelSet.begin(), labelSet.end());

	    if (m_LabelValues.empty()) return;

	    // Fixed-image geometry (for fixed distance maps)
	    const auto & fxSize = this->m_FixedImage->GetLargestPossibleRegion().GetSize();
	    const auto & fxSpc  = this->m_FixedImage->GetSpacing();
	    const auto & fxOrg  = this->m_FixedImage->GetOrigin();
	    const auto & fxDir  = this->m_FixedImage->GetDirection();

	    // Moving-image geometry (for moving distance maps, kept in moving space)
	    const auto & mvSize = this->m_MovingImage->GetLargestPossibleRegion().GetSize();
	    const auto & mvSpc  = this->m_MovingImage->GetSpacing();
	    const auto & mvOrg  = this->m_MovingImage->GetOrigin();
	    const auto & mvDir  = this->m_MovingImage->GetDirection();

	    std::cout << "[LabelMetric] Initializing distance maps for "
	              << m_LabelValues.size() << " label(s)..." << std::endl;

	    for (auto L : m_LabelValues)
	    {
	        // Fixed binary + distance map (in fixed image space)
	        auto fixedBin  = this->BuildBinaryFromLabel(
	            m_FixedLabelMap,  fxSize, fxSpc, fxOrg, fxDir, L);
	        m_FixedDistMaps[L] = this->ComputeSignedDist(fixedBin);

	        // Moving binary + distance map (in moving image space)
	        auto movingBin = this->BuildBinaryFromLabel(
	            m_MovingLabelMap, mvSize, mvSpc, mvOrg, mvDir, L);
	        m_MovingDistMaps[L]     = this->ComputeSignedDist(movingBin);
	        m_MovingDistGradMaps[L] = this->ComputeGradient(m_MovingDistMaps[L]);

	        // Interpolators for moving distance map and its gradient
	        auto di = DistInterpType::New();
	        di->SetInputImage(m_MovingDistMaps[L]);
	        m_MovingDistInterps[L] = di;

	        auto gi = GradInterpType::New();
	        gi->SetInputImage(m_MovingDistGradMaps[L]);
	        m_MovingDistGradInterps[L] = gi;

	        std::cout << "[LabelMetric]   Label " << static_cast<int>(L) << " done." << std::endl;
	    }
	}

	template <class TFixedImage, class TMovingImage>
	typename Mplus<TFixedImage, TMovingImage>::MeasureType
	Mplus<TFixedImage, TMovingImage>::GetKappaValue(const ParametersType & parameters) const
	{
	    if (m_LabelValues.empty()) return 0.0;

	    // Ensure transform has the current parameters
	    this->m_Transform->SetParameters(parameters);

	    const unsigned int nLabels = static_cast<unsigned int>(m_LabelValues.size());
	    const double flatWeight    = 1.0 / static_cast<double>(nLabels);

	    // Stride to hit approximately m_LabelNumberOfSamples sample points
	    const unsigned long totalPix = this->m_FixedImage
	        ->GetLargestPossibleRegion().GetNumberOfPixels();
	    unsigned int stride = 1;
	    if (m_LabelNumberOfSamples < totalPix)
	    {
	        stride = static_cast<unsigned int>(
	            std::ceil(std::pow(static_cast<double>(totalPix) /
	                               static_cast<double>(m_LabelNumberOfSamples),
	                               1.0 / TFixedImage::ImageDimension)));
	        if (stride < 1) stride = 1;
	    }

	    double totalValue = 0.0;

	    for (auto L : m_LabelValues)
	    {
	        const double kappaL = (m_LabelKappaWeights.count(L) > 0)
	            ? m_LabelKappaWeights.at(L) : flatWeight;
	        if (kappaL == 0.0) continue;

	        const auto & fixedDistMap    = m_FixedDistMaps.at(L);
	        const auto & movingDistInterp= m_MovingDistInterps.at(L);

	        double sumSqDiff = 0.0;
	        long long countF = 0, countM = 0, countIntersect = 0;
	        unsigned int n = 0, pixIdx = 0;

	        itk::ImageRegionConstIteratorWithIndex<TFixedImage> it(
	            fixedDistMap, fixedDistMap->GetLargestPossibleRegion());

	        for (it.GoToBegin(); !it.IsAtEnd(); ++it, ++pixIdx)
	        {
	            if (pixIdx % stride != 0) continue;

	            typename TFixedImage::PointType fixedPt;
	            fixedDistMap->TransformIndexToPhysicalPoint(it.GetIndex(), fixedPt);
	            const auto movingPt = this->m_Transform->TransformPoint(fixedPt);

	            if (!movingDistInterp->IsInsideBuffer(movingPt)) continue;

	            const double dFixed  = static_cast<double>(it.Get());
	            const double dMoving = movingDistInterp->Evaluate(movingPt);
	            const double diff    = dFixed - dMoving;
	            sumSqDiff += diff * diff;
	            ++n;

	            // Also accumulate binary Dice stats (dFixed<0 → inside label)
	            bool inF = (dFixed  <= 0.0);
	            bool inM = (dMoving <= 0.0);
	            if (inF) ++countF;
	            if (inM) ++countM;
	            if (inF && inM) ++countIntersect;
	        }

	        if (n == 0) continue;

	        totalValue += kappaL * sumSqDiff / static_cast<double>(n);

	        // Cache Dice for monitoring
	        double dice = 0.0;
	        if (countF + countM > 0)
	            dice = 2.0 * countIntersect / static_cast<double>(countF + countM);
	        m_LastDice[L] = dice;
	    }

	    return static_cast<MeasureType>(m_LabelKappa * totalValue);
	}

	template <class TFixedImage, class TMovingImage>
	void
	Mplus<TFixedImage, TMovingImage>::GetKappaDerivative(
	    const ParametersType & parameters, DerivativeType & derivative) const
	{
	    const unsigned int nParams = static_cast<unsigned int>(parameters.GetSize());
	    derivative.SetSize(nParams);
	    derivative.Fill(0.0);

	    if (m_LabelValues.empty()) return;

	    this->m_Transform->SetParameters(parameters);

	    const unsigned int nLabels   = static_cast<unsigned int>(m_LabelValues.size());
	    const double       flatWeight= 1.0 / static_cast<double>(nLabels);

	    const unsigned long totalPix = this->m_FixedImage
	        ->GetLargestPossibleRegion().GetNumberOfPixels();
	    unsigned int stride = 1;
	    if (m_LabelNumberOfSamples < totalPix)
	    {
	        stride = static_cast<unsigned int>(
	            std::ceil(std::pow(static_cast<double>(totalPix) /
	                               static_cast<double>(m_LabelNumberOfSamples),
	                               1.0 / TFixedImage::ImageDimension)));
	        if (stride < 1) stride = 1;
	    }

	    for (auto L : m_LabelValues)
	    {
	        const double kappaL = (m_LabelKappaDerivWeights.count(L) > 0)
	            ? m_LabelKappaDerivWeights.at(L) : flatWeight;
	        if (kappaL == 0.0) continue;

	        const auto & fixedDistMap    = m_FixedDistMaps.at(L);
	        const auto & movingDistInterp= m_MovingDistInterps.at(L);
	        const auto & gradInterp      = m_MovingDistGradInterps.at(L);

	        // Accumulate per-label contribution into a temporary vector
	        std::vector<double> localDeriv(nParams, 0.0);
	        unsigned int n = 0, pixIdx = 0;

	        itk::ImageRegionConstIteratorWithIndex<TFixedImage> it(
	            fixedDistMap, fixedDistMap->GetLargestPossibleRegion());

	        for (it.GoToBegin(); !it.IsAtEnd(); ++it, ++pixIdx)
	        {
	            if (pixIdx % stride != 0) continue;

	            typename TFixedImage::PointType fixedPt;
	            fixedDistMap->TransformIndexToPhysicalPoint(it.GetIndex(), fixedPt);
	            const auto movingPt = this->m_Transform->TransformPoint(fixedPt);

	            if (!movingDistInterp->IsInsideBuffer(movingPt) ||
	                !gradInterp->IsInsideBuffer(movingPt)) continue;

	            const double dFixed   = static_cast<double>(it.Get());
	            const double dMoving  = movingDistInterp->Evaluate(movingPt);
	            const double residual = dFixed - dMoving;

	            // Gradient of moving distance map at the transformed point
	            const auto gradVec = gradInterp->Evaluate(movingPt);

	            // Transform Jacobian at the fixed point (dim × nParams)
	            TransformJacobianType jac(TFixedImage::ImageDimension, nParams);
	            this->m_Transform->ComputeJacobianWithRespectToParameters(fixedPt, jac);

	            // dE/dθ_j = -2 * residual * (∇dM · J[:,j])
	            for (unsigned int j = 0; j < nParams; ++j)
	            {
	                double dot = 0.0;
	                for (unsigned int d = 0; d < TFixedImage::ImageDimension; ++d)
	                    dot += static_cast<double>(gradVec[d]) * jac(d, j);
	                localDeriv[j] += -2.0 * residual * dot;
	            }
	            ++n;
	        }

	        if (n == 0) continue;

	        // Scale by m_LabelKappa * kappaL / N — must include m_LabelKappa to
	        // match GetKappaValue which returns m_LabelKappa * totalValue.
	        const double scale = m_LabelKappa * kappaL / static_cast<double>(n);
	        #pragma omp parallel for
	        for (unsigned int j = 0; j < nParams; ++j)
	            derivative[j] += scale * localDeriv[j];
	    }
	}

	template <class TFixedImage, class TMovingImage>
	void
	Mplus<TFixedImage, TMovingImage>::GetKappaValueAndDerivative(
	    const ParametersType & parameters, MeasureType & value,
	    DerivativeType & derivative) const
	{
	    const unsigned int nParams = static_cast<unsigned int>(parameters.GetSize());
	    derivative.SetSize(nParams);
	    derivative.Fill(0.0);
	    value = 0.0;

	    if (m_LabelValues.empty()) return;

	    this->m_Transform->SetParameters(parameters);

	    const unsigned int nLabels    = static_cast<unsigned int>(m_LabelValues.size());
	    const double       flatWeight = 1.0 / static_cast<double>(nLabels);

	    const unsigned long totalPix = this->m_FixedImage
	        ->GetLargestPossibleRegion().GetNumberOfPixels();
	    unsigned int stride = 1;
	    if (m_LabelNumberOfSamples < totalPix) {
	        stride = static_cast<unsigned int>(
	            std::ceil(std::pow(static_cast<double>(totalPix) /
	                               static_cast<double>(m_LabelNumberOfSamples),
	                               1.0 / TFixedImage::ImageDimension)));
	        if (stride < 1) stride = 1;
	    }

	    double totalValue = 0.0;

	    for (auto L : m_LabelValues)
	    {
	        const double kappaV = (m_LabelKappaWeights.count(L) > 0)
	            ? m_LabelKappaWeights.at(L) : flatWeight;
	        const double kappaD = (m_LabelKappaDerivWeights.count(L) > 0)
	            ? m_LabelKappaDerivWeights.at(L) : flatWeight;

	        const auto & fixedDistMap     = m_FixedDistMaps.at(L);
	        const auto & movingDistInterp = m_MovingDistInterps.at(L);
	        const auto & gradInterp       = m_MovingDistGradInterps.at(L);

	        double sumSqDiff = 0.0;
	        long long countF = 0, countM = 0, countIntersect = 0;
	        std::vector<double> localDeriv(nParams, 0.0);
	        unsigned int n = 0, pixIdx = 0;

	        itk::ImageRegionConstIteratorWithIndex<TFixedImage> it(
	            fixedDistMap, fixedDistMap->GetLargestPossibleRegion());

	        for (it.GoToBegin(); !it.IsAtEnd(); ++it, ++pixIdx)
	        {
	            if (pixIdx % stride != 0) continue;

	            typename TFixedImage::PointType fixedPt;
	            fixedDistMap->TransformIndexToPhysicalPoint(it.GetIndex(), fixedPt);
	            const auto movingPt = this->m_Transform->TransformPoint(fixedPt);

	            if (!movingDistInterp->IsInsideBuffer(movingPt)) continue;

	            const double dFixed  = static_cast<double>(it.Get());
	            const double dMoving = movingDistInterp->Evaluate(movingPt);
	            const double residual = dFixed - dMoving;

	            // Value accumulation
	            sumSqDiff += residual * residual;

	            // Derivative accumulation
	            if (kappaD != 0.0 && gradInterp->IsInsideBuffer(movingPt)) {
	                const auto gradVec = gradInterp->Evaluate(movingPt);
	                TransformJacobianType jac(TFixedImage::ImageDimension, nParams);
	                this->m_Transform->ComputeJacobianWithRespectToParameters(fixedPt, jac);
	                for (unsigned int j = 0; j < nParams; ++j) {
	                    double dot = 0.0;
	                    for (unsigned int d = 0; d < TFixedImage::ImageDimension; ++d)
	                        dot += static_cast<double>(gradVec[d]) * jac(d, j);
	                    localDeriv[j] += -2.0 * residual * dot;
	                }
	            }

	            ++n;

	            bool inF = (dFixed  <= 0.0);
	            bool inM = (dMoving <= 0.0);
	            if (inF) ++countF;
	            if (inM) ++countM;
	            if (inF && inM) ++countIntersect;
	        }

	        if (n == 0) continue;

	        totalValue += kappaV * sumSqDiff / static_cast<double>(n);

        // Include m_LabelKappa to keep value/derivative consistent.
        const double scaleD = m_LabelKappa * kappaD / static_cast<double>(n);
        #pragma omp parallel for
        for (unsigned int j = 0; j < nParams; ++j)
            derivative[j] += scaleD * localDeriv[j];

        double dice = 0.0;
        if (countF + countM > 0)
            dice = 2.0 * countIntersect / static_cast<double>(countF + countM);
        m_LastDice[L] = dice;
    }

	    value = static_cast<MeasureType>(m_LabelKappa * totalValue);
	}

} // end namespace itk

#endif
