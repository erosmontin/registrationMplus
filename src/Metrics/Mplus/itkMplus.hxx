#ifndef __itkMplus_hxx
#define __itkMplus_hxx

#include "itkMplus.h"
#include "../NGF/NGFImageMetric/NGFImageToImageMetric/Code/itkNGFMetricKernel.h"
#include "itkMinimumMaximumImageCalculator.h"

// <<< add for auto-eta
#include "itkGradientMagnitudeImageFilter.h"
#include "itkImageRegionIterator.h"
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>
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
#include <iterator>
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
			m_LabelDistanceMax     = 20.0;
			m_LabelUseNarrowBand   = false;
			m_LabelNarrowBandWidth = 5.0;
			m_LabelUseHuber        = false;
			m_LabelHuberDelta      = 0.25;

		m_NormalizeMSE = false;
		m_MSEIntensityRangeSquared = 1.0;
		m_NormalizeGD = false;
		m_GDNormalizationFactor = 1.0;
		m_NCDegenerate = false;

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

			// Cache fixed-image intensity range squared for NormalizeMSE
			if (this->m_NormalizeMSE)
			{
				typedef itk::MinimumMaximumImageCalculator<TFixedImage> MinMaxCalcType;
				typename MinMaxCalcType::Pointer calc = MinMaxCalcType::New();
				calc->SetImage(this->m_FixedImage);
				calc->Compute();
				const double range = static_cast<double>(calc->GetMaximum())
				                   - static_cast<double>(calc->GetMinimum());
				m_MSEIntensityRangeSquared = (range > 1e-6) ? (range * range) : 1.0;
			}
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

			// Safety: NC divides by moving-image std-dev.  If the moving image is
			// constant (e.g. pre-warped background) the denominator is zero and the
			// derivative buffer fills with NaN → access violation.  Check variance
			// over 5% of voxels and auto-disable if the image is flat.
			{
				constexpr double kMinStd = 1e-6;
				constexpr double kSampleFrac = 0.05;
				const unsigned long totalPix = this->m_MovingImage->GetLargestPossibleRegion().GetNumberOfPixels();
				const size_t maxSamp = std::max(static_cast<size_t>(5000),
				                                static_cast<size_t>(totalPix * kSampleFrac));
				double sum = 0.0, sum2 = 0.0;
				size_t cnt = 0;
				itk::ImageRegionConstIterator<MovingImageType> it(
					this->m_MovingImage, this->m_MovingImage->GetLargestPossibleRegion());
				for (; !it.IsAtEnd() && cnt < maxSamp; ++it, ++cnt)
				{
					double v = static_cast<double>(it.Get());
					sum  += v;
					sum2 += v * v;
				}
				double mean = sum / cnt;
				double var  = sum2 / cnt - mean * mean;
				if (std::sqrt(var) < kMinStd)
				{
					std::cout << "  [SAFETY] NC: moving image is constant (std="
					          << std::scientific << std::setprecision(2)
					          << std::sqrt(var) << std::defaultfloat
					          << ") — auto-disabling NC to prevent division-by-zero crash." << std::endl;
					this->m_NCDegenerate = true;
				}
			}
		}

		// ── Gradient Difference (GD) ──────────────────────────────────────────
		// NOTE: GD's Initialize() calls a ResampleImageFilter that calls
		// SetInputImage() on the interpolator, leaving it in a modified state.
		// Give GD its own private interpolator clone to avoid corrupting the
		// shared interpolator used by MA, NGF, MSE, and NC.
		if (this->m_Rho != 0.0 || this->m_RhoDerivative != 0.0)
		{
			typename InterpolatorType::Pointer gdInterp =
				dynamic_cast<InterpolatorType*>(this->GetInterpolator()->CreateAnother().GetPointer());
			if (!gdInterp) gdInterp = this->GetInterpolator();
			m_GD = GDType::New();
			m_GD->SetFixedImage(this->GetFixedImage());
			m_GD->SetMovingImage(this->GetMovingImage());
			m_GD->SetTransform(this->GetTransform());
			m_GD->SetInterpolator(gdInterp);
			m_GD->SetFixedImageRegion(overlap);
			m_GD->SetDerivativeDelta(0.001);
			m_GD->SetNumberOfThreads(this->GetNumberOfThreads());
			m_GD->SetUseCachingOfBSplineWeights(false);
			m_GD->Initialize();

			// Cache 1/N_overlap_voxels for NormalizeGD.  GD iterates the full
			// fixed-image (overlap) region and accumulates per-voxel terms in
			// [0,1], so dividing by the voxel count converts the sum into a
			// mean comparable in magnitude to MI/NGF/NC.
			if (this->m_NormalizeGD)
			{
				const double n = static_cast<double>(overlap.GetNumberOfPixels());
				m_GDNormalizationFactor = (n > 0.0) ? (1.0 / n) : 1.0;
			}
		}

		// ── Normalized Mutual Information (NMI) ──────────────────────────────
		// Give NMI its own interpolator clone for the same reason as GD above.
		if (this->m_Sigma != 0.0 || this->m_SigmaDerivative != 0.0)
		{
			typename InterpolatorType::Pointer nmiInterp =
				dynamic_cast<InterpolatorType*>(this->GetInterpolator()->CreateAnother().GetPointer());
			if (!nmiInterp) nmiInterp = this->GetInterpolator();
			m_NMI = NMIType::New();
			m_NMI->SetFixedImage(this->GetFixedImage());
			m_NMI->SetMovingImage(this->GetMovingImage());
			m_NMI->SetTransform(this->GetTransform());
			m_NMI->SetInterpolator(nmiInterp);
			m_NMI->SetFixedImageRegion(overlap);
			m_NMI->SetNumberOfThreads(this->GetNumberOfThreads());
			m_NMI->SetUseCachingOfBSplineWeights(false);
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
				// Foreground-aware estimator:
				//   1) determine a foreground intensity threshold from the
				//      intensity range (5% above the min);
				//   2) restrict the gradient-magnitude population to voxels
				//      whose intensity is above that threshold;
				//   3) use the MEDIAN of those gradients as η (robust to
				//      outliers and to large background regions caused by
				//      pre-warping / FOV mismatch).
				constexpr double kPercentile      = 0.50;
				constexpr double kFgIntensityFrac = 0.05;
				constexpr double kSampleFraction  = 0.05;

				bool etaMAtFloor = false;

				// ── Fixed image η ──────────────────────────────────────────
				{
					using GradFilterType = itk::GradientMagnitudeImageFilter<FixedImageType, FixedImageType>;
					typename GradFilterType::Pointer gradFilter = GradFilterType::New();
					gradFilter->SetInput(this->m_FixedImage);
					gradFilter->Update();

					const unsigned long totalPix = this->m_FixedImage->GetLargestPossibleRegion().GetNumberOfPixels();
					const size_t maxSamples = std::max(static_cast<size_t>(20000),
					                                   static_cast<size_t>(totalPix * kSampleFraction));

					double minVal =  std::numeric_limits<double>::max();
					double maxVal = -std::numeric_limits<double>::max();
					{
						itk::ImageRegionConstIterator<FixedImageType> it(
							this->m_FixedImage, this->m_FixedImage->GetLargestPossibleRegion());
						size_t cnt = 0;
						for (; !it.IsAtEnd() && cnt < maxSamples; ++it, ++cnt)
						{
							double v = static_cast<double>(it.Get());
							if (v < minVal) minVal = v;
							if (v > maxVal) maxVal = v;
						}
					}
					const double range = maxVal - minVal;
					const double fgThreshold = (range > 0.0)
					    ? (minVal + kFgIntensityFrac * range) : minVal;

					std::vector<double> mags;
					mags.reserve(maxSamples);
					{
						itk::ImageRegionConstIterator<FixedImageType> itImg(
							this->m_FixedImage, this->m_FixedImage->GetLargestPossibleRegion());
						itk::ImageRegionConstIterator<FixedImageType> itGrad(
							gradFilter->GetOutput(),
							gradFilter->GetOutput()->GetLargestPossibleRegion());
						size_t cnt = 0;
						for (; !itImg.IsAtEnd() && mags.size() < maxSamples;
						     ++itImg, ++itGrad, ++cnt)
						{
							if (static_cast<double>(itImg.Get()) > fgThreshold)
								mags.push_back(static_cast<double>(itGrad.Get()));
						}
					}

					double etaF;
					bool atFloor = false;
					if (mags.empty()) { etaF = kMinEta; atFloor = true; }
					else
					{
						std::sort(mags.begin(), mags.end());
						etaF = mags[static_cast<size_t>(kPercentile * mags.size())];
						if (etaF < kMinEta) { etaF = kMinEta; atFloor = true; }
					}
					this->SetFixedEta(etaF);
					std::cout << "  Fixed  η: "
					          << std::scientific << std::setprecision(6) << etaF
					          << std::defaultfloat
					          << "   (fg voxels: " << mags.size() << "/" << maxSamples
					          << ", range [" << minVal << ", " << maxVal
					          << "], fg thr " << fgThreshold << ")"
					          << (atFloor ? " [WARNING: floor — degenerate fixed image]" : "")
					          << std::endl;
				}

				// ── Moving image η ─────────────────────────────────────────
				{
					using GradFilterType = itk::GradientMagnitudeImageFilter<MovingImageType, MovingImageType>;
					typename GradFilterType::Pointer gradFilter = GradFilterType::New();
					gradFilter->SetInput(this->m_MovingImage);
					gradFilter->Update();

					const unsigned long totalPix = this->m_MovingImage->GetLargestPossibleRegion().GetNumberOfPixels();
					const size_t maxSamples = std::max(static_cast<size_t>(20000),
					                                   static_cast<size_t>(totalPix * kSampleFraction));

					double minVal =  std::numeric_limits<double>::max();
					double maxVal = -std::numeric_limits<double>::max();
					{
						itk::ImageRegionConstIterator<MovingImageType> it(
							this->m_MovingImage, this->m_MovingImage->GetLargestPossibleRegion());
						size_t cnt = 0;
						for (; !it.IsAtEnd() && cnt < maxSamples; ++it, ++cnt)
						{
							double v = static_cast<double>(it.Get());
							if (v < minVal) minVal = v;
							if (v > maxVal) maxVal = v;
						}
					}
					const double range = maxVal - minVal;
					const double fgThreshold = (range > 0.0)
					    ? (minVal + kFgIntensityFrac * range) : minVal;

					std::vector<double> mags;
					mags.reserve(maxSamples);
					{
						itk::ImageRegionConstIterator<MovingImageType> itImg(
							this->m_MovingImage, this->m_MovingImage->GetLargestPossibleRegion());
						itk::ImageRegionConstIterator<MovingImageType> itGrad(
							gradFilter->GetOutput(),
							gradFilter->GetOutput()->GetLargestPossibleRegion());
						size_t cnt = 0;
						for (; !itImg.IsAtEnd() && mags.size() < maxSamples;
						     ++itImg, ++itGrad, ++cnt)
						{
							if (static_cast<double>(itImg.Get()) > fgThreshold)
								mags.push_back(static_cast<double>(itGrad.Get()));
						}
					}

					double etaM;
					if (mags.empty()) { etaM = kMinEta; etaMAtFloor = true; }
					else
					{
						std::sort(mags.begin(), mags.end());
						etaM = mags[static_cast<size_t>(kPercentile * mags.size())];
						if (etaM < kMinEta) { etaM = kMinEta; etaMAtFloor = true; }
					}
					this->SetMovingEta(etaM);
					std::cout << "  Moving η: "
					          << std::scientific << std::setprecision(6) << etaM
					          << std::defaultfloat
					          << "   (fg voxels: " << mags.size() << "/" << maxSamples
					          << ", range [" << minVal << ", " << maxVal
					          << "], fg thr " << fgThreshold << ")"
					          << (etaMAtFloor ? " [WARNING: floor — degenerate moving image]" : "")
					          << std::endl;
				}

				// Safety: precomputing the moving NGF with a degenerate η produces
				// huge ratios in the gradient buffer that lead to NaN/Inf and
				// access violations on the first GetDerivative call.  Auto-disable.
				if (etaMAtFloor && this->m_NGFPrecomputeGradient)
				{
					std::cout << "  [SAFETY] Auto-disabling --ngfprecompute because moving η is at floor." << std::endl;
					this->m_NGFPrecomputeGradient = false;
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
		if (this->m_Yota != 0.0 && !this->m_NCDegenerate)
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
		double raw = m_MSE->GetValue(parameters);

		if (this->m_NormalizeMSE)
		{
			// Normalise by the fixed-image intensity range squared so that the
			// contribution sits in [0,1], comparable to MI/NGF/NC values.
			// Range is pre-cached in Initialize() to avoid per-call overhead.
			raw /= this->m_MSEIntensityRangeSquared;
		}

		return static_cast<MeasureType>(raw * this->m_Nu);
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
		double raw = m_GD->GetValue(parameters);
		if (this->m_NormalizeGD)
			raw *= this->m_GDNormalizationFactor;
		return static_cast<MeasureType>(raw * this->m_Rho);
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
		if (this->m_YotaDerivative != 0.0 && !this->m_NCDegenerate)
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
		// NMI derivative: see GetValueAndDerivative — finite-difference NMI
		// derivative is infeasible and unstable for B-spline transforms.
		// Skip with a one-shot warning so NMI contributes value-only.
		if (this->m_SigmaDerivative != 0.0)
		{
			static bool warned = false;
			if (!warned)
			{
				std::cout << "[Mplus] WARNING: NMI derivative requested but skipped "
				             "(finite-difference NMI derivative is infeasible and "
				             "unstable for B-spline transforms; NMI contributes to "
				             "the value only)." << std::endl;
				warned = true;
			}
			g.Fill(0.0);
		}
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
			for (long long p = 0; p < static_cast<long long>(derivative.GetSize()); ++p)
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
					for (long long p = 0; p < static_cast<long long>(derivative.GetSize()); ++p)
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
			for (long long p = 0; p < static_cast<long long>(derivative.GetSize()); ++p)
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
			for (long long p = 0; p < static_cast<long long>(derivative.GetSize()); ++p)
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
	for (long long i = 0; i < static_cast<long long>(der.size()); ++i) {
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
		if (this->m_NormalizeGD)
		{
			const double f = this->m_GDNormalizationFactor;
			for (unsigned int i = 0; i < derivative.size(); ++i)
				derivative[i] *= f;
		}
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
		if (this->m_NormalizeMSE)
		{
			const double inv = 1.0 / this->m_MSEIntensityRangeSquared;
			for (unsigned int i = 0; i < derivative.size(); ++i)
				derivative[i] *= inv;
		}
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
		bool labelValueAlreadyWeighted = false;

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
		// Apply MSE normalisation: divide both value and derivative by range^2
		if (this->m_NormalizeMSE)
		{
			const double inv = 1.0 / this->m_MSEIntensityRangeSquared;
			rawValC *= inv;
			for (unsigned int i = 0; i < nParams; ++i)
				rawDerC[i] *= inv;
		}

		// GD (Gradient Difference)
		// NOTE: GradientDifferenceImageToImageMetric does not override
		// GetValueAndDerivative - the base-class default uses threaded sampling
		// with pre-allocated BSpline weight arrays sized for the metric's own
		// thread count.  Calling GetValueAndDerivative can crash if that count
		// differs from the registration thread count.  Use separate GetValue +
		// GetDerivative calls instead (GD::GetDerivative uses finite differences
		// so it already calls GetValue internally and is always safe).
		if (this->m_Rho != 0.0)
			rawValD = m_GD->GetValue(parameters);
		if (this->m_RhoDerivative != 0.0)
			m_GD->GetDerivative(parameters, rawDerD);
		// Apply GD normalisation: divide both value and derivative by overlap
		// voxel count so GD becomes a mean in [0,1] (comparable to MI/NGF/NC).
		if (this->m_NormalizeGD)
		{
			const double f = this->m_GDNormalizationFactor;
			rawValD *= f;
			for (unsigned int i = 0; i < nParams; ++i)
				rawDerD[i] *= f;
		}

		// NC
		if (!this->m_NCDegenerate)
		{
			if (this->m_Yota != 0.0 && this->m_YotaDerivative != 0.0)
				m_NC->GetValueAndDerivative(parameters, rawValE, rawDerE);
			else if (this->m_Yota != 0.0)
				rawValE = m_NC->GetValue(parameters);
			else if (this->m_YotaDerivative != 0.0)
				m_NC->GetDerivative(parameters, rawDerE);
		}

		// Kappa (label metric)
		if ((this->m_LabelKappa != 0.0 || this->m_LabelKappaDerivative != 0.0)
		    && m_FixedLabelMap && m_MovingLabelMap) {
			if (this->m_LabelKappa != 0.0 && this->m_LabelKappaDerivative != 0.0)
				this->GetKappaValueAndDerivative(parameters, rawValF, rawDerF);
			else if (this->m_LabelKappa != 0.0) {
				rawValF = this->GetKappaValue(parameters);
				labelValueAlreadyWeighted = true;
			}
			else
				this->GetKappaDerivative(parameters, rawDerF);
		}

		// NMI (Normalized Mutual Information).
		//
		// HistogramImageToImageMetric::GetDerivative uses central finite
		// differences: for each of N parameters it perturbs the SHARED transform,
		// recomputes the full histogram, and restores the parameters.  For a
		// B-spline transform with thousands of parameters this is both
		// computationally infeasible (hours per evaluation) and crash-prone
		// because the shared-transform parameter perturbations corrupt state
		// observed by the other (multi-threaded) sub-metrics on the next call.
		// We therefore evaluate NMI as a VALUE-ONLY term in B-spline mode and
		// leave its derivative contribution at zero, even when SigmaDerivative
		// is set.  A one-shot warning is emitted so users know.
		if (this->m_Sigma != 0.0)
			rawValG = m_NMI->GetValue(parameters);
		if (this->m_SigmaDerivative != 0.0)
		{
			static bool warned = false;
			if (!warned)
			{
				std::cout << "[Mplus] WARNING: NMI derivative requested but skipped "
				             "(finite-difference NMI derivative is infeasible and "
				             "unstable for B-spline transforms; NMI contributes to "
				             "the value only)." << std::endl;
				warned = true;
			}
			rawDerG.Fill(0.0);
		}

		// ── Cache weighted per-sub-metric contributions ───────────────────
		this->m_LastValMI    = this->m_Alpha  * rawValA;
		this->m_LastValNGF   = this->m_Lambda * rawValB;
		this->m_LastValMSE   = this->m_Nu     * rawValC;
		this->m_LastValGD    = this->m_Rho    * rawValD;
		this->m_LastValNC    = this->m_Yota   * rawValE;
		this->m_LastValLabel = labelValueAlreadyWeighted ? rawValF
		                                                : this->m_LabelKappa * rawValF;
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
			for (long long p = 0; p < static_cast<long long>(nParams); ++p)
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
					for (long long p = 0; p < static_cast<long long>(nParams); ++p)
						Derivative[p] *= scale;
				}
			}

				Value = this->m_Alpha * rawValA + this->m_Lambda * rawValB
				      + this->m_Nu * rawValC + this->m_Rho * rawValD
				      + this->m_Yota * rawValE
				      + (labelValueAlreadyWeighted ? rawValF : this->m_LabelKappa * rawValF)
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
			for (long long p = 0; p < static_cast<long long>(nParams); ++p)
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
				      + this->m_ScaleLabel
				        * (labelValueAlreadyWeighted ? rawValF : this->m_LabelKappa * rawValF)
				      + this->m_ScaleNMI * this->m_Sigma  * rawValG;
		}
		else
		{
			// Mode 0: consistent weighted sum
			Derivative.SetSize(nParams);
#pragma omp parallel for
			for (long long p = 0; p < static_cast<long long>(nParams); ++p)
				Derivative[p] = this->m_AlphaDerivative  * rawDerA[p]
				              + this->m_LambdaDerivative * rawDerB[p]
				              + this->m_NuDerivative     * rawDerC[p]
				              + this->m_RhoDerivative    * rawDerD[p]
				              + this->m_YotaDerivative   * rawDerE[p]
				              + this->m_LabelKappaDerivative * rawDerF[p]
				              + this->m_SigmaDerivative  * rawDerG[p];

				Value = this->m_Alpha * rawValA + this->m_Lambda * rawValB
				      + this->m_Nu * rawValC + this->m_Rho * rawValD
				      + this->m_Yota * rawValE
				      + (labelValueAlreadyWeighted ? rawValF : this->m_LabelKappa * rawValF)
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

	    // Guard against an all-zero binary mask. SignedMaurerDistanceMapImageFilter
	    // fills the output with NumericTraits::max() (~3.4e+38) when no foreground
	    // is present, which then squares into ~1e+76 in the metric and freezes the
	    // optimizer. Return a finite, large-but-bounded distance map instead.
	    bool hasForeground = false;
	    {
	        itk::ImageRegionConstIterator<TFixedImage> it(
	            binaryImage, binaryImage->GetLargestPossibleRegion());
	        for (; !it.IsAtEnd(); ++it)
	            if (it.Get() != 0.0f) { hasForeground = true; break; }
	    }
	    if (!hasForeground)
	    {
	        std::cout << "[LabelMetric] WARNING: binary mask is empty after "
	                     "resampling onto the reference geometry — returning a "
	                     "constant distance map (label term contributes nothing "
	                     "for this label). Check that the reference image FOV "
	                     "covers the label."
	                  << std::endl;
	        typename TFixedImage::Pointer out = TFixedImage::New();
	        out->CopyInformation(binaryImage);
	        out->SetRegions(binaryImage->GetLargestPossibleRegion());
	        out->Allocate();
	        // Use the user-configured clamp distance so the residual normalises
	        // to ±1 and produces zero gradient (constant field).
	        const double dmax = std::max(1e-6, m_LabelDistanceMax);
	        out->FillBuffer(static_cast<typename TFixedImage::PixelType>(dmax));
	        return out;
	    }

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
	bool
	Mplus<TFixedImage, TMovingImage>::LabelSampleInBand(
	    double dFixed, double dMoving) const
	{
	    if (!m_LabelUseNarrowBand) return true;
	    const double w = std::max(0.0, m_LabelNarrowBandWidth);
	    return (std::abs(dFixed) <= w) || (std::abs(dMoving) <= w);
	}

	template <class TFixedImage, class TMovingImage>
	double
	Mplus<TFixedImage, TMovingImage>::LabelClampDistance(double d) const
	{
	    const double dmax = std::max(1e-6, m_LabelDistanceMax);
	    return std::max(-dmax, std::min(dmax, d));
	}

	template <class TFixedImage, class TMovingImage>
	double
	Mplus<TFixedImage, TMovingImage>::LabelLoss(double residualNorm) const
	{
	    if (!m_LabelUseHuber) return residualNorm * residualNorm;

	    const double delta = std::max(1e-8, m_LabelHuberDelta);
	    const double a = std::abs(residualNorm);
	    if (a <= delta) return residualNorm * residualNorm;
	    return 2.0 * delta * a - delta * delta;
	}

	template <class TFixedImage, class TMovingImage>
	double
	Mplus<TFixedImage, TMovingImage>::LabelLossDerivative(double residualNorm) const
	{
	    if (!m_LabelUseHuber) return 2.0 * residualNorm;

	    const double delta = std::max(1e-8, m_LabelHuberDelta);
	    const double a = std::abs(residualNorm);
	    if (a <= delta) return 2.0 * residualNorm;
	    return (residualNorm >= 0.0) ? 2.0 * delta : -2.0 * delta;
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

	    // Collect unique non-zero labels from each map and count support.  The
	    // metric only evaluates labels present in both maps; one-sided labels
	    // produce invalid/sentinel distance maps and overwhelm the loss.
	    std::map<LabelPixelType, unsigned long long> fixedCounts;
	    std::map<LabelPixelType, unsigned long long> movingCounts;
	    std::set<LabelPixelType> fixedLabels;
	    std::set<LabelPixelType> movingLabels;
	    {
	        itk::ImageRegionConstIterator<LabelImageType> it(
	            m_FixedLabelMap, m_FixedLabelMap->GetLargestPossibleRegion());
	        for (; !it.IsAtEnd(); ++it)
	        {
	            const LabelPixelType v = it.Get();
	            if (v == 0) continue;
	            fixedLabels.insert(v);
	            ++fixedCounts[v];
	        }
	    }
	    {
	        itk::ImageRegionConstIterator<LabelImageType> it(
	            m_MovingLabelMap, m_MovingLabelMap->GetLargestPossibleRegion());
	        for (; !it.IsAtEnd(); ++it)
	        {
	            const LabelPixelType v = it.Get();
	            if (v == 0) continue;
	            movingLabels.insert(v);
	            ++movingCounts[v];
	        }
	    }

	    std::set<LabelPixelType> labelSet;
	    std::set_intersection(
	        fixedLabels.begin(),  fixedLabels.end(),
	        movingLabels.begin(), movingLabels.end(),
	        std::inserter(labelSet, labelSet.begin()));

	    // Report any labels we are dropping so cropping issues are visible.
	    auto reportMissing = [](const std::set<LabelPixelType> & a,
	                            const std::set<LabelPixelType> & b,
	                            const char * sideMissing)
	    {
	        for (auto v : a)
	            if (b.find(v) == b.end())
	                std::cout << "[LabelMetric] WARNING: label "
	                          << static_cast<int>(v)
	                          << " present in " << sideMissing
	                          << " label map only — skipped (no counterpart)."
	                          << std::endl;
	    };
	    reportMissing(fixedLabels,  movingLabels, "fixed");
	    reportMissing(movingLabels, fixedLabels,  "moving");

	    if (labelSet.empty())
	    {
	        std::cout << "[LabelMetric] WARNING: no labels are common to fixed "
	                     "and moving label maps — label term disabled."
	                  << std::endl;
	        return;
	    }

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
	              << labelSet.size() << " label(s)..." << std::endl;

	    const bool needLabelDerivatives = (m_LabelKappaDerivative != 0.0);

	    for (auto L : labelSet)
	    {
	        const unsigned long long nFixed =
	            (fixedCounts.count(L) > 0) ? fixedCounts[L] : 0ull;
	        const unsigned long long nMoving =
	            (movingCounts.count(L) > 0) ? movingCounts[L] : 0ull;
	        if (nFixed == 0ull || nMoving == 0ull)
	        {
	            std::cout << "[LabelMetric]   Label " << static_cast<int>(L)
	                      << " skipped (missing in "
	                      << ((nFixed == 0ull) ? "fixed" : "moving")
	                      << " map)." << std::endl;
	            continue;
	        }

	        // Fixed binary + distance map (in fixed image space)
	        auto fixedBin  = this->BuildBinaryFromLabel(
	            m_FixedLabelMap,  fxSize, fxSpc, fxOrg, fxDir, L);
	        m_FixedDistMaps[L] = this->ComputeSignedDist(fixedBin);

	        // Moving binary + distance map (in moving image space)
	        auto movingBin = this->BuildBinaryFromLabel(
	            m_MovingLabelMap, mvSize, mvSpc, mvOrg, mvDir, L);
	        m_MovingDistMaps[L]     = this->ComputeSignedDist(movingBin);

	        // Interpolators for moving distance map and its gradient
	        auto di = DistInterpType::New();
	        di->SetInputImage(m_MovingDistMaps[L]);
	        m_MovingDistInterps[L] = di;

	        if (needLabelDerivatives)
	        {
	            m_MovingDistGradMaps[L] = this->ComputeGradient(m_MovingDistMaps[L]);
	            auto gi = GradInterpType::New();
	            gi->SetInputImage(m_MovingDistGradMaps[L]);
	            m_MovingDistGradInterps[L] = gi;
	        }

	        m_LabelValues.push_back(L);
	        std::cout << "[LabelMetric]   Label " << static_cast<int>(L) << " done." << std::endl;
	    }

	    if (m_LabelValues.empty())
	        std::cout << "[LabelMetric]   No shared non-zero labels found; label metric disabled."
	                  << std::endl;
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
	    using SizeValueType = itk::SizeValueType;
	    const SizeValueType totalPix = this->m_FixedImage
	        ->GetLargestPossibleRegion().GetNumberOfPixels();
	    const SizeValueType sampleTarget = std::max<SizeValueType>(
	        1, static_cast<SizeValueType>(m_LabelNumberOfSamples));
	    const double distNorm = std::max(1e-6, m_LabelDistanceMax);
	    SizeValueType stride = 1;
	    if (sampleTarget < totalPix)
	    {
	        stride = totalPix / sampleTarget;
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
	        unsigned int n = 0; SizeValueType pixIdx = 0;

	        itk::ImageRegionConstIteratorWithIndex<TFixedImage> it(
	            fixedDistMap, fixedDistMap->GetLargestPossibleRegion());

	        for (it.GoToBegin(); !it.IsAtEnd(); ++it, ++pixIdx)
	        {
	            if (pixIdx % stride != 0) continue;

	            typename TFixedImage::PointType fixedPt;
	            fixedDistMap->TransformIndexToPhysicalPoint(it.GetIndex(), fixedPt);
	            const auto movingPt = this->m_Transform->TransformPoint(fixedPt);

	            if (!movingDistInterp->IsInsideBuffer(movingPt)) continue;

	            const double dFixedRaw  = static_cast<double>(it.Get());
	            const double dMovingRaw = movingDistInterp->Evaluate(movingPt);
	            if (!std::isfinite(dFixedRaw) || !std::isfinite(dMovingRaw)) continue;
	            if (!this->LabelSampleInBand(dFixedRaw, dMovingRaw)) continue;

	            const double dFixed = this->LabelClampDistance(dFixedRaw);
	            const double dMoving = this->LabelClampDistance(dMovingRaw);
	            const double residualNorm = (dFixed - dMoving) / distNorm;
	            if (!std::isfinite(residualNorm)) continue;
	            const double loss = this->LabelLoss(residualNorm);
	            if (!std::isfinite(loss)) continue;

	            sumSqDiff += loss;
	            ++n;

	            // Also accumulate binary Dice stats (dFixed<0 → inside label)
	            bool inF = (dFixedRaw  <= 0.0);
	            bool inM = (dMovingRaw <= 0.0);
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

	    using SizeValueType = itk::SizeValueType;
	    const SizeValueType totalPix = this->m_FixedImage
	        ->GetLargestPossibleRegion().GetNumberOfPixels();
	    const SizeValueType sampleTarget = std::max<SizeValueType>(
	        1, static_cast<SizeValueType>(m_LabelNumberOfSamples));
	    const double distNorm = std::max(1e-6, m_LabelDistanceMax);
	    const double dmax = std::max(1e-6, m_LabelDistanceMax);
	    SizeValueType stride = 1;
	    if (sampleTarget < totalPix)
	    {
	        stride = totalPix / sampleTarget;
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
	        unsigned int n = 0; SizeValueType pixIdx = 0;

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

	            const double dFixedRaw  = static_cast<double>(it.Get());
	            const double dMovingRaw = movingDistInterp->Evaluate(movingPt);
	            if (!std::isfinite(dFixedRaw) || !std::isfinite(dMovingRaw)) continue;
	            if (!this->LabelSampleInBand(dFixedRaw, dMovingRaw)) continue;

	            const double dFixed = this->LabelClampDistance(dFixedRaw);
	            const double dMoving = this->LabelClampDistance(dMovingRaw);
	            const double residualNorm = (dFixed - dMoving) / distNorm;
	            if (!std::isfinite(residualNorm)) continue;
	            const double dLossdResidualNorm = this->LabelLossDerivative(residualNorm);
	            if (!std::isfinite(dLossdResidualNorm)) continue;
	            const double clampMovingFactor = (std::abs(dMovingRaw) <= dmax) ? 1.0 : 0.0;
	            if (clampMovingFactor == 0.0) continue;

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
	                const double contrib = -(dLossdResidualNorm / distNorm)
	                                     * clampMovingFactor * dot;
	                if (std::isfinite(contrib))
	                    localDeriv[j] += contrib;
	            }
	            ++n;
	        }

	        if (n == 0) continue;

        // Per-label scale only. The global label weight (m_LabelKappaDerivative)
        // is applied by the caller (GetDerivative / GetValueAndDerivative),
        // mirroring how every other sub-metric exposes its raw derivative.
        const double scale = kappaL / static_cast<double>(n);
        #pragma omp parallel for
        for (long long j = 0; j < static_cast<long long>(nParams); ++j)
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

	    using SizeValueType = itk::SizeValueType;
	    const SizeValueType totalPix = this->m_FixedImage
	        ->GetLargestPossibleRegion().GetNumberOfPixels();
	    const SizeValueType sampleTarget = std::max<SizeValueType>(
	        1, static_cast<SizeValueType>(m_LabelNumberOfSamples));
	    const double distNorm = std::max(1e-6, m_LabelDistanceMax);
	    const double dmax = std::max(1e-6, m_LabelDistanceMax);
	    SizeValueType stride = 1;
	    if (sampleTarget < totalPix) {
	        stride = totalPix / sampleTarget;
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
	        unsigned int n = 0; SizeValueType pixIdx = 0;

	        itk::ImageRegionConstIteratorWithIndex<TFixedImage> it(
	            fixedDistMap, fixedDistMap->GetLargestPossibleRegion());

	        for (it.GoToBegin(); !it.IsAtEnd(); ++it, ++pixIdx)
	        {
	            if (pixIdx % stride != 0) continue;

	            typename TFixedImage::PointType fixedPt;
	            fixedDistMap->TransformIndexToPhysicalPoint(it.GetIndex(), fixedPt);
	            const auto movingPt = this->m_Transform->TransformPoint(fixedPt);

	            if (!movingDistInterp->IsInsideBuffer(movingPt)) continue;

	            const double dFixedRaw  = static_cast<double>(it.Get());
	            const double dMovingRaw = movingDistInterp->Evaluate(movingPt);
	            if (!std::isfinite(dFixedRaw) || !std::isfinite(dMovingRaw)) continue;
	            if (!this->LabelSampleInBand(dFixedRaw, dMovingRaw)) continue;

	            const double dFixed = this->LabelClampDistance(dFixedRaw);
	            const double dMoving = this->LabelClampDistance(dMovingRaw);
	            const double residualNorm = (dFixed - dMoving) / distNorm;
	            if (!std::isfinite(residualNorm)) continue;

	            // Value accumulation
	            const double loss = this->LabelLoss(residualNorm);
	            if (!std::isfinite(loss)) continue;
	            sumSqDiff += loss;

	            // Derivative accumulation
	            if (kappaD != 0.0 && gradInterp->IsInsideBuffer(movingPt)) {
	                const double dLossdResidualNorm = this->LabelLossDerivative(residualNorm);
	                if (!std::isfinite(dLossdResidualNorm)) continue;
	                const double clampMovingFactor = (std::abs(dMovingRaw) <= dmax) ? 1.0 : 0.0;
	                if (clampMovingFactor == 0.0) continue;

	                const auto gradVec = gradInterp->Evaluate(movingPt);
	                TransformJacobianType jac(TFixedImage::ImageDimension, nParams);
	                this->m_Transform->ComputeJacobianWithRespectToParameters(fixedPt, jac);
	                for (unsigned int j = 0; j < nParams; ++j) {
	                    double dot = 0.0;
	                    for (unsigned int d = 0; d < TFixedImage::ImageDimension; ++d)
	                        dot += static_cast<double>(gradVec[d]) * jac(d, j);
	                    const double contrib = -(dLossdResidualNorm / distNorm)
	                                         * clampMovingFactor * dot;
	                    if (std::isfinite(contrib))
	                        localDeriv[j] += contrib;
	                }
	            }

	            ++n;

	            bool inF = (dFixedRaw  <= 0.0);
	            bool inM = (dMovingRaw <= 0.0);
	            if (inF) ++countF;
	            if (inM) ++countM;
	            if (inF && inM) ++countIntersect;
	        }

	        if (n == 0) continue;

	        totalValue += kappaV * sumSqDiff / static_cast<double>(n);

        // Per-label scale only. Global m_LabelKappaDerivative is applied by
        // the caller (GetValueAndDerivative), matching every other sub-metric.
        const double scaleD = kappaD / static_cast<double>(n);
        #pragma omp parallel for
        for (long long j = 0; j < static_cast<long long>(nParams); ++j)
            derivative[j] += scaleD * localDeriv[j];

        double dice = 0.0;
        if (countF + countM > 0)
            dice = 2.0 * countIntersect / static_cast<double>(countF + countM);
        m_LastDice[L] = dice;
    }

	    // Return unweighted value; caller applies the global label weight.
	    value = static_cast<MeasureType>(totalValue);
	}

} // end namespace itk

#endif

