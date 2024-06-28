#ifndef __itkMANGFMSE_hxx
#define __itkMANGFMSE_hxx

#include "itkMANGFMSE.h"
#include "../NGF/NGFImageMetric/NGFImageToImageMetric/Code/itkNGFMetricKernel.h"

namespace itk
{
	/**
	 * Constructor
	 */
	template <class TFixedImage, class TMovingImage>
	MANGFMSE<TFixedImage, TMovingImage>::MANGFMSE()
	{
		m_MA = MattesType::New();
		m_NGF = NGFType::New();
		m_MSE = MSEType::New();
		m_INTERNALL_interpolator = LFType::New();
		m_Lambda = 1.0;
		m_LambdaDerivative = m_Lambda;
		m_NGFNumberOfSamples = 20000;
		m_MANumberOfSamples = 20000;
		m_MSENumberOfSamples = 20000;
		m_BinNumbers = 50;
		m_NumberOfThreads = 1;
		m_Evaluator = 0;
		m_FixedEta = 1;
		m_MovingEta = 1;
		m_Alpha = 1.0;
		m_AlphaDerivative = 1.0;
		m_Nu = 1.0;
		m_NuDerivative = 1.0;

		m_UseCachingOfBSplineWeights = true;
		m_UseExplicitPDFDerivatives = true;
		m_NormalizeDerivatives = true;
	}
	template <class TFixedImage, class TMovingImage>
	void
	MANGFMSE<TFixedImage, TMovingImage>::NormalizeDerivative(DerivativeType &derivative) const
	{
		if (m_NormalizeDerivatives)
		{

double minVal = *std::min_element(derivative.begin(), derivative.end());
double maxVal = *std::max_element(derivative.begin(), derivative.end());
// Check if maxVal and minVal are not equal to avoid division by zero
if (maxVal != minVal)
{
	#pragma omp parallel for
	for (unsigned int i = 0; i < derivative.size(); ++i)
	{
		derivative[i] = 2 * (derivative[i] - minVal) / (maxVal - minVal) - 1;
	}
}
		}
	}

	template <class TFixedImage, class TMovingImage>
	MANGFMSE<TFixedImage, TMovingImage>::~MANGFMSE() {}

	/**
	 * Initialize
	 */

	template <class TFixedImage, class TMovingImage>
	void
	MANGFMSE<TFixedImage, TMovingImage>::Initialize(void)
	{

		Superclass::Initialize();

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

		m_MSE->SetFixedImage(this->m_FixedImage);
		m_MSE->SetMovingImage(this->m_MovingImage);
		m_MSE->SetInterpolator(this->m_Interpolator); // m_MSE->SetInterpolator(this->m_INTERNALL_interpolator);
		m_MSE->SetTransform(this->m_Transform);
		m_MSE->SetFixedImageRegion(this->m_FixedImage->GetRequestedRegion());
		m_MSE->UseAllPixelsOff();
		m_MSE->SetNumberOfThreads(this->m_NumberOfThreads);
		m_MSE->SetUseCachingOfBSplineWeights(this->m_UseCachingOfBSplineWeights);
		m_MSE->SetNumberOfSpatialSamples(this->m_MSENumberOfSamples);
		m_MSE->ReinitializeSeed();
		m_MSE->Initialize();

		// add a resampling filter to the NGF metric

		m_NGF->SetFixedImage(this->m_FixedImage);
		m_NGF->SetMovingImage(this->m_MovingImage);
		m_NGF->SetInterpolator(this->m_INTERNALL_interpolator);
		m_NGF->SetTransform(this->m_Transform);
		m_NGF->SetNumberOfSpatialSamples(this->m_NGFNumberOfSamples);
		m_NGF->SetNumberOfThreads(this->m_NumberOfThreads);
		m_NGF->SetUseCachingOfBSplineWeights(this->m_UseCachingOfBSplineWeights);
		m_NGF->UseAllPixelsOff();
		RegionType inputRegion = this->m_FixedImage->GetRequestedRegion();
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

	template <class TFixedImage, class TMovingImage>
	typename MANGFMSE<TFixedImage, TMovingImage>::MeasureType
	MANGFMSE<TFixedImage, TMovingImage>::GetValue(const ParametersType &parameters) const
	{
		double a, b, c;

		a = 0.0;
		if (this->m_Alpha != 0.0)
			a = this->GetMAValue(parameters);
		b = 0.0;
		if (this->m_Lambda != 0.0)
			b = this->GetNGFValue(parameters);
		c = 0.0;
		if (this->m_Nu != 0.0)
			c = this->GetMSEValue(parameters);
		return a + b + c;
	}

	template <class TFixedImage, class TMovingImage>
	typename MANGFMSE<TFixedImage, TMovingImage>::MeasureType
	MANGFMSE<TFixedImage, TMovingImage>::GetNGFValue(const ParametersType &parameters) const
	{
		return static_cast<MeasureType>(m_NGF->GetValue(parameters) * this->m_Lambda);
	}

	template <class TFixedImage, class TMovingImage>
	typename MANGFMSE<TFixedImage, TMovingImage>::MeasureType
	MANGFMSE<TFixedImage, TMovingImage>::GetMAValue(const ParametersType &parameters) const
	{
		return static_cast<MeasureType>(m_MA->GetValue(parameters) * this->m_Alpha);
	}

	template <class TFixedImage, class TMovingImage>
	typename MANGFMSE<TFixedImage, TMovingImage>::MeasureType
	MANGFMSE<TFixedImage, TMovingImage>::GetMSEValue(const ParametersType &parameters) const
	{
		return static_cast<MeasureType>(m_MSE->GetValue(parameters) * this->m_Nu);
	}

	template <class TFixedImage, class TMovingImage>
	void
	MANGFMSE<TFixedImage, TMovingImage>::GetDerivative(const ParametersType &parameters, DerivativeType &derivative) const
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

		derivative = a;
double derivativeSum = this->m_AlphaDerivative + this->m_LambdaDerivative + this->m_NuDerivative;
#pragma omp parallel for
		for (long unsigned int p = 0; p < derivative.GetSize(); p++)
		{	
			derivative[p] = a[p] * this->m_AlphaDerivative + this->m_LambdaDerivative * b[p] + this->m_NuDerivative * c[p];
			derivative[p] = derivative[p] / derivativeSum;	
		}
	}

	template <class TFixedImage, class TMovingImage>
	void
	MANGFMSE<TFixedImage, TMovingImage>::GetMADerivative(const ParametersType &parameters, DerivativeType &derivative) const
	{

		m_MA->GetDerivative(parameters, derivative);
		// normalize the derivative
		this->NormalizeDerivative(derivative);


	}

	template <class TFixedImage, class TMovingImage>
	void
	MANGFMSE<TFixedImage, TMovingImage>::GetNGFDerivative(const ParametersType &parameters, DerivativeType &derivative) const
	{

		m_NGF->GetDerivative(parameters, derivative);

		this->NormalizeDerivative(derivative);

	}

	template <class TFixedImage, class TMovingImage>
	void
	MANGFMSE<TFixedImage, TMovingImage>::GetMSEDerivative(const ParametersType &parameters, DerivativeType &derivative) const
	{
		m_MSE->GetDerivative(parameters, derivative);

		this->NormalizeDerivative(derivative);

	}

	template <class TFixedImage, class TMovingImage>
	void
	MANGFMSE<TFixedImage, TMovingImage>::GetValueAndDerivative(const ParametersType &parameters, MeasureType &Value, DerivativeType &Derivative) const
	{
		Value = this->GetValue(parameters);
		this->GetDerivative(parameters, Derivative);
	}

	template <class TImageType, class TMovingImage>
	void
	MANGFMSE<TImageType, TMovingImage>::
		PrintSelf(std::ostream &os, Indent indent) const
	{
		Superclass::PrintSelf(os, indent);
	}

} // end namespace itk

#endif
