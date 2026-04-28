/*=========================================================================
Eta is defined as the Habe rdefinition of NGF, different by the itk implemntation (downloaded and add).
 *=========================================================================*/
#ifndef __itkMplus_h
#define __itkMplus_h

#include "itkImageToImageMetric.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "../NGF/NGFImageMetric/NGFImageToImageMetric/Code/itkNormalizedGradientFieldImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkMutualInformationHistogramImageToImageMetric.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkGradientDifferenceImageToImageMetric.h"
#include "itkNormalizedMutualInformationHistogramImageToImageMetric.h"
#include <map>
#include <vector>

namespace itk
{
template <class TFixedImage, class TMovingImage>
class ITK_EXPORT Mplus:
public ImageToImageMetric<TFixedImage, TMovingImage>
{
public:

	/** Standard class typedefs. */
	typedef Mplus     Self;
	typedef ImageToImageMetric<TFixedImage, TMovingImage> Superclass;
	typedef SmartPointer<Self>                            Pointer;
	typedef SmartPointer<const Self>                      ConstPointer;



	/** Method for creation through the object factory. */
	itkNewMacro(Self);

	/** Run-time type information (and related methods). */
	itkTypeMacro(Mplus,
			ImageToImageMetric);

	itkGetMacro( Lambda, double);
	itkSetMacro( Lambda, double);
	
	itkGetMacro( LambdaDerivative, double);
	itkSetMacro( LambdaDerivative, double);

	itkGetMacro( BinNumbers, int);
	itkSetMacro( BinNumbers, int);

	itkGetMacro( MANumberOfSamples, unsigned int);
	itkSetMacro( MANumberOfSamples, unsigned int);
	
	itkGetMacro( MSENumberOfSamples, unsigned int);
	itkSetMacro( MSENumberOfSamples, unsigned int);



	itkGetMacro( UseCachingOfBSplineWeights, bool);
	itkSetMacro( UseCachingOfBSplineWeights, bool);

	itkGetMacro( UseExplicitPDFDerivatives, bool);
	itkSetMacro( UseExplicitPDFDerivatives, bool);

	/** Derivative merge mode:
	 *   0 = consistent weighted sum (safe for LBFGS-B, default)
	 *   1 = normalize + rescale   (RSGD only)
	 *   2 = main-metric adaptive scaling (LBFGS-B safe) */
	itkGetMacro( DerivativeMode, int);
	itkSetMacro( DerivativeMode, int);

	/** Backward-compatible wrapper – true → mode 1, false → mode 0 */
	void SetNormalizeDerivatives(bool v) { m_DerivativeMode = v ? 1 : 0; }
	bool GetNormalizeDerivatives() const { return m_DerivativeMode == 1; }

	/** Index of the "main" metric for mode 2:
	 *   0 = MI (default), 1 = NGF, 2 = MSE, 3 = NC, 4 = Label, 5 = GD, 6 = NMI */
	itkGetMacro( MainMetricIndex, int);
	itkSetMacro( MainMetricIndex, int);	
	
	itkGetMacro( NGFNumberOfSamples, unsigned int);
	itkSetMacro( NGFNumberOfSamples, unsigned int);

	itkSetMacro( NCNumberOfSamples, unsigned int);
	itkGetMacro( NCNumberOfSamples, unsigned int);

	itkSetMacro( GDNumberOfSamples, unsigned int);
	itkGetMacro( GDNumberOfSamples, unsigned int);

	itkSetMacro( NMINumberOfSamples, unsigned int);
	itkGetMacro( NMINumberOfSamples, unsigned int);

	itkSetMacro( NMIBinNumbers, int);
	itkGetMacro( NMIBinNumbers, int);

	// ── Cached per-sub-metric values from last evaluation (weighted contributions) ─
	itkGetMacro( LastValMI,    double);
	itkGetMacro( LastValNGF,   double);
	itkGetMacro( LastValMSE,   double);
	itkGetMacro( LastValNC,    double);
	itkGetMacro( LastValGD,    double);
	itkGetMacro( LastValNMI,   double);
	itkGetMacro( LastValLabel, double);
	itkGetMacro( LastValTotal, double);


	itkGetMacro( FixedEta, double);
	itkSetMacro( FixedEta, double);

	itkGetMacro( MovingEta, double);
	itkSetMacro( MovingEta, double);


	itkGetMacro( AutoEstimateEta, bool );
	itkSetMacro( AutoEstimateEta, bool );

	itkGetMacro( NumberOfThreads, unsigned int);
	itkSetMacro( NumberOfThreads, unsigned int);
	
	
	itkGetMacro( Evaluator, char);
	itkSetMacro( Evaluator, char);
	
	itkGetMacro( Alpha, double);
	itkSetMacro( Alpha, double);
	
	itkGetMacro( AlphaDerivative, double);
	itkSetMacro( AlphaDerivative, double);
	
	itkGetMacro( Nu, double);
	itkSetMacro( Nu, double);
	
	itkGetMacro( NuDerivative, double);
	itkSetMacro( NuDerivative, double);

	/** When true, the raw MSE value is divided by (intensity_range)^2 before
	 *  applying the Nu weight.  This normalises MSE to [0,1] so it stays
	 *  comparable in magnitude to MI / NGF / NC.  Intensity range is taken
	 *  from the fixed image min/max.  Default: false. */
	itkGetMacro( NormalizeMSE, bool);
	itkSetMacro( NormalizeMSE, bool);
	itkBooleanMacro( NormalizeMSE);

	itkGetMacro( Yota, double);
	itkSetMacro( Yota, double);

	itkGetMacro( YotaDerivative, double);
	itkSetMacro( YotaDerivative, double);

	itkGetMacro( ComputeOverlap, bool );
	itkSetMacro( ComputeOverlap, bool );

	itkGetMacro( OverlapPadding, unsigned int);
	itkSetMacro( OverlapPadding, unsigned int);

	/** Forward an intensity threshold to all active sub-metrics during Initialize(). */
	void SetFixedImageThreshold(double t) { m_FixedImageThreshold = t; m_UseFixedImageThreshold = true; }
	double GetFixedImageThreshold() const { return m_FixedImageThreshold; }
	bool GetUseFixedImageThreshold() const { return m_UseFixedImageThreshold; }

	itkGetMacro( Rho, double);
	itkSetMacro( Rho, double);

	itkGetMacro( RhoDerivative, double);
	itkSetMacro( RhoDerivative, double);

	itkGetMacro( Sigma, double);
	itkSetMacro( Sigma, double);

	itkGetMacro( SigmaDerivative, double);
	itkSetMacro( SigmaDerivative, double);

	
	void SetNGFSpacing(const typename TFixedImage::SpacingType& spacing) { m_NGFSpacing = spacing; }
	typename TFixedImage::SpacingType GetNGFSpacing() const { return m_NGFSpacing; }

	/** When true, the NGF of the moving image is precomputed once and
	 *  its vector field is resampled each iteration instead of recomputing
	 *  the gradient from the resampled scalar image.  Faster but approximate
	 *  (ignores transform-Jacobian rotation of gradient vectors). */
	itkGetMacro( NGFPrecomputeGradient, bool);
	itkSetMacro( NGFPrecomputeGradient, bool);
	itkBooleanMacro( NGFPrecomputeGradient);

	// ── Label-map / ROI metric (kappa term) ──────────────────────────────────
	/** Integer label pixel type.  Short accommodates up to 32767 structures. */
	typedef short                                                        LabelPixelType;
	typedef itk::Image<LabelPixelType, TFixedImage::ImageDimension>      LabelImageType;
	typedef typename LabelImageType::ConstPointer                        LabelImageConstPointer;

	/** Gradient image type used internally for distance-map gradients. */
	typedef itk::CovariantVector<float, TFixedImage::ImageDimension>     GradientPixelType;
	typedef itk::Image<GradientPixelType, TFixedImage::ImageDimension>   GradientImageType;
	typedef typename GradientImageType::Pointer                          GradientImagePointer;

	/** Per-label kappa weight map: key = label integer value, value = weight. */
	typedef std::map<LabelPixelType, double>                             LabelWeightMapType;

	void SetFixedLabelMap(LabelImageConstPointer img)  { m_FixedLabelMap  = img; }
	void SetMovingLabelMap(LabelImageConstPointer img) { m_MovingLabelMap = img; }
	LabelImageConstPointer GetFixedLabelMap()  const   { return m_FixedLabelMap;  }
	LabelImageConstPointer GetMovingLabelMap() const   { return m_MovingLabelMap; }

	/** Global kappa weight – the distance-map MSE is multiplied by this.
	 *  Set to 0 (default) to disable the label metric entirely. */
	itkSetMacro(LabelKappa, double);
	itkGetMacro(LabelKappa, double);
	itkSetMacro(LabelKappaDerivative, double);
	itkGetMacro(LabelKappaDerivative, double);

	/** Per-label overrides.  Each label value maps to a weight ∈ [0,1].
	 *  Labels absent from the map receive a flat weight of 1/nLabels. */
	void SetLabelKappaWeights(const LabelWeightMapType & w)          { m_LabelKappaWeights      = w; }
	void SetLabelKappaDerivativeWeights(const LabelWeightMapType & w){ m_LabelKappaDerivWeights  = w; }
	const LabelWeightMapType & GetLabelKappaWeights()          const  { return m_LabelKappaWeights; }
	const LabelWeightMapType & GetLabelKappaDerivativeWeights() const { return m_LabelKappaDerivWeights; }

		itkSetMacro(LabelNumberOfSamples, unsigned int);
		itkGetMacro(LabelNumberOfSamples, unsigned int);

		/** Clamp signed distances to [-LabelDistanceMax, LabelDistanceMax] before
		 *  evaluating the label loss. Units are millimetres when spacing is in mm. */
		itkSetMacro(LabelDistanceMax, double);
		itkGetMacro(LabelDistanceMax, double);

		/** Optional narrow-band mode: only voxels close to either boundary contribute. */
		itkSetMacro(LabelUseNarrowBand, bool);
		itkGetMacro(LabelUseNarrowBand, bool);
		itkBooleanMacro(LabelUseNarrowBand);
		itkSetMacro(LabelNarrowBandWidth, double);
		itkGetMacro(LabelNarrowBandWidth, double);

		/** Optional robust Huber loss on normalized residuals. */
		itkSetMacro(LabelUseHuber, bool);
		itkGetMacro(LabelUseHuber, bool);
		itkBooleanMacro(LabelUseHuber);
		itkSetMacro(LabelHuberDelta, double);
		itkGetMacro(LabelHuberDelta, double);

	/** Read-only access to the per-label Dice coefficients stored after
	 *  the most recent GetKappaValue() call. Key = label value, value ∈ [0,1]. */
	const std::map<LabelPixelType, double> & GetLastDice() const { return m_LastDice; }








	/** Types inherited from Superclass. */
	typedef typename Superclass::TransformType                TransformType;
	typedef typename Superclass::TransformPointer             TransformPointer;
	typedef typename Superclass::TransformJacobianType        TransformJacobianType;
	typedef typename Superclass::InterpolatorType             InterpolatorType;
	typedef typename Superclass::MeasureType                  MeasureType;
	typedef typename Superclass::DerivativeType               DerivativeType;
	typedef typename Superclass::ParametersType               ParametersType;
	typedef typename Superclass::FixedImageType               FixedImageType;
	typedef typename Superclass::MovingImageType              MovingImageType;
	typedef typename Superclass::MovingImagePointType         MovingImagePointType;
	typedef typename Superclass::FixedImageConstPointer       FixedImageConstPointer;
	typedef typename Superclass::MovingImageConstPointer      MovingImageConstPointer;
	typedef typename Superclass::CoordinateRepresentationType CoordinateRepresentationType;
	typedef typename Superclass::FixedImageSampleContainer    FixedImageSampleContainer;
	typedef typename Superclass::ImageDerivativesType         ImageDerivativesType;
	typedef typename Superclass::WeightsValueType             WeightsValueType;
	typedef typename Superclass::IndexValueType               IndexValueType;



	// Needed for evaluation of Jacobian.
	typedef typename Superclass::FixedImagePointType FixedImagePointType;

	/** The moving image dimension. */
	itkStaticConstMacro(MovingImageDimension, unsigned int,
			MovingImageType::ImageDimension);

	typedef typename    TFixedImage::RegionType      RegionType;
	typedef typename    TFixedImage::SizeType        SizeType;
	typedef typename    TFixedImage::IndexType       IndexType;


	virtual void Initialize(void);
	void print() const;
	MeasureType GetValue(const ParametersType & parameters) const;
	MeasureType GetNGFValue(const ParametersType & parameters) const;
	MeasureType GetMAValue(const ParametersType & parameters) const;
	MeasureType GetMSEValue(const ParametersType & parameters) const;
	MeasureType GetNCValue(const ParametersType & parameters) const;
	MeasureType GetGDValue(const ParametersType & parameters) const;
	MeasureType GetNMIValue(const ParametersType & parameters) const;
	MeasureType GetKappaValue(const ParametersType & parameters) const;
	void        GetKappaDerivative(const ParametersType & parameters,
	                               DerivativeType & derivative) const;
	void        GetKappaValueAndDerivative(const ParametersType & parameters,
	                               MeasureType & value,
	                               DerivativeType & derivative) const;

	/** Get the derivatives of the match measure. */
	void GetDerivative(const ParametersType & parameters,
			DerivativeType & Derivative) const;
	void GetMADerivative(const ParametersType & parameters,
			DerivativeType & Derivative) const;
	void GetNGFDerivative(const ParametersType & parameters,
			DerivativeType & Derivative) const;
	void GetMSEDerivative(const ParametersType & parameters,
			DerivativeType & Derivative) const;
	void GetNCDerivative(const ParametersType & parameters, DerivativeType & Derivative) const;
	void GetGDDerivative(const ParametersType & parameters, DerivativeType & Derivative) const;
	void GetNMIDerivative(const ParametersType & parameters, DerivativeType & Derivative) const;
			


	void NormalizeDerivative(DerivativeType & Derivative) const;
	/**  Get the value and derivatives for single valued optimizers. */
	void GetValueAndDerivative(const ParametersType & parameters,MeasureType & Value,DerivativeType & Derivative) const;

	//void SetRegularizationTerm(double s);
	void NormalizeComponents(DerivativeType & derivative) const
	{
			double norm = 0.0;
	#pragma omp parallel for reduction(+:norm)
	for (long long i = 0; i < static_cast<long long>(derivative.size()); ++i) {
		norm += derivative[i] * derivative[i];
	}
	norm = std::sqrt(norm);

	// Check if norm is large enough to avoid division by zero
	if (norm > 1.0e-10) {
		#pragma omp parallel for
		for (long long i = 0; i < static_cast<long long>(derivative.size()); ++i) {
			derivative[i] /= norm;
		}
	}
	}

	void RescaleComponents(DerivativeType & derivative)
	{
			{

double minVal = *std::min_element(derivative.begin(), derivative.end());
double maxVal = *std::max_element(derivative.begin(), derivative.end());
// Check if maxVal and minVal are not equal to avoid division by zero
if (maxVal != minVal)
{
	#pragma omp parallel for
	for (long long i = 0; i < static_cast<long long>(derivative.size()); ++i)
	{
		derivative[i] = 2 * (derivative[i] - minVal) / (maxVal - minVal) - 1;
	}
}
		}
	}

protected:
	int m_BinNumbers;
	unsigned int m_MANumberOfSamples;
	unsigned int m_NGFNumberOfSamples;
	unsigned int m_MSENumberOfSamples;
	unsigned int m_NCNumberOfSamples;
	unsigned int m_GDNumberOfSamples;
	unsigned int m_NMINumberOfSamples;
	int m_NMIBinNumbers;
	double m_Rho;
	double m_RhoDerivative;
	double m_Sigma;
	double m_SigmaDerivative;
	double m_FixedEta;
	double m_MovingEta;
	double m_Lambda;
	double m_LambdaDerivative;
	double m_Yota;
	double m_YotaDerivative;

	unsigned int m_NumberOfThreads;
	char m_Evaluator;
	double m_Alpha;
	double m_AlphaDerivative;
	double m_Nu;
	double m_NuDerivative;
	bool m_UseCachingOfBSplineWeights;
	bool m_UseExplicitPDFDerivatives;
	bool m_NormalizeMSE;              // divide raw MSE by intensity-range^2 before Nu weighting
	double m_MSEIntensityRangeSquared; // cached (range^2), computed once in Initialize()
	int  m_DerivativeMode;     // 0=consistent, 1=normalized, 2=main-metric
	int  m_MainMetricIndex;    // 0=MI,1=NGF,2=MSE,3=NC,4=Label,5=GD,6=NMI
	bool   m_AutoEstimateEta;
	bool m_ComputeOverlap;
	bool m_NGFPrecomputeGradient;
	unsigned int m_OverlapPadding;
	double m_FixedImageThreshold;
	bool   m_UseFixedImageThreshold;

	/** Cached per-metric derivative-norm scale factors (set in GetDerivative mode 2,
	 *  consumed in GetValue mode 2 for value/derivative consistency). */
	mutable double m_ScaleMA, m_ScaleNGF, m_ScaleMSE, m_ScaleNC, m_ScaleLabel, m_ScaleGD, m_ScaleNMI;

	/** Cached weighted per-sub-metric contributions from the last GetValue /
	 *  GetValueAndDerivative call. Updated every evaluation — zero-overhead since
	 *  the values are already computed as part of the normal metric evaluation. */
	mutable double m_LastValMI;
	mutable double m_LastValNGF;
	mutable double m_LastValMSE;
	mutable double m_LastValNC;
	mutable double m_LastValGD;
	mutable double m_LastValNMI;
	mutable double m_LastValLabel;
	mutable double m_LastValTotal;

	// ── label metric members ──────────────────────────────────────────────────
	LabelImageConstPointer  m_FixedLabelMap;
	LabelImageConstPointer  m_MovingLabelMap;
		double                  m_LabelKappa;
		double                  m_LabelKappaDerivative;
		LabelWeightMapType      m_LabelKappaWeights;
		LabelWeightMapType      m_LabelKappaDerivWeights;
		unsigned int            m_LabelNumberOfSamples;
		double                  m_LabelDistanceMax;
		bool                    m_LabelUseNarrowBand;
		double                  m_LabelNarrowBandWidth;
		bool                    m_LabelUseHuber;
		double                  m_LabelHuberDelta;

	// Distance maps and interpolators per label (computed once in Initialize)
	typedef std::map<LabelPixelType, typename TFixedImage::Pointer>        DistMapContainer;
	typedef std::map<LabelPixelType, GradientImagePointer>                 GradMapContainer;
	typedef itk::LinearInterpolateImageFunction<TFixedImage, double>       DistInterpType;
	typedef itk::VectorLinearInterpolateImageFunction<
	    GradientImageType, double>                                          GradInterpType;
	typedef std::map<LabelPixelType, typename DistInterpType::Pointer>     DistInterpContainer;
	typedef std::map<LabelPixelType, typename GradInterpType::Pointer>     GradInterpContainer;

	mutable DistMapContainer       m_FixedDistMaps;
	mutable DistMapContainer       m_MovingDistMaps;
	mutable GradMapContainer       m_MovingDistGradMaps;
	mutable DistInterpContainer    m_MovingDistInterps;
	mutable GradInterpContainer    m_MovingDistGradInterps;
	mutable std::vector<LabelPixelType>              m_LabelValues;
	mutable std::map<LabelPixelType, double>         m_LastDice;

protected:

	Mplus();
	virtual ~Mplus();
	void PrintSelf(std::ostream & os, Indent indent) const;
	typedef MattesMutualInformationImageToImageMetric<FixedImageType,MovingImageType> MattesType;
	typedef NormalizedGradientFieldImageToImageMetric<FixedImageType,MovingImageType> NGFType;
	typedef LinearInterpolateImageFunction<FixedImageType,double > LFType;
	typedef MeanSquaresImageToImageMetric<FixedImageType,MovingImageType> MSEType;
	typedef NormalizedCorrelationImageToImageMetric<FixedImageType,MovingImageType>    NCType;
	typedef GradientDifferenceImageToImageMetric<FixedImageType,MovingImageType>       GDType;
	typedef NormalizedMutualInformationHistogramImageToImageMetric<FixedImageType,MovingImageType> NMIType;



	typedef typename NGFType::MovingNGFType MovingNGFType;
	typedef typename NGFType::FixedNGFType  FixedNGFType;

	/** Compute max–min over one derivative vector */
	double ComputeDerivativeRange(const DerivativeType & der) const;

	/** Compute mean of one derivative vector */
	double ComputeDerivativeMean(const DerivativeType & der) const;

	/** Compute standard deviation of one derivative vector */
	double ComputeDerivativeStdDev(const DerivativeType & der) const;

	/** cached range of the *last* sub‐metric derivative */
	mutable double m_RangeDerivatives;

	mutable double m_STDDerivatives;
	mutable double m_MeanDerivatives;

	double  ComputeDerivativeNorm(const DerivativeType & derivative) const;
	mutable double m_LastComponentNorm;



private:
	// purposely not implemented
	Mplus(const Self &);
	// purposely not implemented
	void operator=(const Self &);

	typename MattesType::Pointer m_MA;
	typename NGFType::Pointer m_NGF;
	typename MSEType::Pointer m_MSE;
	typename NCType::Pointer m_NC;
	typename GDType::Pointer m_GD;
	typename NMIType::Pointer m_NMI;

	typename LFType::Pointer m_INTERNALL_interpolator;
	

	typename TFixedImage::SpacingType m_NGFSpacing;

	// ── label metric private helpers ──────────────────────────────────────────
	void InitializeLabelMetric();

	/** Build a float binary image (1.0 where label==L, 0.0 elsewhere) on the
	 *  coordinate grid of @p refSpacing / @p refOrigin / @p refDirection. */
	typename TFixedImage::Pointer
	BuildBinaryFromLabel(
	    LabelImageConstPointer                          labelMap,
	    const typename TFixedImage::SizeType          & refSize,
	    const typename TFixedImage::SpacingType       & refSpacing,
	    const typename TFixedImage::PointType         & refOrigin,
	    const typename TFixedImage::DirectionType     & refDirection,
	    LabelPixelType                                  L) const;

	/** Signed Euclidean distance transform: negative inside label, positive outside. */
	typename TFixedImage::Pointer
	ComputeSignedDist(const typename TFixedImage::Pointer & binaryImage) const;

	/** Gradient of a scalar float image (via GradientRecursiveGaussianImageFilter). */
		GradientImagePointer
		ComputeGradient(const typename TFixedImage::Pointer & image) const;

		bool LabelSampleInBand(double dFixed, double dMoving) const;
		double LabelClampDistance(double d) const;
		double LabelLoss(double residualNorm) const;
		double LabelLossDerivative(double residualNorm) const;




	};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMplus.hxx"
#endif

#endif
