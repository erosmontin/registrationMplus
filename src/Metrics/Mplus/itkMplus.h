/*=========================================================================
Eta is defined as the Habe rdefinition of NGF, different by the itk implemntation (downloaded and add).
 *=========================================================================*/
#ifndef __itkMplus_h
#define __itkMplus_h

#include "itkImageToImageMetric.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "../NGF/NGFImageMetric/NGFImageToImageMetric/Code/itkNormalizedGradientFieldImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkMutualInformationHistogramImageToImageMetric.h"

#include "itkNormalizedMutualInformationHistogramImageToImageMetric.h"
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

	itkGetMacro( NormalizeDerivatives, bool);
	itkSetMacro( NormalizeDerivatives, bool);	
	
	itkGetMacro( NGFNumberOfSamples, unsigned int);
	itkSetMacro( NGFNumberOfSamples, unsigned int);

	itkSetMacro(NMINumberOfSamples, unsigned int);
	itkGetMacro( NMINumberOfSamples, unsigned int);
	itkSetMacro( HMINumberOfSamples, unsigned int);
	itkGetMacro( HMINumberOfSamples, unsigned int);


	itkGetMacro( FixedEta, double);
	itkSetMacro( FixedEta, double);

	itkGetMacro( MovingEta, double);
	itkSetMacro( MovingEta, double);

// <<< new: toggle auto-estimate η for NGF
  itkGetMacro( AutoEstimateEta, bool );
  itkSetMacro( AutoEstimateEta, bool );
// >>>

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

	itkGetMacro( Yota, double);
	itkSetMacro( Yota, double);

	itkGetMacro( YotaDerivative, double);
	itkSetMacro( YotaDerivative, double);


	itkGetMacro( Rho, double);
	itkSetMacro( Rho, double);

	itkGetMacro( RhoDerivative, double);
	itkSetMacro( RhoDerivative, double);

	
	void SetNGFSpacing(const typename TFixedImage::SpacingType& spacing) { m_NGFSpacing = spacing; }
	typename TFixedImage::SpacingType GetNGFSpacing() const { return m_NGFSpacing; }







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
	MeasureType GetHMIValue(const ParametersType & parameters) const;        
	MeasureType GetNMIValue(const ParametersType & parameters) const;

	/** Get the derivatives of the match measure. */
	void GetDerivative(const ParametersType & parameters,
			DerivativeType & Derivative) const;
	void GetMADerivative(const ParametersType & parameters,
			DerivativeType & Derivative) const;
	void GetNGFDerivative(const ParametersType & parameters,
			DerivativeType & Derivative) const;
	void GetMSEDerivative(const ParametersType & parameters,
			DerivativeType & Derivative) const;
	void        GetHMIDerivative(const ParametersType & parameters, 
                              DerivativeType & Derivative) const;       
	void GetNMIDerivative(const ParametersType & parameters,
			DerivativeType & Derivative) const;
	void NormalizeDerivative(DerivativeType & Derivative) const;
	/**  Get the value and derivatives for single valued optimizers. */
	void GetValueAndDerivative(const ParametersType & parameters,MeasureType & Value,DerivativeType & Derivative) const;

	//void SetRegularizationTerm(double s);
	void NormalizeComponents(DerivativeType & derivative)
	{
			double norm = 0.0;
	#pragma omp parallel for reduction(+:norm)
	for (unsigned int i = 0; i < derivative.size(); ++i) {
		norm += derivative[i] * derivative[i];
	}
	norm = std::sqrt(norm);

	// Check if norm is not zero to avoid division by zero
	if (norm != 1.0e-10) {
		#pragma omp parallel for
		for (unsigned int i = 0; i < derivative.size(); ++i) {
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
	for (unsigned int i = 0; i < derivative.size(); ++i)
	{
		derivative[i] = 2 * (derivative[i] - minVal) / (maxVal - minVal) - 1;
	}
}
		}
	}

	int m_BinNumbers;
	unsigned int m_MANumberOfSamples;
	unsigned int m_NGFNumberOfSamples;
	unsigned int m_MSENumberOfSamples;
	unsigned int m_HMINumberOfSamples;
	unsigned int m_NMINumberOfSamples;
	double m_FixedEta;
	double m_MovingEta;
	double m_Lambda;
	double m_LambdaDerivative;
	double m_Yota;
	double m_YotaDerivative;
	double m_Rho;
	double m_RhoDerivative;
	unsigned int m_NumberOfThreads;
	char m_Evaluator;
	double m_Alpha;
	double m_AlphaDerivative;
	double m_Nu;
	double m_NuDerivative;
	bool m_UseCachingOfBSplineWeights;
	bool m_UseExplicitPDFDerivatives;
	bool m_NormalizeDerivatives;
	// <<< new member
	bool   m_AutoEstimateEta;
	// >>>

protected:

	Mplus();
	virtual ~Mplus();
	void PrintSelf(std::ostream & os, Indent indent) const;
	typedef MattesMutualInformationImageToImageMetric<FixedImageType,MovingImageType> MattesType;
	typedef NormalizedGradientFieldImageToImageMetric<FixedImageType,MovingImageType> NGFType;
	typedef LinearInterpolateImageFunction<FixedImageType,double > LFType;
	typedef MeanSquaresImageToImageMetric<FixedImageType,MovingImageType> MSEType;
	typedef MutualInformationHistogramImageToImageMetric<
            FixedImageType,MovingImageType>                HMIType;
	typedef NormalizedMutualInformationHistogramImageToImageMetric<FixedImageType,MovingImageType> NMIType;



	typedef typename NGFType::MovingNGFType MovingNGFType;
	typedef typename NGFType::FixedNGFType  FixedNGFType;




private:
	// purposely not implemented
	Mplus(const Self &);
	// purposely not implemented
	void operator=(const Self &);

	typename MattesType::Pointer m_MA;
	typename NGFType::Pointer m_NGF;
	typename MSEType::Pointer m_MSE;
	typename HMIType::Pointer m_HMI;                                         
	typename NMIType::Pointer m_NMI;
	typename LFType::Pointer m_INTERNALL_interpolator;
	
	

	typename TFixedImage::SpacingType m_NGFSpacing;




};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMplus.hxx"
#endif

#endif
