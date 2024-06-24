/*=========================================================================
Eta is defined as the Habe rdefinition of NGF, different by the itk implemntation (downloaded and add).
 *=========================================================================*/
#ifndef __itkMANGF_h
#define __itkMANGF_h

#include "itkImageToImageMetric.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "../NGF/NGFImageMetric/NGFImageToImageMetric/Code/itkNormalizedGradientFieldImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"

namespace itk
{
template <class TFixedImage, class TMovingImage>
class ITK_EXPORT MANGF:
public ImageToImageMetric<TFixedImage, TMovingImage>
{
public:

	/** Standard class typedefs. */
	typedef MANGF     Self;
	typedef ImageToImageMetric<TFixedImage, TMovingImage> Superclass;
	typedef SmartPointer<Self>                            Pointer;
	typedef SmartPointer<const Self>                      ConstPointer;



	/** Method for creation through the object factory. */
	itkNewMacro(Self);

	/** Run-time type information (and related methods). */
	itkTypeMacro(MANGF,
			ImageToImageMetric);

	itkGetMacro( lambda, double);
	itkSetMacro( lambda, double);
	
	itkGetMacro( lambdaDerivative, double);
	itkSetMacro( lambdaDerivative, double);

	itkGetMacro( BinNumbers, int);
	itkSetMacro( BinNumbers, int);
	itkGetMacro( MAnumberOfSamples, unsigned int);
	itkSetMacro( MAnumberOfSamples, unsigned int);

	itkGetMacro(lambdaMI, double);
	itkSetMacro(lambdaMI, double);

	itkGetMacro(lambdaDerivativeMI, double);
	itkSetMacro(lambdaDerivativeMI, double);



	itkGetMacro( NGFnumberOfSamples, unsigned int);
	itkSetMacro( NGFnumberOfSamples, unsigned int);

	itkGetMacro( FixedEta, double);
	itkSetMacro( FixedEta, double);

	itkGetMacro( MovingEta, double);
	itkSetMacro( MovingEta, double);

	itkGetMacro( NumberOfThreads, unsigned int);
	itkSetMacro( NumberOfThreads, unsigned int);
	
	
	itkGetMacro( Evaluator, char);
	itkSetMacro( Evaluator, char);
	



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
	typedef typename    TFixedImage::SpacingType     SpacingType;


	virtual void Initialize(void)
	throw ( ExceptionObject );

	MeasureType GetValue(const ParametersType & parameters) const;
	MeasureType GetNGFValue(const ParametersType & parameters) const;
	MeasureType GetMAValue(const ParametersType & parameters) const;

	/** Get the derivatives of the match measure. */
	void GetDerivative(const ParametersType & parameters,
			DerivativeType & Derivative) const;
	void GetMADerivative(const ParametersType & parameters,
			DerivativeType & Derivative) const;
	void GetNGFDerivative(const ParametersType & parameters,
			DerivativeType & Derivative) const;

	/**  Get the value and derivatives for single valued optimizers. */
	void GetValueAndDerivative(const ParametersType & parameters,MeasureType & Value,DerivativeType & Derivative) const;

	//void SetRegularizationTerm(double s);


	int m_BinNumbers;
	unsigned int m_MAnumberOfSamples;
	unsigned int m_NGFnumberOfSamples;
	double m_FixedEta;
	double m_MovingEta;
	double m_lambda;
	double m_lambdaDerivative;
	unsigned int m_NumberOfThreads;
	char m_Evaluator;
	double m_lambdaMI;
	double m_lambdaDerivativeMI;



protected:

	MANGF();
	virtual ~MANGF();
	void PrintSelf(std::ostream & os, Indent indent) const;
	typedef MattesMutualInformationImageToImageMetric<FixedImageType,MovingImageType> MattesType;
	typedef NormalizedGradientFieldImageToImageMetric<FixedImageType,MovingImageType> NGFType;
	typedef LinearInterpolateImageFunction<FixedImageType,double > LFType;



	typedef typename NGFType::MovingNGFType MovingNGFType;
	typedef typename NGFType::FixedNGFType  FixedNGFType;




private:
	// purposely not implemented
	MANGF(const Self &);
	// purposely not implemented
	void operator=(const Self &);


	typename MattesType::Pointer m_MA;
	typename NGFType::Pointer m_NGF;
	typename LFType::Pointer m_INTERNALL_interpolator;






};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMANGF.hxx"
#endif

#endif
