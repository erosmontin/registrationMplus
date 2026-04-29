/*=========================================================================
 * itkNMIFromMattes.h
 *
 * Normalized Mutual Information metric with analytical derivatives, built
 * on top of MattesMutualInformationImageToImageMetric.
 *
 * NMI = ( H(F) + H(M) ) / H(F,M)
 *
 * where H(F), H(M), H(F,M) are the marginal and joint Shannon entropies
 * estimated with the same Parzen-window joint PDF that Mattes MI uses.
 *
 * Derivative:
 *   grad NMI = ( grad H(M)  -  NMI * grad H(F,M) ) / H(F,M)
 *
 * H(F) is constant (fixed image does not change) so grad H(F) = 0.
 * grad H(M) and grad H(F,M) are derived analytically from the joint PDF
 * and its parameter derivatives — the same arrays already computed by the
 * Mattes base class.
 *
 * This gives LBFGS-B-compatible value/gradient pairs at the same cost as
 * standard Mattes MI (one forward pass, same Parzen weights).
 *=========================================================================*/
#ifndef __itkNMIFromMattes_h
#define __itkNMIFromMattes_h

#include "itkMattesMutualInformationImageToImageMetric.h"

namespace itk
{

template <class TFixedImage, class TMovingImage>
class ITK_EXPORT NMIFromMattesMetric
    : public MattesMutualInformationImageToImageMetric<TFixedImage, TMovingImage>
{
public:
    typedef NMIFromMattesMetric                                          Self;
    typedef MattesMutualInformationImageToImageMetric<TFixedImage,
                                                     TMovingImage>      Superclass;
    typedef SmartPointer<Self>                                           Pointer;
    typedef SmartPointer<const Self>                                     ConstPointer;

    itkNewMacro(Self);
    itkTypeMacro(NMIFromMattesMetric,
                 MattesMutualInformationImageToImageMetric);

    typedef typename Superclass::MeasureType    MeasureType;
    typedef typename Superclass::DerivativeType DerivativeType;
    typedef typename Superclass::ParametersType ParametersType;
    typedef typename Superclass::PDFValueType   PDFValueType;

    /** Compute NMI value only. */
    MeasureType GetValue(const ParametersType & parameters) const ITK_OVERRIDE;

    /** Compute NMI value and its analytical derivative simultaneously. */
    void GetValueAndDerivative(const ParametersType & parameters,
                               MeasureType          & value,
                               DerivativeType       & derivative) const ITK_OVERRIDE;

    /** Derivative-only convenience wrapper. */
    void GetDerivative(const ParametersType & parameters,
                       DerivativeType       & derivative) const ITK_OVERRIDE;

protected:
    NMIFromMattesMetric() {}
    virtual ~NMIFromMattesMetric() {}

private:
    NMIFromMattesMetric(const Self &);
    void operator=(const Self &);

    /** Compute H(F), H(M), H(F,M) from the joint PDF that Mattes has already
     *  filled.  Called after the base-class pass. */
    void ComputeEntropies(double & HF, double & HM, double & HFM) const;

    /** Compute grad H(M) and grad H(F,M) from the joint PDF derivatives.
     *  Called after GetValueAndDerivative on the base class. */
    void ComputeEntropyDerivatives(double             HFM,
                                   double             HM,
                                   DerivativeType   & gradHM,
                                   DerivativeType   & gradHFM) const;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkNMIFromMattes.hxx"
#endif

#endif // __itkNMIFromMattes_h
