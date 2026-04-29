/*=========================================================================
 * itkNMIFromMattes.hxx
 *=========================================================================*/
#ifndef __itkNMIFromMattes_hxx
#define __itkNMIFromMattes_hxx

#include "itkNMIFromMattes.h"
#include "itkImageRegionConstIterator.h"
#include <cmath>
#include <limits>

namespace itk
{

// ─────────────────────────────────────────────────────────────────────────────
// ComputeEntropies
// Reads the joint PDF already built by Mattes and computes marginal entropies.
// ─────────────────────────────────────────────────────────────────────────────
template <class TFixedImage, class TMovingImage>
void
NMIFromMattesMetric<TFixedImage, TMovingImage>
::ComputeEntropies(double & HF, double & HM, double & HFM) const
{
    constexpr double kEps = std::numeric_limits<PDFValueType>::epsilon();

    typename Superclass::JointPDFType::Pointer jointPDF = this->GetJointPDF();
    if (!jointPDF)
    {
        HF = HM = HFM = 0.0;
        return;
    }

    const unsigned int nBins = static_cast<unsigned int>(this->GetNumberOfHistogramBins());

    // Accumulate joint entropy H(F,M) and both marginals from the public joint PDF.
    // p(f) = sum_m p(f,m),  p(m) = sum_f p(f,m)  — no private member access needed.
    HFM = 0.0;
    HM  = 0.0;
    HF  = 0.0;

    std::vector<double> movingMarginal(nBins, 0.0);
    std::vector<double> fixedMarginal(nBins, 0.0);

    // Joint PDF index layout: idx[0]=movingBin, idx[1]=fixedBin
    typename Superclass::JointPDFType::IndexType jidx;
    for (unsigned int mBin = 0; mBin < nBins; ++mBin)
    {
        jidx[0] = static_cast<long>(mBin);
        for (unsigned int fBin = 0; fBin < nBins; ++fBin)
        {
            jidx[1] = static_cast<long>(fBin);
            const double pjoint = static_cast<double>(jointPDF->GetPixel(jidx));
            if (pjoint > kEps)
                HFM -= pjoint * std::log(pjoint);
            movingMarginal[mBin] += pjoint;
            fixedMarginal[fBin]  += pjoint;
        }
    }

    for (unsigned int mBin = 0; mBin < nBins; ++mBin)
    {
        const double pm = movingMarginal[mBin];
        if (pm > kEps) HM -= pm * std::log(pm);
    }
    for (unsigned int fBin = 0; fBin < nBins; ++fBin)
    {
        const double pf = fixedMarginal[fBin];
        if (pf > kEps) HF -= pf * std::log(pf);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ComputeEntropyDerivatives
// Uses the joint PDF derivatives (3-D array: param × moving-bin × fixed-bin)
// already computed by Mattes to derive grad H(M) and grad H(F,M).
//
//   grad_k H(F,M) = -sum_{f,m} (1 + log p(f,m)) * d/dk p(f,m)
//   grad_k H(M)   = -sum_m    (1 + log p_M(m))  * sum_f d/dk p(f,m)
//
// Note: d/dk p(f,m) is exactly the joint PDF derivative stored by Mattes.
// ─────────────────────────────────────────────────────────────────────────────
template <class TFixedImage, class TMovingImage>
void
NMIFromMattesMetric<TFixedImage, TMovingImage>
::ComputeEntropyDerivatives(double           HFM,
                            double           HM,
                            DerivativeType & gradHM,
                            DerivativeType & gradHFM) const
{
    constexpr double kEps = std::numeric_limits<PDFValueType>::epsilon();

    typename Superclass::JointPDFType::Pointer jointPDF = this->GetJointPDF();
    typename Superclass::JointPDFDerivativesType::Pointer jointPDFDeriv =
        this->GetJointPDFDerivatives();

    const unsigned int nParams = gradHFM.Size();
    const unsigned int nBins   = static_cast<unsigned int>(this->GetNumberOfHistogramBins());

    gradHFM.Fill(0.0);
    gradHM.Fill(0.0);

    // JointPDFDerivatives is only populated when UseExplicitPDFDerivatives=true.
    // If not available, leave derivatives at zero (value-only mode).
    if (!jointPDFDeriv || !jointPDF)
        return;

    // Pre-compute moving marginal and its log
    std::vector<double> movingMarginal(nBins, 0.0);
    {
        itk::ImageRegionConstIterator<typename Superclass::JointPDFType>
            it(jointPDF, jointPDF->GetBufferedRegion());
        it.GoToBegin();
        for (unsigned int mBin = 0; mBin < nBins; ++mBin)
            for (unsigned int fBin = 0; fBin < nBins; ++fBin, ++it)
                movingMarginal[mBin] += static_cast<double>(it.Get());
    }

    std::vector<double> logMoving(nBins, 0.0);
    for (unsigned int mBin = 0; mBin < nBins; ++mBin)
        if (movingMarginal[mBin] > kEps)
            logMoving[mBin] = std::log(movingMarginal[mBin]);

    // Joint PDF log table
    // Access via index type: [param, moving, fixed]
    typename Superclass::JointPDFDerivativesType::IndexType idx;

    for (unsigned int p = 0; p < nParams; ++p)
    {
        idx[0] = static_cast<long>(p);
        double dHFM = 0.0;
        double dHM  = 0.0;

        for (unsigned int mBin = 0; mBin < nBins; ++mBin)
        {
            idx[1] = static_cast<long>(mBin);
            double sumOverFixed = 0.0;  // sum_f d/dp p(f,m)

            for (unsigned int fBin = 0; fBin < nBins; ++fBin)
            {
                idx[2] = static_cast<long>(fBin);
                const double dpjoint = static_cast<double>(
                    jointPDFDeriv->GetPixel(idx));
                if (dpjoint == 0.0) continue;

                // grad H(F,M) term
                typename Superclass::JointPDFType::IndexType jidx;
                jidx[0] = static_cast<long>(mBin);
                jidx[1] = static_cast<long>(fBin);
                const double pjoint = static_cast<double>(jointPDF->GetPixel(jidx));
                if (pjoint > kEps)
                    dHFM -= (1.0 + std::log(pjoint)) * dpjoint;
                else
                    dHFM -= dpjoint;  // lim p->0: (1+log p)*dp → dp

                sumOverFixed += dpjoint;
            }

            // grad H(M) term
            if (movingMarginal[mBin] > kEps)
                dHM -= (1.0 + logMoving[mBin]) * sumOverFixed;
            else
                dHM -= sumOverFixed;
        }

        gradHFM[p] = dHFM;
        gradHM[p]  = dHM;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// GetValue
// ─────────────────────────────────────────────────────────────────────────────
template <class TFixedImage, class TMovingImage>
typename NMIFromMattesMetric<TFixedImage, TMovingImage>::MeasureType
NMIFromMattesMetric<TFixedImage, TMovingImage>
::GetValue(const ParametersType & parameters) const
{
    // Let Mattes fill the joint PDF
    this->Superclass::GetValue(parameters);

    double HF, HM, HFM;
    this->ComputeEntropies(HF, HM, HFM);

    if (HFM < 1e-12)
        return static_cast<MeasureType>(0.0);

    // NMI is maximised → negate so LBFGS-B minimises
    const double nmi = (HF + HM) / HFM;
    return static_cast<MeasureType>(-nmi);
}

// ─────────────────────────────────────────────────────────────────────────────
// GetValueAndDerivative
// ─────────────────────────────────────────────────────────────────────────────
template <class TFixedImage, class TMovingImage>
void
NMIFromMattesMetric<TFixedImage, TMovingImage>
::GetValueAndDerivative(const ParametersType & parameters,
                        MeasureType          & value,
                        DerivativeType       & derivative) const
{
    // Let Mattes fill the joint PDF and its derivatives in one pass
    MeasureType    unusedMI;
    DerivativeType unusedDeriv(parameters.Size());
    this->Superclass::GetValueAndDerivative(parameters, unusedMI, unusedDeriv);

    // Compute entropies from the now-filled joint PDF
    double HF, HM, HFM;
    this->ComputeEntropies(HF, HM, HFM);

    if (HFM < 1e-12)
    {
        value = static_cast<MeasureType>(0.0);
        derivative.Fill(0.0);
        return;
    }

    const double nmi = (HF + HM) / HFM;
    value = static_cast<MeasureType>(-nmi);

    // Compute entropy derivatives
    DerivativeType gradHM(parameters.Size());
    DerivativeType gradHFM(parameters.Size());
    this->ComputeEntropyDerivatives(HFM, HM, gradHM, gradHFM);

    // grad NMI = ( grad H(M) - NMI * grad H(F,M) ) / H(F,M)
    // Negate because we minimise -NMI
    derivative.SetSize(parameters.Size());
    for (unsigned int p = 0; p < parameters.Size(); ++p)
        derivative[p] = static_cast<typename DerivativeType::ValueType>(
            -(gradHM[p] - nmi * gradHFM[p]) / HFM);
}

// ─────────────────────────────────────────────────────────────────────────────
// GetDerivative
// ─────────────────────────────────────────────────────────────────────────────
template <class TFixedImage, class TMovingImage>
void
NMIFromMattesMetric<TFixedImage, TMovingImage>
::GetDerivative(const ParametersType & parameters,
                DerivativeType       & derivative) const
{
    MeasureType value;
    this->GetValueAndDerivative(parameters, value, derivative);
}

} // end namespace itk

#endif // __itkNMIFromMattes_hxx
