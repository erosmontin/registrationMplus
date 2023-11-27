/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkHistogramToTextureFeaturesFilterExtended_h
#define __itkHistogramToTextureFeaturesFilterExtended_h

#include "itkHistogram.h"
#include "itkMacro.h"
#include "itkProcessObject.h"
#include "itkSimpleDataObjectDecorator.h"

namespace itk
{
namespace Statistics
{
/** \class HistogramToTextureFeaturesFilter
*  \brief This class computes texture feature coefficients from a grey level
* co-occurrence matrix.
*
* This class computes features that summarize image texture, given a grey level
* co-occurrence matrix (generated by a ScalarImageToCooccurrenceMatrixFilter
* or related class).
*
* The features calculated are as follows (where \f$ g(i, j) \f$ is the element in
* cell i, j of a a normalized GLCM):
*
* "Energy" \f$ = f_1 = \sum_{i,j}g(i, j)^2 \f$
*
* "Entropy" \f$ = f_2 = -\sum_{i,j}g(i, j) \log_2 g(i, j)\f$, or 0 if \f$g(i, j) = 0\f$
*
* "Correlation" \f$ = f_3 = \sum_{i,j}\frac{(i - \mu)(j - \mu)g(i, j)}{\sigma^2} \f$
*
* "Difference Moment" \f$= f_4 = \sum_{i,j}\frac{1}{1 + (i - j)^2}g(i, j) \f$
*
* "Inertia" \f$ = f_5 = \sum_{i,j}(i - j)^2g(i, j) \f$ (sometimes called "contrast.")
*
* "Cluster Shade" \f$ = f_6 = \sum_{i,j}((i - \mu) + (j - \mu))^3 g(i, j) \f$
*
* "Cluster Prominence" \f$ = f_7 = \sum_{i,j}((i - \mu) + (j - \mu))^4 g(i, j) \f$
*
* "Haralick's Correlation" \f$ = f_8 = \frac{\sum_{i,j}(i, j) g(i, j) -\mu_t^2}{\sigma_t^2} \f$
* where \f$\mu_t\f$ and \f$\sigma_t\f$ are the mean and standard deviation of the row
* (or column, due to symmetry) sums.
*
* Above, \f$ \mu =  \f$ (weighted pixel average) \f$ = \sum_{i,j}i \cdot g(i, j) =
* \sum_{i,j}j \cdot g(i, j) \f$ (due to matrix summetry), and
*
* \f$ \sigma =  \f$ (weighted pixel variance) \f$ = \sum_{i,j}(i - \mu)^2 \cdot g(i, j) =
* \sum_{i,j}(j - \mu)^2 \cdot g(i, j)  \f$  (due to matrix summetry)
*
* A good texture feature set to use is the Conners, Trivedi and Harlow set:
* features 1, 2, 4, 5, 6, and 7. There is some correlation between the various
* features, so using all of them at the same time is not necessarialy a good idea.
*
* NOTA BENE: The input histogram will be forcably normalized!
* This algorithm takes three passes through the input
* histogram if the histogram was already normalized, and four if not.
*
* Web references:
*
* http://www.cssip.uq.edu.au/meastex/www/algs/algs/algs.html
* http://www.ucalgary.ca/~mhallbey/texture/texture_tutorial.html
*
* Print references:
*
* Haralick, R.M., K. Shanmugam and I. Dinstein. 1973.  Textural Features for
* Image Classification. IEEE Transactions on Systems, Man and Cybernetics.
* SMC-3(6):610-620.
*
* Haralick, R.M. 1979. Statistical and Structural Approaches to Texture.
* Proceedings of the IEEE, 67:786-804.
*
* R.W. Conners and C.A. Harlow. A Theoretical Comaprison of Texture Algorithms.
* IEEE Transactions on Pattern Analysis and Machine Intelligence,  2:204-222, 1980.
*
* R.W. Conners, M.M. Trivedi, and C.A. Harlow. Segmentation of a High-Resolution
* Urban Scene using Texture  Operators. Computer Vision, Graphics and Image
* Processing, 25:273-310,  1984.
*
* \sa ScalarImageToCooccurrenceMatrixFilter
* \sa ScalarImageToTextureFeaturesFilter
*
* Author: Zachary Pincus
* \ingroup ITKStatistics
*/

template< typename THistogram >
class HistogramToTextureFeaturesFilterExtended:public ProcessObject
{
public:
  /** Standard typedefs */
  typedef HistogramToTextureFeaturesFilterExtended Self;
  typedef ProcessObject                    Superclass;
  typedef SmartPointer< Self >             Pointer;
  typedef SmartPointer< const Self >       ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(HistogramToTextureFeaturesFilterExtended, ProcessObject);

  /** standard New() method support */
  itkNewMacro(Self);



  typedef THistogram                                    HistogramType;
  typedef typename HistogramType::Pointer               HistogramPointer;
  typedef typename HistogramType::ConstPointer          HistogramConstPointer;
  typedef typename HistogramType::MeasurementType       MeasurementType;
  typedef typename HistogramType::MeasurementVectorType MeasurementVectorType;
  typedef typename HistogramType::IndexType             IndexType;
  typedef typename HistogramType::AbsoluteFrequencyType AbsoluteFrequencyType;
  typedef typename HistogramType::RelativeFrequencyType RelativeFrequencyType;

  typedef typename HistogramType::TotalAbsoluteFrequencyType
  TotalAbsoluteFrequencyType;

  typedef typename HistogramType::TotalRelativeFrequencyType
  TotalRelativeFrequencyType;

  /** Container to hold relative frequencies of the histogram */
  typedef std::vector< RelativeFrequencyType > RelativeFrequencyContainerType;

  /** Method to Set/Get the input Histogram */
  using Superclass::SetInput;
  void SetInput(const HistogramType *histogram);

  const HistogramType * GetInput() const;

  /** Smart Pointer type to a DataObject. */
  typedef DataObject::Pointer DataObjectPointer;

  /** Type of DataObjects used for scalar outputs */
  typedef SimpleDataObjectDecorator< MeasurementType > MeasurementObjectType;

  /** Return energy texture value. */
  MeasurementType GetEnergy() const;

  const MeasurementObjectType * GetEnergyOutput() const;

  /** Return entropy texture value. */
  MeasurementType GetEntropy() const;

  const MeasurementObjectType * GetEntropyOutput() const;

  /** return correlation texture value. */
  MeasurementType GetCorrelation() const;

  const MeasurementObjectType * GetCorrelationOutput() const;

  /** Return inverse difference moment texture value. */
  MeasurementType GetInverseDifferenceMoment() const;

  const MeasurementObjectType * GetInverseDifferenceMomentOutput() const;

  /** Return inertia texture value. */
  MeasurementType GetInertia() const;

  const MeasurementObjectType * GetInertiaOutput() const;

  /** Return cluster shade texture value. */
  MeasurementType GetClusterShade() const;

  const MeasurementObjectType * GetClusterShadeOutput() const;

  /** Return cluster prominence texture value. */
  MeasurementType GetClusterProminence() const;

  const MeasurementObjectType * GetClusterProminenceOutput() const;

  /** Return Haralick correlation texture value. */
  MeasurementType GetHaralickCorrelation() const;

  const MeasurementObjectType * GetHaralickCorrelationOutput() const;

  /**The new part */
  /** Return Haralick correlation texture value. */
  MeasurementType GetContrast() const;

  const MeasurementObjectType * GetContrastOutput() const;

  /** Return Autocorrelation texture value. */
  MeasurementType GetAutocorrelation() const;

  const MeasurementObjectType * GetAutocorrelationOutput() const;

  /** Return ClusterTendency texture value. */
  MeasurementType GetClusterTendency() const;

  const MeasurementObjectType * GetClusterTendencyOutput() const;


  /** Return Dissimilarity texture value. */
  MeasurementType GetDissimilarity() const;

  const MeasurementObjectType * GetDissimilarityOutput() const;

  /** Return Homogeneity texture value. */
  MeasurementType GetHomogeneity() const;

  const MeasurementObjectType * GetHomogeneityOutput() const;

  /** Return Homogeneity2 texture value. */
  MeasurementType GetHomogeneity2() const;

  const MeasurementObjectType * GetHomogeneity2Output() const;


  /** Return InverseDifferenceMoment2 texture value. */
  MeasurementType GetInverseDifferenceMoment2() const;

  const MeasurementObjectType * GetInverseDifferenceMoment2Output() const;

  /** Return InverseDifferenceMoment2 texture value. */
  MeasurementType GetDifferenceEntropy() const;

  const MeasurementObjectType * GetDifferenceEntropyOutput() const;

  /** Return InverseDifferenceMoment2 texture value. */
  MeasurementType GetInverseVariance() const;

  const MeasurementObjectType * GetInverseVarianceOutput() const;


  /** Return InverseDifferenceMoment2 texture value. */
  MeasurementType GetInformationalMeasureOfCorrelation1() const;

  const MeasurementObjectType * GetInformationalMeasureOfCorrelation1Output() const;

  /** Return InverseDifferenceMoment2 texture value. */
   MeasurementType GetInformationalMeasureOfCorrelation2() const;

   const MeasurementObjectType * GetInformationalMeasureOfCorrelation2Output() const;

   /** Return InverseDifferenceMoment2 texture value. */
   MeasurementType GetMaxProbability() const;
   const MeasurementObjectType * GetMaxProbabilityOutput() const;

   /** Return SumAverage texture value. */
   MeasurementType GetSumAverage() const;
   const MeasurementObjectType * GetSumAverageOutput() const;


   /** Return SumEntropy texture value. */
   MeasurementType GetSumEntropy() const;
   const MeasurementObjectType * GetSumEntropyOutput() const;

//   /** Return SumVariance texture value. */
//   MeasurementType GetSumVariance() const;
//   const MeasurementObjectType * GetSumVarianceOutput() const;
//
//   /** Return Variance texture value. */
//   MeasurementType GetVariance() const;
//   const MeasurementObjectType * GetVarianceOutput() const;


  /** Texture feature types */
  typedef enum {
    Energy,
    Entropy,
    Correlation,
    InverseDifferenceMoment,
    Inertia,
    ClusterShade,
    ClusterProminence,
    HaralickCorrelation,
    Contrast,
    Autocorrelation,
    ClusterTendency,
    Dissimilarity,
    Homogeneity,
    Homogeneity2,
    InverseDifferenceMoment2,
    DifferenceEntropy,
    InverseVariance,
    InformationalMeasureOfCorrelation1,
    InformationalMeasureOfCorrelation2,
    MaxProbability,
    SumAverage,
    SumEntropy,
    InvalidFeatureName
    }  TextureFeatureName;

    //SumVariance,
     //Variance,

  /** convenience method to access the texture values */
  MeasurementType GetFeature(TextureFeatureName name);

protected:
  HistogramToTextureFeaturesFilterExtended();
  ~HistogramToTextureFeaturesFilterExtended() {}
  virtual void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  /** Make a DataObject to be used for output output. */
  typedef ProcessObject::DataObjectPointerArraySizeType DataObjectPointerArraySizeType;
  using Superclass::MakeOutput;
  virtual DataObjectPointer MakeOutput(DataObjectPointerArraySizeType) ITK_OVERRIDE;

  virtual void GenerateData() ITK_OVERRIDE;

private:
  HistogramToTextureFeaturesFilterExtended(const Self &); //purposely not implemented
  void operator=(const Self &);                   //purposely not implemented

  void ComputeMeansAndVariances(double & pixelMean, double & marginalMean,
                                double & marginalDevSquared, double & pixelVariance);



  RelativeFrequencyContainerType m_RelativeFrequencyContainer;
};
} // end of namespace Statistics
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkHistogramToTextureFeaturesFilterExtended.hxx"
#endif

#endif
