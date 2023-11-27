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
#ifndef __itkHistogramToTextureFeaturesFilterExtended_hxx
#define __itkHistogramToTextureFeaturesFilterExtended_hxx

#include "itkHistogramToTextureFeaturesFilterExtended.h"

#include "itkNumericTraits.h"
#include "vnl/vnl_math.h"
#include "itkMath.h"


#include <iostream>     // std::cout
#include <algorithm>    // std::fill
#include <vector>       // std::vector

//#include "mymath.h"

namespace itk
{
namespace Statistics
{
//constructor
template< typename THistogram >
HistogramToTextureFeaturesFilterExtended< THistogram >::HistogramToTextureFeaturesFilterExtended(void)
{


	this->ProcessObject::SetNumberOfRequiredInputs(1);


	// allocate the data objects for the outputs which are
	// just decorators real types +2 rispetto a quelle che miprendo
	for ( int i = 0; i < 23; ++i )
	{
		this->ProcessObject::SetNthOutput( i, this->MakeOutput(i) );
	}
}

template< typename THistogram >
void
HistogramToTextureFeaturesFilterExtended< THistogram >
::SetInput(const HistogramType *histogram)
 {
	this->ProcessObject::SetNthInput( 0, const_cast< HistogramType * >( histogram ) );
 }

template< typename THistogram >
const typename
HistogramToTextureFeaturesFilterExtended< THistogram >::HistogramType *
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetInput() const
 {
	return itkDynamicCastInDebugMode< const HistogramType * >( this->GetPrimaryInput() );

 }

template< typename THistogram >
typename
HistogramToTextureFeaturesFilterExtended< THistogram >::DataObjectPointer
HistogramToTextureFeaturesFilterExtended< THistogram >
::MakeOutput( DataObjectPointerArraySizeType itkNotUsed(idx) )
 {
	return MeasurementObjectType::New().GetPointer();
 }

template< typename THistogram >
void
HistogramToTextureFeaturesFilterExtended< THistogram >::GenerateData(void)
{
	typedef typename HistogramType::ConstIterator HistogramIterator;

	const HistogramType *inputHistogram = this->GetInput();

	//Normalize the absolute frequencies and populate the relative frequency
	//container
	TotalRelativeFrequencyType totalFrequency =
			static_cast< TotalRelativeFrequencyType >( inputHistogram->GetTotalFrequency() );

	m_RelativeFrequencyContainer.clear();

	for ( HistogramIterator hit = inputHistogram->Begin();
			hit != inputHistogram->End(); ++hit )
	{
		AbsoluteFrequencyType frequency = hit.GetFrequency();
		RelativeFrequencyType relativeFrequency =  frequency / totalFrequency;
		m_RelativeFrequencyContainer.push_back(relativeFrequency);
	}

	// Now get the various means and variances. This is takes two passes
	// through the histogram.
	double pixelMean;
	double marginalMean;
	double marginalDevSquared;
	double pixelVariance;

	double m_Pmean0;
	double m_Pmean1;
	double m_Pstd0;
	double m_Pstd1;
	double m_H0,m_H1,m_H,m_HXY1,m_HXY2;

	this->ComputeMeansAndVariances(pixelMean, marginalMean, marginalDevSquared,pixelVariance);





	// Finally compute the texture features. Another one pass.
	MeasurementType energy      = NumericTraits< MeasurementType >::ZeroValue();
	MeasurementType entropy     = NumericTraits< MeasurementType >::ZeroValue();
	MeasurementType correlation = NumericTraits< MeasurementType >::ZeroValue();

	MeasurementType inverseDifferenceMoment      =
			NumericTraits< MeasurementType >::ZeroValue();

	MeasurementType inertia             = NumericTraits< MeasurementType >::ZeroValue();
	MeasurementType clusterShade        = NumericTraits< MeasurementType >::ZeroValue();
	MeasurementType clusterProminence   = NumericTraits< MeasurementType >::ZeroValue();
	MeasurementType haralickCorrelation = NumericTraits< MeasurementType >::ZeroValue();

	MeasurementType contrast = NumericTraits< MeasurementType >::ZeroValue();
	MeasurementType autocorrelation = NumericTraits< MeasurementType >::ZeroValue();
	MeasurementType clusterTendency = NumericTraits< MeasurementType >::ZeroValue();

	MeasurementType dissimilarity = NumericTraits< MeasurementType >::ZeroValue();
	MeasurementType homogeneity = NumericTraits< MeasurementType >::ZeroValue();
	MeasurementType homogeneity2 = NumericTraits< MeasurementType >::ZeroValue();

	MeasurementType inverseDifferenceMoment2 = NumericTraits< MeasurementType >::ZeroValue();
	MeasurementType differenceEntropy = NumericTraits< MeasurementType >::ZeroValue();

	MeasurementType inverseVariance = NumericTraits< MeasurementType >::ZeroValue();

	MeasurementType informationalMeasureOfCorrelation1 = NumericTraits< MeasurementType >::ZeroValue();
	MeasurementType informationalMeasureOfCorrelation2 = NumericTraits< MeasurementType >::ZeroValue();


MeasurementType    maxProbability = NumericTraits< MeasurementType >::ZeroValue();
MeasurementType    sumAverage = NumericTraits< MeasurementType >::ZeroValue();
MeasurementType    sumEntropy = NumericTraits< MeasurementType >::ZeroValue();
//MeasurementType    sumVariance = NumericTraits< MeasurementType >::ZeroValue();
//MeasurementType    variance = NumericTraits< MeasurementType >::ZeroValue();




	double pixelVarianceSquared = pixelVariance * pixelVariance;
	// Variance is only used in correlation. If variance is 0, then
	//   (index[0] - pixelMean) * (index[1] - pixelMean)
	// should be zero as well. In this case, set the variance to 1. in
	// order to avoid NaN correlation.
	if( Math::FloatAlmostEqual( pixelVarianceSquared, 0.0, 4, 2*NumericTraits<double>::epsilon() ) )
	{
		pixelVarianceSquared = 1.;
	}
	const double log2 = std::log(2.0);

	typename RelativeFrequencyContainerType::const_iterator rFreqIterator =
			m_RelativeFrequencyContainer.begin();

	typename RelativeFrequencyContainerType::const_iterator r2FreqIterator =
			m_RelativeFrequencyContainer.begin();

	typename RelativeFrequencyContainerType::const_iterator r3FreqIterator =
			m_RelativeFrequencyContainer.begin();

	double N= inputHistogram->GetSize(0)*inputHistogram->GetSize(1);

	std::vector<double> P (inputHistogram->GetSize(0)+inputHistogram->GetSize(1));
	std::vector<double> Pxy (inputHistogram->GetSize(0)+inputHistogram->GetSize(1));
	std::vector<double> Px (std::sqrt(N));
	std::vector<double> Py (std::sqrt(N));


	std::vector<double> Pxmy (std::sqrt(N));
	std::vector<double> Pxpy (2 * std::sqrt(N));

	std::fill(P.begin(),P.end(),0);
	std::fill(Pxy.begin(),Pxy.end(),0);
	std::fill(Px.begin(),Px.end(),0);
	std::fill(Py.begin(),Py.end(),0);
	std::fill(Pxmy.begin(),Pxmy.end(),0);
	std::fill(Pxpy.begin(),Pxpy.end(),0);

	//**intialize the values*/
	IndexType o;

	for ( HistogramIterator hit = inputHistogram->Begin();
			hit != inputHistogram->End(); ++hit )
	{

		RelativeFrequencyType frequency = *r2FreqIterator;
		++r2FreqIterator;

		o=hit.GetIndex();
		//		Px.at(P[0])+=hit.GetFrequency()/totalFrequency;
		//		Py.at(P[1])+=hit.GetFrequency()/totalFrequency;
		P.push_back(frequency);
		Px.at(o[1])+=frequency;
		Py.at(o[0])+=frequency;
		Pxmy.at(std::abs(o[0]-o[1]))+=frequency;
		Pxpy.at(o[0]+o[1])+=frequency;

	}

	double d;
	for ( HistogramIterator hit = inputHistogram->Begin();
			hit != inputHistogram->End(); ++hit )
	{

		RelativeFrequencyType frequency = *r3FreqIterator;
		++r3FreqIterator;

		o=hit.GetIndex();
		d=Px.at(o[0]) * Py.at(o[1]);
		Pxy.push_back(d);
		m_HXY1 += ( d > 1e-6 ) ?  - frequency * std::log( d ) / std::log( 2.0 ) :0;

	}

	m_Pmean0=CALC::mean(Px);
	m_Pmean1=CALC::mean(Py);
	m_Pstd0=CALC::std(Px);
	m_Pstd1=CALC::std(Py);

	m_H0=CALC::entropy(Px);
	m_H1=CALC::entropy(Py);
	m_H=CALC::entropy(P);
	m_HXY2=CALC::entropy(Pxy);

	//	std::cout<<"\nMean0:"<<m_H1<<"\nMean1"<<m_H1<<"\nH"<<m_H<<"\nXY2"<<m_HXY2<<"\nXY1"<<m_HXY1<<std::endl;
	//	std::cout<<"\nMean0:"<<m_Pmean0<<"\nMean1"<<m_Pmean1<<"\nstd0"<<m_Pstd0<<"\nstd1"<<m_Pstd1<<std::endl;

	double CS=0.0;

	for ( HistogramIterator hit = inputHistogram->Begin();
			hit != inputHistogram->End(); ++hit )
	{
		RelativeFrequencyType frequency = *rFreqIterator;
		++rFreqIterator;
		if ( frequency == 0 )
		{
			continue; // no use doing these calculations if we're just multiplying by
			//zero.
		}

		IndexType index = inputHistogram->GetIndex( hit.GetInstanceIdentifier() );

		haralickCorrelation += index[0] * index[1] * frequency;

		autocorrelation=haralickCorrelation;

		contrast+= std::pow(std::abs(index[0]-index[1]),2) * frequency;

		CS=( index[0] +index[1] - m_Pmean0  - m_Pmean1 );

		clusterShade += std::pow(CS , 3.0 ) * frequency;
		clusterProminence +=std::pow( CS, 4.0 ) * frequency;
		clusterTendency += std::pow(CS, 2.0 ) * frequency;

		dissimilarity+=std::abs(index[0]-index[1])*frequency;

		energy += frequency * frequency;

		//		entropy -= ( frequency > 1e-6 ) ? frequency *std::log(frequency) / log2:0;

		entropy -= ( frequency > 1e-6 ) ? frequency *std::log(frequency) / log2 :0;

		homogeneity+= frequency / ( 1.0 + std::abs(index[0] - index[1]));

		homogeneity2+=frequency / ( 1.0 + std::pow(std::abs( index[0] - index[1] ) , 2));
		inertia += ( index[0] - index[1] ) * ( index[0] - index[1] ) * frequency;

		inverseDifferenceMoment += frequency/ ( 1.0 + std::pow(std::abs(index[0] - index[1])/(N*N),2) );

		inverseDifferenceMoment2 += frequency/ ( 1.0 + std::abs(index[0] - index[1])/(N) );

		differenceEntropy += frequency	/ ( 1.0 + ( ( ( index[0] - index[1] ) * ( index[0] - index[1] ) )/(N*N)) );

		double q=  std::pow( std::abs(index[0] - index[1]),2);
		inverseVariance += (index[0]== index[1]) ? 0 : frequency/q;

	}

	haralickCorrelation = ( haralickCorrelation - marginalMean * marginalMean )
                        																/ marginalDevSquared;

	correlation = (autocorrelation - (m_Pmean0 * m_Pmean1))/(m_Pstd0 *m_Pstd1);

	differenceEntropy = CALC::entropy(Pxmy);

	std::vector<double> ff(2);
	ff.push_back(m_H0);ff.push_back(m_H1);
	informationalMeasureOfCorrelation1=m_H-m_HXY1/ CALC::max(ff);

	informationalMeasureOfCorrelation2= std::sqrt(1- std::exp(-2*(m_HXY2-m_H)));


	maxProbability = CALC::max(P);

	std::vector<double>::iterator cP;
	int u=0;
	for ( cP = Pxpy.begin(); cP != Pxpy.end(); ++cP ){
			sumAverage+= *cP * (double)u;
			sumEntropy += ( *cP > 1e-6 ) ?  - *cP * std::log( *cP ) / std::log( 2.0 ) :0;
			u++;
			};



	MeasurementObjectType *energyOutputObject =
			static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(0) );
	energyOutputObject->Set(energy);

	MeasurementObjectType *entropyOutputObject =
			static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(1) );
	entropyOutputObject->Set(entropy);

	MeasurementObjectType *correlationOutputObject =
			static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(2) );
	correlationOutputObject->Set(correlation);

	MeasurementObjectType *inverseDifferenceMomentOutputObject =
			static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(3) );
	inverseDifferenceMomentOutputObject->Set(inverseDifferenceMoment);

	MeasurementObjectType *inertiaOutputObject =
			static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(4) );
	inertiaOutputObject->Set(inertia);

	MeasurementObjectType *clusterShadeOutputObject =
			static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(5) );
	clusterShadeOutputObject->Set(clusterShade);

	MeasurementObjectType *clusterProminenceOutputObject =
			static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(6) );
	clusterProminenceOutputObject->Set(clusterProminence);

	MeasurementObjectType *haralickCorrelationOutputObject =
			static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(7) );
	haralickCorrelationOutputObject->Set(haralickCorrelation);

	/*news*/
	MeasurementObjectType *contrastOutputObject =
			static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(8) );
	contrastOutputObject->Set(contrast);

	MeasurementObjectType *autocorrelationOutputObject =
			static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(9) );
	autocorrelationOutputObject->Set(autocorrelation);

	MeasurementObjectType *clusterTendencyOutputObject =
			static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(10) );
	clusterTendencyOutputObject->Set(clusterTendency);

	MeasurementObjectType *dissimilarityOutputObject =
			static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(11) );
	dissimilarityOutputObject->Set(dissimilarity);


	MeasurementObjectType *homogeneityOutputObject =
			static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(12) );
	homogeneityOutputObject->Set(homogeneity);

	MeasurementObjectType *homogeneity2OutputObject =
			static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(13) );
	homogeneity2OutputObject->Set(homogeneity2);

	MeasurementObjectType *inverseDifferenceMoment2OutputObject =
			static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(14) );
	inverseDifferenceMoment2OutputObject->Set(inverseDifferenceMoment2);

	MeasurementObjectType *differenceEntropyOutputObject =
			static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(15) );
	differenceEntropyOutputObject->Set(differenceEntropy);

	MeasurementObjectType *inverseVarianceOutputObject =
			static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(16) );
	inverseVarianceOutputObject->Set(inverseVariance);

	MeasurementObjectType *informationalMeasureOfCorrelation1OutputObject =
			static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(17) );
	informationalMeasureOfCorrelation1OutputObject->Set(informationalMeasureOfCorrelation1);

	MeasurementObjectType *informationalMeasureOfCorrelation2OutputObject =
			static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(18) );
	informationalMeasureOfCorrelation2OutputObject->Set(informationalMeasureOfCorrelation2);


	MeasurementObjectType *maxProbabilityOutputObject =
			static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(19) );
	maxProbabilityOutputObject->Set(maxProbability);

	MeasurementObjectType *sumAverageOutputObject =
			static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(20) );
	sumAverageOutputObject->Set(sumAverage);

	MeasurementObjectType *sumEntropyOutputObject =
			static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(21) );
	sumEntropyOutputObject->Set(sumEntropy);



}

template< typename THistogram >
void
HistogramToTextureFeaturesFilterExtended< THistogram >::ComputeMeansAndVariances(double & pixelMean,
		double & marginalMean,
		double & marginalDevSquared,
		double & pixelVariance)
		{
	// This function takes two passes through the histogram and two passes through
	// an array of the same length as a histogram axis. This could probably be
	// cleverly compressed to one pass, but it's not clear that that's necessary.

	typedef typename HistogramType::ConstIterator HistogramIterator;

	const HistogramType *inputHistogram =  this->GetInput();

	// Initialize everything
	typename HistogramType::SizeValueType binsPerAxis = inputHistogram->GetSize(0);
	double *marginalSums = new double[binsPerAxis];
	for ( double *ms_It = marginalSums;
			ms_It < marginalSums + binsPerAxis; ms_It++ )
	{
		*ms_It = 0;
	}
	pixelMean = 0;

	typename RelativeFrequencyContainerType::const_iterator rFreqIterator =
			m_RelativeFrequencyContainer.begin();

	// Ok, now do the first pass through the histogram to get the marginal sums
	// and compute the pixel mean
	HistogramIterator hit = inputHistogram->Begin();
	while ( hit != inputHistogram->End() )
	{
		RelativeFrequencyType frequency = *rFreqIterator;
		IndexType             index = inputHistogram->GetIndex( hit.GetInstanceIdentifier() );
		pixelMean += index[0] * frequency;
		marginalSums[index[0]] += frequency;
		++hit;
		++rFreqIterator;
	}

	/*  Now get the mean and deviaton of the marginal sums.
      Compute incremental mean and SD, a la Knuth, "The  Art of Computer
      Programming, Volume 2: Seminumerical Algorithms",  section 4.2.2.
      Compute mean and standard deviation using the recurrence relation:
      M(1) = x(1), M(k) = M(k-1) + (x(k) - M(k-1) ) / k
      S(1) = 0, S(k) = S(k-1) + (x(k) - M(k-1)) * (x(k) - M(k))
      for 2 <= k <= n, then
      sigma = std::sqrt(S(n) / n) (or divide by n-1 for sample SD instead of
      population SD).
	 */
	marginalMean = marginalSums[0];
	marginalDevSquared = 0;
	for ( unsigned int arrayIndex = 1; arrayIndex < binsPerAxis; arrayIndex++ )
	{
		int    k = arrayIndex + 1;
		double M_k_minus_1 = marginalMean;
		double S_k_minus_1 = marginalDevSquared;
		double x_k = marginalSums[arrayIndex];

		double M_k = M_k_minus_1 + ( x_k - M_k_minus_1 ) / k;
		double S_k = S_k_minus_1 + ( x_k - M_k_minus_1 ) * ( x_k - M_k );

		marginalMean = M_k;
		marginalDevSquared = S_k;
	}
	marginalDevSquared = marginalDevSquared / binsPerAxis;

	rFreqIterator = m_RelativeFrequencyContainer.begin();
	// OK, now compute the pixel variances.
	pixelVariance = 0;
	for ( hit = inputHistogram->Begin(); hit != inputHistogram->End(); ++hit )
	{
		RelativeFrequencyType frequency = *rFreqIterator;
		IndexType             index = inputHistogram->GetIndex( hit.GetInstanceIdentifier() );
		pixelVariance += ( index[0] - pixelMean ) * ( index[0] - pixelMean ) * frequency;
		++rFreqIterator;
	}

	delete[] marginalSums;
		}

template< typename THistogram >
const
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementObjectType *
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetEnergyOutput() const
 {
	return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(0) );
 }

template< typename THistogram >
const
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementObjectType *
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetEntropyOutput() const
 {
	return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(1) );
 }

template< typename THistogram >
const
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementObjectType *
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetCorrelationOutput() const
 {
	return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(2) );
 }

template< typename THistogram >
const
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementObjectType *
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetInverseDifferenceMomentOutput() const
 {
	return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(3) );
 }

template< typename THistogram >
const
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementObjectType *
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetInertiaOutput() const
 {
	return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(4) );
 }

template< typename THistogram >
const
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementObjectType *
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetClusterShadeOutput() const
 {
	return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(5) );
 }

template< typename THistogram >
const
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementObjectType *
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetClusterProminenceOutput() const
 {
	return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(6) );
 }

template< typename THistogram >
const
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementObjectType *
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetHaralickCorrelationOutput() const
 {
	return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(7) );
 }

/** news */
template< typename THistogram >
const
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementObjectType *
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetContrastOutput() const
 {
	return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(8) );
 }

template< typename THistogram >
const
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementObjectType *
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetAutocorrelationOutput() const
 {
	return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(9) );
 }

template< typename THistogram >
const
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementObjectType *
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetClusterTendencyOutput() const
 {
	return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(10) );
 }


template< typename THistogram >
const
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementObjectType *
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetDissimilarityOutput() const
 {
	return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(11) );
 }


template< typename THistogram >
const
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementObjectType *
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetHomogeneityOutput() const
 {
	return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(12) );
 }

template< typename THistogram >
const
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementObjectType *
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetHomogeneity2Output() const
 {
	return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(13) );
 }



template< typename THistogram >
const
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementObjectType *
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetInverseDifferenceMoment2Output() const
 {
	return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(14) );
 }

template< typename THistogram >
const
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementObjectType *
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetDifferenceEntropyOutput() const
 {
	return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(15) );
 }

template< typename THistogram >
const
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementObjectType *
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetInverseVarianceOutput() const
 {
	return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(16) );
 }



template< typename THistogram >
const
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementObjectType *
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetInformationalMeasureOfCorrelation1Output() const
 {
	return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(17) );
 }

template< typename THistogram >
const
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementObjectType *
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetInformationalMeasureOfCorrelation2Output() const
 {
	return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(18) );
 }


template< typename THistogram >
const
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementObjectType *
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetMaxProbabilityOutput() const
 {
	return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(19) );
 }

template< typename THistogram >
const
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementObjectType *
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetSumAverageOutput() const
 {
	return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(20) );
 }

template< typename THistogram >
const
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementObjectType *
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetSumEntropyOutput() const
 {
	return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(21) );
 }








template< typename THistogram >
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementType
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetEnergy() const
 {
	return this->GetEnergyOutput()->Get();
 }

template< typename THistogram >
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementType
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetEntropy() const
 {
	return this->GetEntropyOutput()->Get();
 }

template< typename THistogram >
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementType
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetCorrelation() const
 {
	return this->GetCorrelationOutput()->Get();
 }

template< typename THistogram >
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementType
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetInverseDifferenceMoment() const
 {
	return this->GetInverseDifferenceMomentOutput()->Get();
 }

template< typename THistogram >
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementType
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetInertia() const
 {
	return this->GetInertiaOutput()->Get();
 }

template< typename THistogram >
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementType
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetClusterShade() const
 {
	return this->GetClusterShadeOutput()->Get();
 }

template< typename THistogram >
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementType
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetClusterProminence() const
 {
	return this->GetClusterProminenceOutput()->Get();
 }

template< typename THistogram >
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementType
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetHaralickCorrelation() const
 {
	return this->GetHaralickCorrelationOutput()->Get();
 }

/** news */

template< typename THistogram >
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementType
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetContrast() const
 {
	return this->GetContrastOutput()->Get();
 }


template< typename THistogram >
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementType
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetAutocorrelation() const
 {
	return this->GetAutocorrelationOutput()->Get();
 }


template< typename THistogram >
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementType
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetClusterTendency() const
 {
	return this->GetClusterTendencyOutput()->Get();
 }

template< typename THistogram >
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementType
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetDissimilarity() const
 {
	return this->GetDissimilarityOutput()->Get();
 }


template< typename THistogram >
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementType
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetHomogeneity() const
 {
	return this->GetHomogeneityOutput()->Get();
 }

template< typename THistogram >
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementType
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetHomogeneity2() const
 {
	return this->GetHomogeneity2Output()->Get();
 }


template< typename THistogram >
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementType
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetInverseDifferenceMoment2() const
 {
	return this->GetInverseDifferenceMoment2Output()->Get();
 }

template< typename THistogram >
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementType
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetDifferenceEntropy() const
 {
	return this->GetDifferenceEntropyOutput()->Get();
 }

template< typename THistogram >
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementType
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetInverseVariance() const
 {
	return this->GetInverseVarianceOutput()->Get();
 }

template< typename THistogram >
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementType
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetInformationalMeasureOfCorrelation1() const
 {
	return this->GetInformationalMeasureOfCorrelation1Output()->Get();
 }

template< typename THistogram >
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementType
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetInformationalMeasureOfCorrelation2() const
 {
	return this->GetInformationalMeasureOfCorrelation2Output()->Get();
 }

template< typename THistogram >
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementType
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetMaxProbability() const
 {
	return this->GetMaxProbabilityOutput()->Get();
 }


template< typename THistogram >
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementType
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetSumAverage() const
 {
	return this->GetSumAverageOutput()->Get();
 }

template< typename THistogram >
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementType
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetSumEntropy() const
 {
	return this->GetSumEntropyOutput()->Get();
 }





template< typename THistogram >
typename HistogramToTextureFeaturesFilterExtended< THistogram >::MeasurementType
HistogramToTextureFeaturesFilterExtended< THistogram >
::GetFeature(TextureFeatureName feature)
 {
	switch ( feature )
	{
	case Energy :
		return this->GetEnergy();
	case Entropy:
		return this->GetEntropy();
	case Correlation:
		return this->GetCorrelation();
	case InverseDifferenceMoment:
		return this->GetInverseDifferenceMoment();
	case Inertia:
		return this->GetInertia();
	case ClusterShade:
		return this->GetClusterShade();
	case ClusterProminence:
		return this->GetClusterProminence();
	case HaralickCorrelation:
		return this->GetHaralickCorrelation();
	case Contrast:
		return this->GetContrast();
	case Autocorrelation:
		return this->GetAutocorrelation();
	case ClusterTendency:
		return this->GetClusterTendency();
	case Dissimilarity:
		return this->GetDissimilarity();
	case Homogeneity:
		return this->GetHomogeneity();
	case Homogeneity2:
		return this->GetHomogeneity2();
	case InverseDifferenceMoment2:
		return this->GetInverseDifferenceMoment2();
	case DifferenceEntropy:
		return this->GetDifferenceEntropy();
	case InverseVariance:
		return this->GetInverseVariance();
	case InformationalMeasureOfCorrelation1:
		return this->GetInformationalMeasureOfCorrelation1();
	case InformationalMeasureOfCorrelation2:
		return this->GetInformationalMeasureOfCorrelation2();
		case MaxProbability:
		return this->GetMaxProbability();
		case SumAverage:
		return this->GetSumAverage();
		case SumEntropy:
		return this->GetSumEntropy();
	default:
		return 0;
	}
 }


template< typename THistogram >
void
HistogramToTextureFeaturesFilterExtended< THistogram >
::PrintSelf(std::ostream & os, Indent indent) const
 {
	Superclass::PrintSelf(os, indent);
 }
} // end of namespace Statistics
} // end of namespace itk

#endif
