#include <itkBSplineScatteredDataPointSetToImageFilter.h>
#include <itkPointSet.h>
#include <itkImageRegionIterator.h>

template <typename TTransform>
typename TTransform::Pointer resampleBSplineTransform(
    TTransform* originalTransform,
    TTransform* newTransform
) {
    using DisplacementFieldType = typename TTransform::ImageType;
    using PointSetType = itk::PointSet<typename TTransform::ParametersValueType, DisplacementFieldType::ImageDimension>;
    using BSplineFilterType = itk::BSplineScatteredDataPointSetToImageFilter<PointSetType, DisplacementFieldType>;

    // Create a point set from the displacement field of the original transform
    typename PointSetType::Pointer pointSet = PointSetType::New();
    typename DisplacementFieldType::Pointer displacementField = originalTransform->GetCoefficientImages();
    itk::ImageRegionIterator<DisplacementFieldType> it(displacementField, displacementField->GetLargestPossibleRegion());
    for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
        typename PointSetType::PointType point;
        displacementField->TransformIndexToPhysicalPoint(it.GetIndex(), point);
        pointSet->SetPoint(it.GetIndex(), point);
        pointSet->SetPointData(it.GetIndex(), it.Get());
    }

    // Fit a B-spline to the point set
    typename BSplineFilterType::Pointer bsplineFilter = BSplineFilterType::New();
    bsplineFilter->SetInput(pointSet);
    bsplineFilter->SetSplineOrder(originalTransform->GetSplineOrder());
    bsplineFilter->SetNumberOfControlPoints(newTransform->GetNumberOfParameters());  // Adjust this to your needs
    bsplineFilter->SetSize(newTransform->GetTransformDomainMeshSize());
    bsplineFilter->SetOrigin(newTransform->GetTransformDomainOrigin());
    bsplineFilter->SetSpacing(newTransform->GetTransformDomainPhysicalDimensions());
    bsplineFilter->SetDirection(newTransform->GetTransformDomainDirection());
    bsplineFilter->Update();

    // The output of the filter is a B-spline that approximates the displacement field
    newTransform->SetCoefficientImages(bsplineFilter->GetOutput());

    return newTransform;
}