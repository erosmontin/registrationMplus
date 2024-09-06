
#include "itkImage.h"
#include <iostream>

template <typename TImage>
bool checkIfSize(typename TImage::Pointer meshImage, typename TImage::Pointer fixedImage) {
    const typename TImage::SizeType& meshSize = meshImage->GetLargestPossibleRegion().GetSize();
    const typename TImage::SizeType& fixedSize = fixedImage->GetLargestPossibleRegion().GetSize();

    if (meshSize != fixedSize) {
        std::cout << "Mesh image and fixed image have different sizes." << std::endl;
        return false;
    }

    return true;
}

template <typename TImage>
bool checkIfResolution(typename TImage::Pointer meshImage, typename TImage::Pointer fixedImage) {
    const typename TImage::SpacingType& meshSpacing = meshImage->GetSpacing();
    const typename TImage::SpacingType& fixedSpacing = fixedImage->GetSpacing();

    if (meshSpacing != fixedSpacing) {
        std::cout << "Mesh image and fixed image have different resolutions." << std::endl;
        return false;
    }

    return true;
}


template <typename TImage>
bool checkIfOrigin(typename TImage::Pointer meshImage, typename TImage::Pointer fixedImage) {
    const typename TImage::PointType& meshOrigin = meshImage->GetOrigin();
    const typename TImage::PointType& fixedOrigin = fixedImage->GetOrigin();

    if (meshOrigin != fixedOrigin) {
        std::cout << "Mesh image and fixed image have different origins." << std::endl;
        return false;
    }

    return true;
}

template <typename TImage>
bool checkIfDirection(typename TImage::Pointer meshImage, typename TImage::Pointer fixedImage) {
    const typename TImage::DirectionType& meshDirection = meshImage->GetDirection();
    const typename TImage::DirectionType& fixedDirection = fixedImage->GetDirection();

    if (meshDirection != fixedDirection) {
        std::cout << "Mesh image and fixed image have different directions." << std::endl;
        return false;
    }

    return true;
}


template <typename TImage>
bool checkIfSameSpace(typename TImage::Pointer image1, typename TImage::Pointer image2) {
    if (!checkIfSize<TImage>(image1, image2)) {
        return false;
    }

    if (!checkIfResolution<TImage>(image1, image2)) {
        return false;
    }

    if (!checkIfOrigin<TImage>(image1, image2)) {
        return false;
    }

    if (!checkIfDirection<TImage>(image1, image2)) {
        return false;
    }

    return true;
}


void deb(std::string s)
{
			std::cout << "-----------------here" <<s<< std::endl;

}
#include "itkImage.h"
#include "itkRandomImageSource.h"

template <typename TImage>
typename TImage::Pointer CreateRandomImage(typename TImage::SizeType size)
{
    using RandomImageSource = itk::RandomImageSource<TImage>;
    typename RandomImageSource::Pointer randomImageSource = RandomImageSource::New();


    randomImageSource->SetSize(size);
    randomImageSource->SetMin(0.0);  // minimum intensity value
    randomImageSource->SetMax(255.0);  // maximum intensity value

    randomImageSource->Update();

    return randomImageSource->GetOutput();
}



template<typename TImage, typename TTransform>
typename TImage::Pointer ApplyTransform(typename TImage::Pointer image, typename TTransform::Pointer transform)
{
    using ResampleFilterType = itk::ResampleImageFilter<TImage, TImage>;
    typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();

    resampler->SetInput(image);
    resampler->SetSize(image->GetLargestPossibleRegion().GetSize());
    resampler->SetOutputOrigin(image->GetOrigin());
    resampler->SetOutputSpacing(image->GetSpacing());
    resampler->SetOutputDirection(image->GetDirection());
    resampler->SetTransform(transform);
    resampler->Update();

    return resampler->GetOutput();
}

template<typename TImage>
void SaveImage(typename TImage::Pointer image, const std::string& filename)
{
    using WriterType = itk::ImageFileWriter<TImage>;
    typename WriterType::Pointer writer = WriterType::New();

    writer->SetFileName(filename);
    writer->SetInput(image);
    writer->Update();
}



#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"


template<typename TDeformationFieldType>
typename TDeformationFieldType::Pointer ReadDeformationField(const std::string& fileName)
{
    typedef itk::ImageFileReader<TDeformationFieldType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(fileName);
    reader->Update();

    return reader->GetOutput();
}
template<typename TDeformationFieldType>
void WriteDeformationField(const std::string& fileName, 
                           typename TDeformationFieldType::Pointer deformationField)
{
    typedef itk::ImageFileWriter<TDeformationFieldType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(fileName);
    writer->SetInput(deformationField);
    writer->Update();
}

