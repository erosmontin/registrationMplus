
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
