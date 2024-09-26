
#include <itkLabelOverlapMeasuresImageFilter.h>

std::map<int, double> calculateDiceCoefficients(itk::Image<itk::LabelMap<2>::LabelType, 2>::Pointer labelMap1, itk::Image<itk::LabelMap<2>::LabelType, 2>::Pointer labelMap2) {
    typedef itk::Image<itk::LabelMap<2>::LabelType, 2> LabelImageType;
    typedef itk::LabelOverlapMeasuresImageFilter<LabelImageType> LabelOverlapMeasuresFilterType;

    LabelOverlapMeasuresFilterType::Pointer labelOverlapMeasuresFilter = LabelOverlapMeasuresFilterType::New();
    labelOverlapMeasuresFilter->SetSourceImage(labelMap1);
    labelOverlapMeasuresFilter->SetTargetImage(labelMap2);
    labelOverlapMeasuresFilter->Update();

    std::map<int, double> diceCoefficients;
    for (unsigned int i = 0; i < labelOverlapMeasuresFilter->GetNumberOfLabels(); i++) {
        int label = labelOverlapMeasuresFilter->GetLabelSet()->GetNthElement(i);
        double diceCoefficient = labelOverlapMeasuresFilter->GetDiceCoefficient(label);
        diceCoefficients[label] = diceCoefficient;
    }

    return diceCoefficients;
}


#include <gtest/gtest.h>
#include <itkImage.h>
#include <itkLabelMap.h>
#include <itkLabelImageToLabelMapFilter.h>

extern std::map<int, double> calculateDiceCoefficients(itk::Image<itk::LabelMap<2>::LabelType, 2>::Pointer labelMap1, itk::Image<itk::LabelMap<2>::LabelType, 2>::Pointer labelMap2);

TEST(CalculateDiceCoefficientsTest, BasicTest) {
    typedef itk::Image<itk::LabelMap<2>::LabelType, 2> LabelImageType;
    typedef itk::LabelMap<2> LabelMapType;
    typedef itk::LabelImageToLabelMapFilter<LabelImageType, LabelMapType> LabelImageToLabelMapFilterType;

    // Create and set the first label map
    LabelImageType::Pointer labelImage1 = LabelImageType::New();
    // ... code to initialize labelImage1 ...
    LabelImageToLabelMapFilterType::Pointer filter1 = LabelImageToLabelMapFilterType::New();
    filter1->SetInput(labelImage1);
    filter1->Update();
    LabelMapType::Pointer labelMap1 = filter1->GetOutput();

    // Create and set the second label map
    LabelImageType::Pointer labelImage2 = LabelImageType::New();
    // ... code to initialize labelImage2 ...
    LabelImageToLabelMapFilterType::Pointer filter2 = LabelImageToLabelMapFilterType::New();
    filter2->SetInput(labelImage2);
    filter2->Update();
    LabelMapType::Pointer labelMap2 = filter2->GetOutput();

    // Calculate the Dice coefficients
    std::map<int, double> diceCoefficients = calculateDiceCoefficients(labelMap1, labelMap2);

    // ... code to check the dice coefficients ...
    EXPECT_EQ(diceCoefficients.size(), 2);
    EXPECT_NEAR(diceCoefficients[1], 0.8, 0.01);
    EXPECT_NEAR(diceCoefficients[2], 0.6, 0.01);
    // ... more checks ...
}
 
}