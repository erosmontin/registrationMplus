#include <itkImage.h>
#include <itkCommand.h>
#include <itkLBFGSBOptimizer.h>
#include <chrono>    // << add this
#include <iomanip>                       // << for std::setprecision

#include <iostream>

class LBFGSBOptimizeCommandIterationUpdate : public itk::Command
{
public:
    using Clock = std::chrono::steady_clock;    // << add this alias
    typedef  LBFGSBOptimizeCommandIterationUpdate   Self;
    typedef  itk::Command             Superclass;
    typedef itk::SmartPointer<Self>   Pointer;
    itkNewMacro( Self );
protected:
    LBFGSBOptimizeCommandIterationUpdate()
        : m_StartTime( Clock::now() )              // << initialize start
    {};
private:
    Clock::time_point m_StartTime;               // << store start
public:
    typedef itk::LBFGSBOptimizer    OptimizerType;
    typedef   const OptimizerType * OptimizerPointer;
    void Execute(itk::Object *caller, const itk::EventObject & event)
    {
        Execute( (const itk::Object *)caller, event);
    }
    void Execute(const itk::Object * object, const itk::EventObject & event) override
    {
        OptimizerPointer optimizer = static_cast< OptimizerPointer >( object );
        if( !(itk::IterationEvent().CheckEvent( &event )) )
        {
            return;
        }

        // compute elapsed
        auto now     = Clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - m_StartTime);
        double secs  = elapsed.count() / 1000.0;

        // print with elapsed time
        std::cout << "\r[" << std::fixed << std::setprecision(2)
                  << secs << std::defaultfloat << " s] Iter: "
                  << optimizer->GetCurrentIteration()
                  << "  Metric: " << optimizer->GetCachedValue()
                  << "  Inf-norm proj-grad: "
                  << optimizer->GetInfinityNormOfProjectedGradient()
                  << std::flush;
    }
};

#include <itkRegularStepGradientDescentOptimizer.h>
class RegularStepGradientDescentOptimizerCommandIterationUpdate : public itk::Command
{
public:
    using Clock = std::chrono::steady_clock;    // << add this alias
    typedef  RegularStepGradientDescentOptimizerCommandIterationUpdate Self;
    typedef  itk::Command                                             Superclass;
    typedef itk::SmartPointer<Self>                                   Pointer;
    itkNewMacro( Self );

    itkSetMacro(ShowGradient, bool);
    itkGetMacro(ShowGradient, bool);

protected:
    RegularStepGradientDescentOptimizerCommandIterationUpdate()
      : m_ShowGradient(false),
        m_StartTime( Clock::now() )               // << initialize start clock
    {};

private:
    bool                        m_ShowGradient;
    Clock::time_point           m_StartTime;       // << store start time

public:
    typedef itk::RegularStepGradientDescentOptimizer    OptimizerType;
    typedef const OptimizerType *                       OptimizerPointer;

    void Execute(itk::Object *caller, const itk::EventObject & event) override
    {
        this->Execute( static_cast<const itk::Object*>(caller), event );
    }

    void Execute(const itk::Object * object, const itk::EventObject & event) override
    {
        OptimizerPointer optimizer = static_cast< OptimizerPointer >( object );
        if( ! itk::IterationEvent().CheckEvent( &event ) )
        {
            return;
        }

        // compute elapsed
        auto now     = Clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - m_StartTime);
        double secs  = elapsed.count() / 1000.0;

        long iteration = optimizer->GetCurrentIteration();

        // overwrite the same line
        std::cout << "\r[" << std::fixed << std::setprecision(2)
                  << secs << std::defaultfloat << " s] Iter: " << iteration
                  << "  Value: " << optimizer->GetValue();

        if (m_ShowGradient)
        {
            std::cout << "  Grad: " << optimizer->GetGradient();
        }

        std::cout << "  GradTol: " << optimizer->GetGradientMagnitudeTolerance()
                  << "  Step: "    << optimizer->GetCurrentStepLength()
                  << std::flush;
    }
};

#include "itkParticleSwarmOptimizer.h"

class ParticleSwarmOptimizeCommandIterationUpdate : public itk::Command
{
public:
    typedef  ParticleSwarmOptimizeCommandIterationUpdate   Self;
    typedef  itk::Command             Superclass;
    typedef itk::SmartPointer<Self>   Pointer;
    itkNewMacro( Self );
protected:
    ParticleSwarmOptimizeCommandIterationUpdate() {};
public:
    typedef itk::ParticleSwarmOptimizer    OptimizerType;
    typedef   const OptimizerType * OptimizerPointer;
    void Execute(itk::Object *caller, const itk::EventObject & event)
    {
        Execute( (const itk::Object *)caller, event);
    }
    void Execute(const itk::Object * object, const itk::EventObject & event)
    {
        OptimizerPointer optimizer = static_cast< OptimizerPointer >( object );
        if( !(itk::IterationEvent().CheckEvent( &event )) )
        {
            return;
        }
    if(optimizer == nullptr)
    {
        std::cout << "Optimizer is null!" << std::endl;
    }
    else
    {
        try
        {
            std::cout << optimizer->GetCurrentPosition() << "   ";
            std::cout << optimizer->GetValue() << std::endl;
        }
        catch(const std::exception& e)
        {
            std::cout << "Caught exception: " << e.what() << std::endl;
        }
    }
    }
};



template <typename TRegistration>
class RegistrationInterfaceCommand : public itk::Command
{
public:
  using Self = RegistrationInterfaceCommand;
  using Superclass = itk::Command;
  using Pointer = itk::SmartPointer<Self>;
  itkNewMacro(Self);
 
protected:
  RegistrationInterfaceCommand() = default;
 
public:
  using RegistrationType = TRegistration;
  using RegistrationPointer = RegistrationType *;
  using OptimizerType = itk::RegularStepGradientDescentOptimizer;
  using OptimizerPointer = OptimizerType *;
  void
  Execute(itk::Object * object, const itk::EventObject & event) override
  {
    if (!(itk::IterationEvent().CheckEvent(&event)))
    {
      return;
    }
    auto registration = static_cast<RegistrationPointer>(object);
    auto optimizer =
      static_cast<OptimizerPointer>(registration->GetModifiableOptimizer());
 
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "MultiResolution Level : " << registration->GetCurrentLevel()
              << std::endl;
    std::cout << std::endl;
 
    if (registration->GetCurrentLevel() == 0)
    {
        // Level 0: keep user-configured step lengths as-is (no hardcoded override)
    }
    else
    {
      optimizer->SetMaximumStepLength(optimizer->GetMaximumStepLength() /
                                      4.0);
      optimizer->SetMinimumStepLength(optimizer->GetMinimumStepLength() /
                                      10.0);
    }
        if (optimizer->GetCurrentIteration() == 0)
    {
        std::cout << "-------------------------------------" << std::endl;
        std::cout << "MultiResolution Level : " << registration->GetCurrentLevel()
                  << std::endl;

        // Print the resolution at the current level
        auto fixedImagePyramid = registration->GetFixedImagePyramid();
        auto movingImagePyramid = registration->GetMovingImagePyramid();
        std::cout << "Fixed Image Resolution : " 
                  << fixedImagePyramid->GetOutput(registration->GetCurrentLevel())->GetSpacing() 
                  << std::endl;
        std::cout << "Moving Image Resolution : " 
                  << movingImagePyramid->GetOutput(registration->GetCurrentLevel())->GetSpacing() 
                  << std::endl;

        std::cout << std::endl;
    }
};

  void
  Execute(const itk::Object *, const itk::EventObject &) override
  {
    return;
  }
};


#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBSplineTransform.h"
#include "itkResampleImageFilter.h"
#include "itkTransformToDisplacementFieldFilter.h"

// NOTE: CreateBSplineFromTransform is intentionally not implemented.
// The approach of directly copying displacement-field voxel values into B-spline
// control-point parameters is conceptually incorrect: B-spline parameters live in
// a coarse control-point grid, not in the image voxel grid, and the mapping
// requires a proper B-spline decomposition (itkBSplineDecompositionImageFilter),
// not a pixel-wise copy.  The return type declared below was also wrong (TImage::Pointer
// instead of TTransform::Pointer).  If this conversion is needed, implement it
// using itkBSplineDecompositionImageFilter on each component of the displacement field.
//
// template <typename TImage, typename TTransform>
// typename TTransform::Pointer CreateBSplineFromTransform(TTransform* transform, TImage* image,
//                                                         typename TImage::SizeType gridSize);


#include "itkTransformToDeformationFieldSource.h"



template<typename TTransformType,typename TMovingImageType,typename TDeformationFieldType>
typename TDeformationFieldType::Pointer TransformToDeformationField(typename TTransformType::Pointer transform, 
                                                                             typename TMovingImageType::Pointer movingImage)
{
    typedef itk::TransformToDeformationFieldSource<TDeformationFieldType, typename TDeformationFieldType::PixelType::ValueType> TransformToDeformationFieldSourceType;
    typename TransformToDeformationFieldSourceType::Pointer td = TransformToDeformationFieldSourceType::New();
    td->SetOutputParametersFromImage(movingImage);
    td->SetTransform(transform);
    td->Update();

    return td->GetOutput();
}

#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"

template<typename TTransformType>
void WriteTransform(const std::string& fileName, typename TTransformType::Pointer transform)
{
    #if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 5) || ITK_VERSION_MAJOR > 4
        typedef itk::TransformFileWriterTemplate<typename TTransformType::ParametersValueType> WriterType;
    #else
        typedef itk::TransformFileWriter WriterType;
    #endif

    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput(transform);
    writer->SetFileName(fileName);
    writer->Update();
}

#include "itkTransformFileReader.h"

itk::TransformBase::Pointer ReadTransformGeneric(const std::string& fileName)
{
    #if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 5) || ITK_VERSION_MAJOR > 4
        typedef itk::TransformFileReaderTemplate<double> ReaderType;
    #else
        typedef itk::TransformFileReader ReaderType;
    #endif

    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(fileName);
    reader->Update();

    return reader->GetTransformList()->front();
}

#include "itkTransformFileReader.h"

template<typename TTransformType>
typename TTransformType::Pointer ReadTransform(const std::string& fileName)
{
    #if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 5) || ITK_VERSION_MAJOR > 4
        typedef itk::TransformFileReaderTemplate<typename TTransformType::ParametersValueType> ReaderType;
    #else
        typedef itk::TransformFileReader ReaderType;
    #endif

    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(fileName);
    reader->Update();

    return dynamic_cast<TTransformType*>(reader->GetTransformList()->front().GetPointer());
}

// ── Per-iteration label-map Dice monitoring observer ──────────────────────────
// Usage:
//   auto obs = LabelMapDiceObserver<TransformType, LabelImageType>::New();
//   obs->SetFixedLabelMap(fixedLabel);
//   obs->SetMovingLabelMap(movingLabel);
//   obs->SetTransform(transform);          // concrete transform (already set by optimizer)
//   obs->SetEvaluateEveryNIterations(5);   // optional – default 1
//   obs->SetStride(4);                     // sub-sample every 4th voxel per dim (default 4)
//   optimizer->AddObserver(itk::IterationEvent(), obs);
//
// At each triggered iteration the observer appends
//   "| Dice: L1=0.85 L2=0.92 mean=0.89"
// to the current console line.

#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include <map>
#include <set>
#include <cmath>

template <typename TTransform, typename TLabelImage>
class LabelMapDiceObserver : public itk::Command
{
public:
    using Self       = LabelMapDiceObserver;
    using Superclass = itk::Command;
    using Pointer    = itk::SmartPointer<Self>;
    itkNewMacro(Self);

    using LabelPixelType  = typename TLabelImage::PixelType;
    using LabelConstPtr   = typename TLabelImage::ConstPointer;
    using TransformPtr    = typename TTransform::Pointer;
    using NNInterpType    = itk::NearestNeighborInterpolateImageFunction<TLabelImage, double>;
    using DiceMap         = std::map<LabelPixelType, double>;

    void SetFixedLabelMap(LabelConstPtr img)  { m_FixedLabelMap  = img; }
    void SetMovingLabelMap(LabelConstPtr img) { m_MovingLabelMap = img; }
    void SetTransform(TransformPtr t)         { m_Transform = t; }

    /** Only evaluate every N optimizer iterations (default = 1). */
    void SetEvaluateEveryNIterations(unsigned int n) { m_Every = (n < 1 ? 1 : n); }

    /** Sub-sampling stride per dimension (default = 4).
     *  stride=4 on a 256³ image → ~4096 samples ≈ fast. */
    void SetStride(unsigned int s) { m_Stride = (s < 1 ? 1 : s); }

    const DiceMap & GetLastDice() const { return m_LastDice; }

protected:
    LabelMapDiceObserver()
        : m_Every(1), m_Stride(4), m_IterCount(0) {}

private:
    LabelConstPtr  m_FixedLabelMap;
    LabelConstPtr  m_MovingLabelMap;
    TransformPtr   m_Transform;
    unsigned int   m_Every;
    unsigned int   m_Stride;
    unsigned long  m_IterCount;
    DiceMap        m_LastDice;

public:
    void Execute(itk::Object * caller, const itk::EventObject & event) override
    {
        this->Execute(static_cast<const itk::Object *>(caller), event);
    }

    void Execute(const itk::Object *, const itk::EventObject & event) override
    {
        if (!itk::IterationEvent().CheckEvent(&event)) return;
        if (!m_FixedLabelMap || !m_MovingLabelMap || !m_Transform) return;

        ++m_IterCount;
        if (m_IterCount % m_Every != 0) return;

        // Set up NN interpolator on the moving label map
        auto nnInterp = NNInterpType::New();
        nnInterp->SetInputImage(m_MovingLabelMap);

        // Collect unique non-zero labels from the fixed map
        std::set<LabelPixelType> labels;
        {
            itk::ImageRegionConstIterator<TLabelImage> it(
                m_FixedLabelMap, m_FixedLabelMap->GetLargestPossibleRegion());
            for (; !it.IsAtEnd(); ++it)
                if (it.Get() != 0) labels.insert(it.Get());
        }

        // Per-label Dice accumulation
        std::map<LabelPixelType, long long> countF, countM, countI;
        for (auto L : labels) { countF[L]=0; countM[L]=0; countI[L]=0; }

        unsigned int pixIdx = 0;
        itk::ImageRegionConstIteratorWithIndex<TLabelImage> it(
            m_FixedLabelMap, m_FixedLabelMap->GetLargestPossibleRegion());

        for (it.GoToBegin(); !it.IsAtEnd(); ++it, ++pixIdx)
        {
            if (pixIdx % m_Stride != 0) continue;

            const LabelPixelType fLabel = it.Get();

            // Physical point in fixed space
            typename TLabelImage::PointType fixedPt;
            m_FixedLabelMap->TransformIndexToPhysicalPoint(it.GetIndex(), fixedPt);

            // Transform to moving space
            const auto movingPt = m_Transform->TransformPoint(fixedPt);

            LabelPixelType mLabel = 0;
            if (nnInterp->IsInsideBuffer(movingPt))
                mLabel = static_cast<LabelPixelType>(
                    std::round(nnInterp->Evaluate(movingPt)));

            for (auto L : labels)
            {
                bool inF = (fLabel  == L);
                bool inM = (mLabel  == L);
                if (inF) ++countF[L];
                if (inM) ++countM[L];
                if (inF && inM) ++countI[L];
            }
        }

        // Compute and cache Dice per label
        m_LastDice.clear();
        double meanDice = 0.0;
        for (auto L : labels)
        {
            double dice = 0.0;
            if (countF[L] + countM[L] > 0)
                dice = 2.0 * countI[L] / static_cast<double>(countF[L] + countM[L]);
            m_LastDice[L] = dice;
            meanDice += dice;
        }
        if (!labels.empty()) meanDice /= static_cast<double>(labels.size());

        // Print inline with the existing metric output
        std::cout << " | Dice:";
        for (auto & kv : m_LastDice)
            std::cout << " L" << static_cast<int>(kv.first)
                      << "=" << std::fixed << std::setprecision(3) << kv.second;
        std::cout << " mean=" << meanDice << std::defaultfloat << std::flush;
    }
};

// ── Per-iteration snapshot observer ───────────────────────────────────────────
// Saves a representative mid-slice PNG (3-panel: fixed | resampled | checkerboard)
// or the full resampled 3D volume (.nii.gz) at every N iterations.
//
// Usage:
//   auto snap = IterationSnapshotObserver<TransformType, ImageType>::New();
//   snap->SetFixedImage(fixedImage);
//   snap->SetMovingImage(movingImage);
//   snap->SetTransform(transform);
//   snap->SetOutputDirectory("/path/to/snapshots");
//   snap->SetSaveEveryNIterations(5);
//   snap->SetSaveStack(false);  // true → full 3D .nii.gz per iteration
//   optimizer->AddObserver(itk::IterationEvent(), snap);

#include "itkCheckerBoardImageFilter.h"
#include "itkTileImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkRGBPixel.h"
#include "itkComposeImageFilter.h"
#include "itksys/SystemTools.hxx"

template <typename TTransform, typename TImage>
class IterationSnapshotObserver : public itk::Command
{
public:
    using Self       = IterationSnapshotObserver;
    using Superclass = itk::Command;
    using Pointer    = itk::SmartPointer<Self>;
    itkNewMacro(Self);

    static constexpr unsigned int Dim = TImage::ImageDimension;
    using PixelType      = typename TImage::PixelType;
    using SliceType      = itk::Image<PixelType, Dim - 1>;
    using UCharPixelType = unsigned char;
    using UCharSliceType = itk::Image<UCharPixelType, Dim - 1>;
    using RGBPixelType   = itk::RGBPixel<UCharPixelType>;
    using RGBSliceType   = itk::Image<RGBPixelType, Dim - 1>;

    void SetFixedImage(const TImage* img)              { m_FixedImage  = img; }
    void SetMovingImage(const TImage* img)             { m_MovingImage = img; }
    void SetTransform(typename TTransform::Pointer t)  { m_Transform   = t; }
    void SetOutputDirectory(const std::string& dir)    { m_OutputDir   = dir; }
    void SetSaveEveryNIterations(unsigned int n)       { m_Every = std::max(1u, n); }
    /** If true, save full 3D resampled volume (.nii.gz).
     *  If false (default), save a 2×2 panel PNG. */
    void SetSaveStack(bool b)                          { m_SaveStack = b; }
    /** If true, add a warped‑grid overlay panel (green lines on anatomical
     *  image) to visualise B‑spline deformation.  Default: true. */
    void SetShowDeformationGrid(bool b)                { m_ShowGrid = b; }
    /** Spacing of the regular grid lines (in voxels).  Default: 20. */
    void SetGridSpacingPixels(unsigned int s)           { m_GridSpacing = std::max(2u, s); }
    /** If true, draw the actual B-spline control-point lattice (warped by
     *  the current transform parameters) instead of a regular pixel grid.
     *  Only has effect when TTransform is itk::BSplineTransform<double,3,3>.
     *  Default: false. */
    void SetShowBSplineMesh(bool b)                    { m_ShowBSplineMesh = b; }

protected:
    IterationSnapshotObserver()
        : m_Every(1), m_SaveStack(false), m_ShowGrid(true),
          m_GridSpacing(20), m_ShowBSplineMesh(false), m_IterCount(0) {}

private:
    typename TImage::ConstPointer      m_FixedImage;
    typename TImage::ConstPointer      m_MovingImage;
    typename TTransform::Pointer       m_Transform;
    std::string                        m_OutputDir;
    unsigned int                       m_Every;
    bool                               m_SaveStack;
    bool                               m_ShowGrid;
    unsigned int                       m_GridSpacing;
    bool                               m_ShowBSplineMesh;
    unsigned long                      m_IterCount;

    // ── helpers ──────────────────────────────────────────────────────────

    typename SliceType::Pointer
    ExtractAxialSlice(const TImage* vol, unsigned int sliceIdx) const
    {
        using ExtractType = itk::ExtractImageFilter<TImage, SliceType>;
        auto ext = ExtractType::New();
        ext->SetDirectionCollapseToSubmatrix();

        auto region = vol->GetLargestPossibleRegion();
        auto sz     = region.GetSize();
        auto idx    = region.GetIndex();
        sz[Dim - 1]  = 0;
        idx[Dim - 1] = sliceIdx;
        typename TImage::RegionType sliceRegion;
        sliceRegion.SetSize(sz);
        sliceRegion.SetIndex(idx);

        ext->SetExtractionRegion(sliceRegion);
        ext->SetInput(vol);
        ext->Update();
        typename SliceType::Pointer out = ext->GetOutput();
        out->DisconnectPipeline();
        return out;
    }

    typename UCharSliceType::Pointer
    ToUChar(const SliceType* slice) const
    {
        using R = itk::RescaleIntensityImageFilter<SliceType, UCharSliceType>;
        auto r = R::New();
        r->SetInput(slice);
        r->SetOutputMinimum(0);
        r->SetOutputMaximum(255);
        r->Update();
        typename UCharSliceType::Pointer out = r->GetOutput();
        out->DisconnectPipeline();
        return out;
    }

    /** Convert a grayscale UChar slice to an RGB slice (all channels equal). */
    typename RGBSliceType::Pointer
    GrayToRGB(const UCharSliceType* gray) const
    {
        auto rgb = RGBSliceType::New();
        rgb->CopyInformation(gray);
        rgb->SetRegions(gray->GetLargestPossibleRegion());
        rgb->Allocate();

        itk::ImageRegionConstIterator<UCharSliceType> gIt(gray, gray->GetLargestPossibleRegion());
        itk::ImageRegionIterator<RGBSliceType>        rIt(rgb,  rgb->GetLargestPossibleRegion());
        for (gIt.GoToBegin(), rIt.GoToBegin(); !gIt.IsAtEnd(); ++gIt, ++rIt)
        {
            RGBPixelType px;
            px.SetRed(gIt.Get());
            px.SetGreen(gIt.Get());
            px.SetBlue(gIt.Get());
            rIt.Set(px);
        }
        rgb->DisconnectPipeline();
        return rgb;
    }

    /** Create a padded RGB panel: the input image is placed at (offX, offY)
     *  on a black canvas of size newW × newH. */
    typename RGBSliceType::Pointer
    PadRGBPanel(const RGBSliceType* input, int newW, int newH,
                int offX, int offY) const
    {
        auto padded = RGBSliceType::New();
        typename RGBSliceType::IndexType startIdx; startIdx.Fill(0);
        typename RGBSliceType::SizeType  padSize;
        padSize[0] = static_cast<unsigned int>(newW);
        padSize[1] = static_cast<unsigned int>(newH);
        typename RGBSliceType::RegionType padRegion;
        padRegion.SetIndex(startIdx);
        padRegion.SetSize(padSize);
        padded->SetRegions(padRegion);
        padded->Allocate();
        RGBPixelType black; black.SetRed(0); black.SetGreen(0); black.SetBlue(0);
        padded->FillBuffer(black);

        itk::ImageRegionConstIterator<RGBSliceType> it(
            input, input->GetLargestPossibleRegion());
        for (it.GoToBegin(); !it.IsAtEnd(); ++it)
        {
            auto si = it.GetIndex();
            typename RGBSliceType::IndexType di;
            di[0] = si[0] + offX;
            di[1] = si[1] + offY;
            if (di[0] >= 0 && di[0] < newW && di[1] >= 0 && di[1] < newH)
                padded->SetPixel(di, it.Get());
        }
        return padded;
    }

    /** Compute the 2D pixel bounding box of the B-spline control lattice
     *  projected onto the fixed image (considering all z-slices).
     *  Returns {minX, minY, maxX, maxY} in fixed-image pixel coords. */
    struct MeshBBox { double minX, minY, maxX, maxY; };
    MeshBBox ComputeMeshBBox2D(const TImage* fixedImg) const
    {
        using BSTransformType = itk::BSplineTransform<double, 3, 3>;
        auto imgSz = fixedImg->GetLargestPossibleRegion().GetSize();
        MeshBBox bb = {0.0, 0.0,
                       static_cast<double>(imgSz[0] - 1),
                       static_cast<double>(imgSz[1] - 1)};
        auto* bst = dynamic_cast<BSTransformType*>(m_Transform.GetPointer());
        if (!bst) return bb;

        using CoeffImageType = typename BSTransformType::ImageType;
        auto coeffImages = bst->GetCoefficientImages();
        auto coeffImg = coeffImages[0];
        auto coeffSize = coeffImg->GetLargestPossibleRegion().GetSize();
        const unsigned int nx = coeffSize[0], ny = coeffSize[1], nz = coeffSize[2];

        for (unsigned int k = 0; k < nz; ++k)
        for (unsigned int j = 0; j < ny; ++j)
        for (unsigned int i = 0; i < nx; ++i)
        {
            typename CoeffImageType::IndexType idx3;
            idx3[0] = i; idx3[1] = j; idx3[2] = k;
            itk::Point<double, 3> physPt;
            coeffImg->TransformIndexToPhysicalPoint(idx3, physPt);
            itk::Point<double, 3> displaced;
            for (unsigned d = 0; d < 3; ++d)
                displaced[d] = physPt[d] + coeffImages[d]->GetPixel(idx3);
            itk::ContinuousIndex<double, 3> ci;
            fixedImg->TransformPhysicalPointToContinuousIndex(displaced, ci);
            bb.minX = std::min(bb.minX, ci[0]);
            bb.minY = std::min(bb.minY, ci[1]);
            bb.maxX = std::max(bb.maxX, ci[0]);
            bb.maxY = std::max(bb.maxY, ci[1]);
        }
        return bb;
    }

    /** Draw a thin 1-pixel border rectangle on the canvas to mark the
     *  original image boundary (drawn in dim yellow). */
    void DrawImageBorder(typename RGBSliceType::Pointer& canvas,
                         int canvasW, int canvasH,
                         int offX, int offY, int imgW, int imgH) const
    {
        // top edge
        DrawLineRGB(canvas.GetPointer(), canvasW, canvasH,
                    offX, offY, offX + imgW - 1, offY, 100, 100, 40);
        // bottom edge
        DrawLineRGB(canvas.GetPointer(), canvasW, canvasH,
                    offX, offY + imgH - 1, offX + imgW - 1, offY + imgH - 1, 100, 100, 40);
        // left edge
        DrawLineRGB(canvas.GetPointer(), canvasW, canvasH,
                    offX, offY, offX, offY + imgH - 1, 100, 100, 40);
        // right edge
        DrawLineRGB(canvas.GetPointer(), canvasW, canvasH,
                    offX + imgW - 1, offY, offX + imgW - 1, offY + imgH - 1, 100, 100, 40);
    }

    /** Create a 3D grid image in *moving-image* space.
     *  Lines every m_GridSpacing voxels; on-line pixels = 1, rest = 0. */
    typename TImage::Pointer
    MakeGridImage(const TImage* ref) const
    {
        auto grid = TImage::New();
        grid->CopyInformation(ref);
        grid->SetRegions(ref->GetLargestPossibleRegion());
        grid->Allocate();
        grid->FillBuffer(0);

        using IteratorType = itk::ImageRegionIteratorWithIndex<TImage>;
        for (IteratorType it(grid, grid->GetLargestPossibleRegion()); !it.IsAtEnd(); ++it)
        {
            auto idx = it.GetIndex();
            bool onLine = false;
            for (unsigned d = 0; d < Dim; ++d)
            {
                if (idx[d] % static_cast<long>(m_GridSpacing) == 0)
                    onLine = true;
            }
            if (onLine)
                it.Set(static_cast<PixelType>(1));
        }
        return grid;
    }

    /** Overlay green grid lines on an anatomical RGB slice.
     *  Where gridMask > 0.5, the pixel becomes bright green;
     *  elsewhere it keeps its original value. */
    typename RGBSliceType::Pointer
    OverlayGreenGrid(const RGBSliceType* anatomy,
                     const UCharSliceType* gridMask) const
    {
        auto out = RGBSliceType::New();
        out->CopyInformation(anatomy);
        out->SetRegions(anatomy->GetLargestPossibleRegion());
        out->Allocate();

        itk::ImageRegionConstIterator<RGBSliceType>   aIt(anatomy,  anatomy->GetLargestPossibleRegion());
        itk::ImageRegionConstIterator<UCharSliceType>  gIt(gridMask, gridMask->GetLargestPossibleRegion());
        itk::ImageRegionIterator<RGBSliceType>         oIt(out,      out->GetLargestPossibleRegion());

        for (aIt.GoToBegin(), gIt.GoToBegin(), oIt.GoToBegin();
             !oIt.IsAtEnd(); ++aIt, ++gIt, ++oIt)
        {
            if (gIt.Get() > 128)
            {
                RGBPixelType px;
                px.SetRed(0);
                px.SetGreen(255);
                px.SetBlue(0);
                oIt.Set(px);
            }
            else
            {
                oIt.Set(aIt.Get());
            }
        }
        out->DisconnectPipeline();
        return out;
    }

    /** Bresenham line draw in green (or any RGB) directly on an RGBSliceType. */
    static void DrawLineRGB(RGBSliceType* img, int W, int H,
                            int x0, int y0, int x1, int y1,
                            unsigned char r, unsigned char g, unsigned char b)
    {
        int dx =  std::abs(x1 - x0), sx = (x0 < x1) ? 1 : -1;
        int dy = -std::abs(y1 - y0), sy = (y0 < y1) ? 1 : -1;
        int err = dx + dy;
        while (true)
        {
            if (x0 >= 0 && x0 < W && y0 >= 0 && y0 < H)
            {
                typename RGBSliceType::IndexType idx;
                idx[0] = x0;  idx[1] = y0;
                RGBPixelType px;  px.SetRed(r);  px.SetGreen(g);  px.SetBlue(b);
                img->SetPixel(idx, px);
            }
            if (x0 == x1 && y0 == y1) break;
            int e2 = 2 * err;
            if (e2 >= dy) { err += dy;  x0 += sx; }
            if (e2 <= dx) { err += dx;  y0 += sy; }
        }
    }

    /** Overlay the actual B-spline control-point lattice on an RGB slice.
     *
     *  Enumerates all control-point nodes, maps each through the current
     *  transform (TransformPoint), projects to 2-D pixel coords on
     *  fixedImg, then draws X- and Y-direction edges for any edge whose
     *  midpoint lies within 'sliceTol' slices of sliceZ.
     *
     *  Falls back silently if TTransform is not
     *  itk::BSplineTransform<double,3,3>.
     */
    void OverlayBSplineMesh(typename RGBSliceType::Pointer& rgbInOut,
                            const TImage* fixedImg,
                            unsigned int sliceZ,
                            int offsetX = 0, int offsetY = 0) const
    {
        using BSTransformType = itk::BSplineTransform<double, 3, 3>;
        auto* bst = dynamic_cast<BSTransformType*>(m_Transform.GetPointer());
        if (!bst) return;

        // Use the coefficient images directly for exact node positions
        // and raw displacement coefficients (no B-spline interpolation).
        using CoeffImageType = typename BSTransformType::ImageType;
        auto coeffImages = bst->GetCoefficientImages();
        auto coeffImg = coeffImages[0];  // reference grid
        auto coeffSize = coeffImg->GetLargestPossibleRegion().GetSize();

        const unsigned int nx = static_cast<unsigned int>(coeffSize[0]);
        const unsigned int ny = static_cast<unsigned int>(coeffSize[1]);
        const unsigned int nz = static_cast<unsigned int>(coeffSize[2]);

        // Use the *canvas* dimensions (which may be padded)
        auto sz2D = rgbInOut->GetLargestPossibleRegion().GetSize();
        const int W = static_cast<int>(sz2D[0]);
        const int H = static_cast<int>(sz2D[1]);

        // Tolerance: half a node spacing in z (in voxel units)
        const double nodeSpacingZ = coeffImg->GetSpacing()[2];
        const double zSpacingVox  = nodeSpacingZ / fixedImg->GetSpacing()[2];
        const double sliceTol     = std::max(1.5, zSpacingVox * 0.6);

        // Precompute displaced positions projected to fixed-image 2D.
        // Displaced position = node physical position + raw coefficient.
        struct Proj { double px, py, pz; };
        std::vector<Proj> proj(nx * ny * nz);
        for (unsigned int k = 0; k < nz; ++k)
        for (unsigned int j = 0; j < ny; ++j)
        for (unsigned int i = 0; i < nx; ++i)
        {
            typename CoeffImageType::IndexType idx3;
            idx3[0] = i;  idx3[1] = j;  idx3[2] = k;

            // Physical position of this coefficient node
            itk::Point<double, 3> physPt;
            coeffImg->TransformIndexToPhysicalPoint(idx3, physPt);

            // Add the raw displacement coefficient
            itk::Point<double, 3> displaced;
            for (unsigned d = 0; d < 3; ++d)
                displaced[d] = physPt[d] + coeffImages[d]->GetPixel(idx3);

            // Project to fixed image continuous index, then add canvas offset
            itk::ContinuousIndex<double, 3> ci;
            fixedImg->TransformPhysicalPointToContinuousIndex(displaced, ci);
            proj[(k*ny + j)*nx + i] = { ci[0] + offsetX,
                                        ci[1] + offsetY,
                                        ci[2] };
        }

        const double sliceZd = static_cast<double>(sliceZ);

        // Draw edge between two projected nodes if their midpoint Z is
        // near sliceZ.  No 2D clipping — the canvas is padded to fit.
        auto drawEdge = [&](const Proj& a, const Proj& b)
        {
            const double midZ = (a.pz + b.pz) * 0.5;
            if (std::abs(midZ - sliceZd) > sliceTol) return;
            DrawLineRGB(rgbInOut.GetPointer(), W, H,
                        static_cast<int>(std::round(a.px)),
                        static_cast<int>(std::round(a.py)),
                        static_cast<int>(std::round(b.px)),
                        static_cast<int>(std::round(b.py)),
                        0, 255, 0);
        };

        for (unsigned int k = 0; k < nz; ++k)
        for (unsigned int j = 0; j < ny; ++j)
        for (unsigned int i = 0; i < nx; ++i)
        {
            const Proj& cur = proj[(k*ny + j)*nx + i];
            if (i + 1 < nx)
                drawEdge(cur, proj[(k*ny + j)*nx + i + 1]);
            if (j + 1 < ny)
                drawEdge(cur, proj[(k*ny + (j+1))*nx + i]);
        }
    }

public:
    void Execute(itk::Object* caller, const itk::EventObject& ev) override
    { this->Execute(static_cast<const itk::Object*>(caller), ev); }

    void Execute(const itk::Object*, const itk::EventObject& ev) override
    {
        if (!itk::IterationEvent().CheckEvent(&ev))  return;
        if (!m_FixedImage || !m_MovingImage || !m_Transform) return;
        if (m_OutputDir.empty()) return;

        ++m_IterCount;
        if (m_IterCount % m_Every != 0) return;

        try
        {
            itksys::SystemTools::MakeDirectory(m_OutputDir);

            // Resample the moving image with the current transform
            using Resample = itk::ResampleImageFilter<TImage, TImage>;
            auto rs = Resample::New();
            rs->SetInput(m_MovingImage);
            rs->SetTransform(m_Transform);
            rs->SetSize(m_FixedImage->GetLargestPossibleRegion().GetSize());
            rs->SetOutputSpacing(m_FixedImage->GetSpacing());
            rs->SetOutputOrigin(m_FixedImage->GetOrigin());
            rs->SetOutputDirection(m_FixedImage->GetDirection());
            rs->SetDefaultPixelValue(0);
            rs->Update();
            typename TImage::Pointer resampled = rs->GetOutput();

            if (m_SaveStack)
            {
                std::ostringstream fn;
                fn << m_OutputDir << "/iter_"
                   << std::setfill('0') << std::setw(4) << m_IterCount
                   << ".nii.gz";

                using W = itk::ImageFileWriter<TImage>;
                auto w = W::New();
                w->SetFileName(fn.str());
                w->SetInput(resampled);
                w->SetUseCompression(true);
                w->Update();
            }
            else
            {
                // ── 2×2 panel PNG ───────────────────────────────────────
                //  (1,1) fixed               | (1,2) registered moving
                //  (2,1) checkerboard         | (2,2) registered + green grid
                auto sz = m_FixedImage->GetLargestPossibleRegion().GetSize();
                unsigned int midZ = sz[Dim - 1] / 2;

                auto fixSlice = ExtractAxialSlice(m_FixedImage.GetPointer(), midZ);
                auto movSlice = ExtractAxialSlice(resampled.GetPointer(),    midZ);

                // Checkerboard
                using CB = itk::CheckerBoardImageFilter<SliceType>;
                auto cb = CB::New();
                cb->SetInput1(fixSlice);
                cb->SetInput2(movSlice);
                typename CB::PatternArrayType pat;
                pat.Fill(8);
                cb->SetCheckerPattern(pat);
                cb->Update();
                typename SliceType::Pointer chkSlice = cb->GetOutput();
                chkSlice->DisconnectPipeline();

                // Convert to unsigned-char [0,255]
                auto fUC  = ToUChar(fixSlice.GetPointer());
                auto mUC  = ToUChar(movSlice.GetPointer());
                auto cUC  = ToUChar(chkSlice.GetPointer());

                // Convert all 3 grayscale panels to RGB
                auto fRGB = GrayToRGB(fUC.GetPointer());
                auto mRGB = GrayToRGB(mUC.GetPointer());
                auto cRGB = GrayToRGB(cUC.GetPointer());

                // ── Compute padding to show full B-spline domain ────────
                int padL = 0, padR = 0, padT = 0, padB = 0;
                auto imgSz2 = m_FixedImage->GetLargestPossibleRegion().GetSize();
                const int imgW = static_cast<int>(imgSz2[0]);
                const int imgH = static_cast<int>(imgSz2[1]);

                if (m_ShowBSplineMesh)
                {
                    auto bb = ComputeMeshBBox2D(m_FixedImage.GetPointer());
                    padL = std::max(0, static_cast<int>(std::ceil(-bb.minX)) + 3);
                    padT = std::max(0, static_cast<int>(std::ceil(-bb.minY)) + 3);
                    padR = std::max(0, static_cast<int>(std::ceil(bb.maxX - (imgW - 1))) + 3);
                    padB = std::max(0, static_cast<int>(std::ceil(bb.maxY - (imgH - 1))) + 3);

                    // Print diagnostic once
                    if (m_IterCount == m_Every)
                    {
                        std::cerr << "[Snapshot] Image size: " << imgW << "x" << imgH
                                  << ", Mesh bbox (px): ["
                                  << bb.minX << ", " << bb.minY << "] -> ["
                                  << bb.maxX << ", " << bb.maxY << "]"
                                  << ", Padding L/R/T/B: "
                                  << padL << "/" << padR << "/" << padT << "/" << padB
                                  << std::endl;
                    }
                }

                bool needPad = (padL > 0 || padR > 0 || padT > 0 || padB > 0);
                int canvasW = imgW + padL + padR;
                int canvasH = imgH + padT + padB;

                // If padding needed, pad all 4 panels to same size
                typename RGBSliceType::Pointer fP, mP, cP;
                if (needPad)
                {
                    fP = PadRGBPanel(fRGB.GetPointer(), canvasW, canvasH, padL, padT);
                    mP = PadRGBPanel(mRGB.GetPointer(), canvasW, canvasH, padL, padT);
                    cP = PadRGBPanel(cRGB.GetPointer(), canvasW, canvasH, padL, padT);
                }
                else
                {
                    fP = fRGB;
                    mP = mRGB;
                    cP = cRGB;
                }

                // Panel (2,2): registered moving with overlay
                //   m_ShowBSplineMesh → real B-spline control-point lattice
                //   m_ShowGrid        → regular pixel-spaced warped grid
                //   neither           → plain registered moving image
                typename RGBSliceType::Pointer gridRGB;
                if (m_ShowBSplineMesh)
                {
                    // Create padded canvas with registered moving image,
                    // then overlay the full B-spline lattice.
                    gridRGB = PadRGBPanel(mRGB.GetPointer(), canvasW, canvasH, padL, padT);
                    OverlayBSplineMesh(gridRGB, m_FixedImage.GetPointer(), midZ, padL, padT);
                    // Draw dim border showing original image extent
                    DrawImageBorder(gridRGB, canvasW, canvasH, padL, padT, imgW, imgH);
                }
                else if (m_ShowGrid)
                {
                    // Build a regular grid in moving-image space, then
                    // resample it into fixed space with the same transform
                    // used for the moving image.  Grid lines that were
                    // originally straight will bend wherever the B-spline
                    // deforms.
                    auto gridImg = MakeGridImage(m_MovingImage.GetPointer());

                    using ResampleGrid = itk::ResampleImageFilter<TImage, TImage>;
                    auto rsg = ResampleGrid::New();
                    rsg->SetInput(gridImg);
                    rsg->SetTransform(m_Transform);
                    rsg->SetSize(m_FixedImage->GetLargestPossibleRegion().GetSize());
                    rsg->SetOutputSpacing(m_FixedImage->GetSpacing());
                    rsg->SetOutputOrigin(m_FixedImage->GetOrigin());
                    rsg->SetOutputDirection(m_FixedImage->GetDirection());
                    rsg->SetDefaultPixelValue(0);
                    rsg->Update();

                    auto warpedGridSlice = ExtractAxialSlice(rsg->GetOutput(), midZ);
                    auto wgUC = ToUChar(warpedGridSlice.GetPointer());

                    gridRGB = OverlayGreenGrid(mRGB.GetPointer(), wgUC.GetPointer());
                }
                else
                {
                    gridRGB = mRGB;
                }

                // Tile 2×2 as RGB:
                //  (0) top-left=fixed  (1) top-right=registered
                //  (2) bot-left=checker (3) bot-right=registered+grid
                using Tile = itk::TileImageFilter<RGBSliceType, RGBSliceType>;
                auto tiler = Tile::New();
                tiler->SetInput(0, fP);         // (1,1) fixed
                tiler->SetInput(1, mP);         // (1,2) registered moving
                tiler->SetInput(2, cP);         // (2,1) checkerboard
                tiler->SetInput(3, gridRGB);    // (2,2) registered + green grid
                itk::FixedArray<unsigned int, Dim - 1> layout;
                layout[0] = 2;   // 2 columns
                layout[1] = 2;   // 2 rows
                tiler->SetLayout(layout);
                tiler->Update();

                // Write PNG
                std::ostringstream fn;
                fn << m_OutputDir << "/iter_"
                   << std::setfill('0') << std::setw(4) << m_IterCount
                   << ".png";

                using W = itk::ImageFileWriter<RGBSliceType>;
                auto w = W::New();
                w->SetFileName(fn.str());
                w->SetInput(tiler->GetOutput());
                w->Update();
            }
        }
        catch (const itk::ExceptionObject& ex)
        {
            std::cerr << "\n[SnapshotObserver] Warning: " << ex.what() << std::endl;
        }
        catch (const std::exception& ex)
        {
            std::cerr << "\n[SnapshotObserver] Warning: " << ex.what() << std::endl;
        }
    }
};