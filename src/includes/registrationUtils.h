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

    void SetFixedImage(const TImage* img)              { m_FixedImage  = img; }
    void SetMovingImage(const TImage* img)             { m_MovingImage = img; }
    void SetTransform(typename TTransform::Pointer t)  { m_Transform   = t; }
    void SetOutputDirectory(const std::string& dir)    { m_OutputDir   = dir; }
    void SetSaveEveryNIterations(unsigned int n)       { m_Every = std::max(1u, n); }
    /** If true, save full 3D resampled volume (.nii.gz).
     *  If false (default), save a mid-axial 3-panel PNG. */
    void SetSaveStack(bool b)                          { m_SaveStack = b; }

protected:
    IterationSnapshotObserver()
        : m_Every(1), m_SaveStack(false), m_IterCount(0) {}

private:
    typename TImage::ConstPointer      m_FixedImage;
    typename TImage::ConstPointer      m_MovingImage;
    typename TTransform::Pointer       m_Transform;
    std::string                        m_OutputDir;
    unsigned int                       m_Every;
    bool                               m_SaveStack;
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
        sz[Dim - 1]  = 0;           // collapse last (axial) dimension
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
            // Ensure output directory exists
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
                // ── full 3D volume ──────────────────────────────────────
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
                // ── mid-axial 3-panel PNG ───────────────────────────────
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
                pat.Fill(8);          // 8×8 checkerboard grid
                cb->SetCheckerPattern(pat);
                cb->Update();
                typename SliceType::Pointer chkSlice = cb->GetOutput();
                chkSlice->DisconnectPipeline();

                // Convert all three to unsigned-char [0,255]
                auto fUC  = ToUChar(fixSlice.GetPointer());
                auto mUC  = ToUChar(movSlice.GetPointer());
                auto cUC  = ToUChar(chkSlice.GetPointer());

                // Tile horizontally: [fixed | resampled | checkerboard]
                using Tile = itk::TileImageFilter<UCharSliceType, UCharSliceType>;
                auto tiler = Tile::New();
                tiler->SetInput(0, fUC);
                tiler->SetInput(1, mUC);
                tiler->SetInput(2, cUC);
                itk::FixedArray<unsigned int, Dim - 1> layout;
                layout[0] = 3;   // 3 columns
                layout[1] = 1;   // 1 row
                tiler->SetLayout(layout);
                tiler->Update();

                // Write PNG
                std::ostringstream fn;
                fn << m_OutputDir << "/iter_"
                   << std::setfill('0') << std::setw(4) << m_IterCount
                   << ".png";

                using W = itk::ImageFileWriter<UCharSliceType>;
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