#include <itkImage.h>
#include <itkCommand.h>
#include <itkLBFGSBOptimizer.h>
#include <chrono>    // << add this
#include <iomanip>                       // << for std::setprecision
#include <algorithm>
#include <vector>

#include <iostream>
#include <functional>
#include <fstream>

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
#include "itkResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"
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
    /** Overlay line width in pixels.  Values <= 0 enable adaptive sizing. */
    void SetOverlayLineWidthPixels(double w)           { m_LineWidthPixels = (w > 0.0) ? w : 0.0; }
    /** If true, draw the actual B-spline control-point lattice (warped by
     *  the current transform parameters) instead of a regular pixel grid.
     *  Only has effect when TTransform is itk::BSplineTransform<double,3,3>.
     *  Default: false. */
    void SetShowBSplineMesh(bool b)                    { m_ShowBSplineMesh = b; }

    /** Provide a callback that returns a name→value map of per-sub-metric
     *  values for the current iteration.  Called at every snapshot; the map
     *  is written to a CSV file (metrics.csv) inside the output directory. */
    void SetMetricValuesGetter(std::function<std::map<std::string,double>()> fn)
    { m_MetricGetter = fn; }

    /** Called once after optimization completes.  Writes convergence.png
     *  into the snapshot directory using accumulated history. */
    void FinalizeConvergencePlot() const
    {
        if (m_OutputDir.empty() || m_IterHistory.empty() || m_SeriesNames.empty()) return;

        const int W = 800, H = 400;
        const int padL = 60, padR = 20, padT = 30, padB = 50;
        const int plotW = W - padL - padR;
        const int plotH = H - padT - padB;

        // White canvas
        auto canvas = RGBSliceType::New();
        typename RGBSliceType::IndexType si0; si0.Fill(0);
        typename RGBSliceType::SizeType  sz0; sz0[0] = W; sz0[1] = H;
        typename RGBSliceType::RegionType reg0; reg0.SetIndex(si0); reg0.SetSize(sz0);
        canvas->SetRegions(reg0);
        canvas->Allocate();
        RGBPixelType white; white.SetRed(255); white.SetGreen(255); white.SetBlue(255);
        canvas->FillBuffer(white);

        // Axes
        DrawLineRGB(canvas.GetPointer(), W, H, padL, padT, padL, padT+plotH, 0, 0, 0);
        DrawLineRGB(canvas.GetPointer(), W, H, padL, padT+plotH, padL+plotW, padT+plotH, 0, 0, 0);

        // Find global y range
        double vMin =  std::numeric_limits<double>::max();
        double vMax = -std::numeric_limits<double>::max();
        for (const auto& sv : m_SeriesValues)
            for (double v : sv) { if (v < vMin) vMin = v; if (v > vMax) vMax = v; }
        if (vMax <= vMin) vMax = vMin + 1.0;

        const double iterRange = std::max(1.0,
            static_cast<double>(m_IterHistory.back() - m_IterHistory.front()));
        const double vRange = vMax - vMin;

        // Color palette per series (C++14-safe, no structured bindings)
        struct RGB3 { unsigned char r, g, b; };
        const RGB3 palette[] = {
            {220,  20,  20},  // red     (Total)
            {  0, 180,   0},  // green   (MI)
            {  0, 100, 255},  // blue    (NGF)
            {255, 140,   0},  // orange  (MSE)
            {170,   0, 220},  // purple  (NC)
            {  0, 200, 200},  // cyan    (GD)
            {200,  20, 200},  // magenta (NMI)
            {100, 180,  30},  // lime    (Label)
        };
        const int nPal = static_cast<int>(sizeof(palette) / sizeof(palette[0]));

        for (size_t s = 0; s < m_SeriesValues.size(); ++s)
        {
            const auto& S = m_SeriesValues[s];
            if (S.size() < 2) continue;
            const RGB3& c = palette[static_cast<int>(s) % nPal];

            // Draw line series
            for (size_t j = 1; j < S.size(); ++j)
            {
                int x0 = padL + static_cast<int>(
                    (m_IterHistory[j-1] - m_IterHistory[0]) / iterRange * plotW);
                int y0 = padT + plotH - static_cast<int>(
                    (S[j-1] - vMin) / vRange * plotH);
                int x1 = padL + static_cast<int>(
                    (m_IterHistory[j] - m_IterHistory[0]) / iterRange * plotW);
                int y1 = padT + plotH - static_cast<int>(
                    (S[j] - vMin) / vRange * plotH);
                DrawLineRGB(canvas.GetPointer(), W, H,
                            x0, y0, x1, y1, c.r, c.g, c.b,
                            2.2, 0.95);
            }

            // Legend swatch + text label
            const int swX = padL + plotW - 130;
            const int swY = padT + 8 + static_cast<int>(s) * 14;
            DrawLineRGB(canvas.GetPointer(), W, H,
                        swX, swY + 1, swX + 22, swY + 1, c.r, c.g, c.b,
                        3.0, 1.0);
            // Draw series name next to swatch
            if (s < m_SeriesNames.size())
                DrawText(canvas.GetPointer(), W, H,
                         swX + 26, swY - 2, m_SeriesNames[s], c.r, c.g, c.b);
        }

        // Axis tick marks (5 ticks each) with value labels
        for (int t = 0; t <= 4; ++t)
        {
            // Y-axis ticks
            int y = padT + plotH - t * plotH / 4;
            DrawLineRGB(canvas.GetPointer(), W, H, padL-5, y, padL, y, 0, 0, 0);
            // Y-axis tick value
            {
                double val = vMin + t * vRange / 4.0;
                std::ostringstream oss;
                if (std::abs(val) < 0.01 && val != 0.0)
                    oss << std::scientific << std::setprecision(1) << val;
                else
                    oss << std::fixed << std::setprecision(2) << val;
                std::string label = oss.str();
                int labelW = static_cast<int>(label.size()) * 6;
                DrawText(canvas.GetPointer(), W, H,
                         padL - 6 - labelW, y - 3, label, 0, 0, 0);
            }

            // X-axis ticks
            int x = padL + t * plotW / 4;
            DrawLineRGB(canvas.GetPointer(), W, H, x, padT+plotH, x, padT+plotH+5, 0, 0, 0);
            // X-axis tick value (iteration number)
            {
                unsigned long iterVal = m_IterHistory.front()
                    + static_cast<unsigned long>(t * iterRange / 4.0);
                std::ostringstream oss;
                oss << iterVal;
                std::string label = oss.str();
                int labelW = static_cast<int>(label.size()) * 6;
                DrawText(canvas.GetPointer(), W, H,
                         x - labelW / 2, padT + plotH + 8, label, 0, 0, 0);
            }
        }

        // Axis titles
        DrawText(canvas.GetPointer(), W, H,
                 padL + plotW / 2 - 30, H - 12, "ITERATION", 0, 0, 0);
        DrawText(canvas.GetPointer(), W, H,
                 2, padT - 12, "METRIC", 0, 0, 0);

        const std::string plotPath = m_OutputDir + "/convergence.png";
        using PW = itk::ImageFileWriter<RGBSliceType>;
        auto pw = PW::New();
        pw->SetFileName(plotPath);
        pw->SetInput(canvas);
        try
        {
            pw->Update();
            std::cout << "[Snapshot] Convergence plot saved: " << plotPath << std::endl;
        }
        catch (const std::exception& ex)
        {
            std::cerr << "[Snapshot] Warning: could not save convergence plot: "
                      << ex.what() << std::endl;
        }
    }

protected:
    IterationSnapshotObserver()
        : m_Every(1), m_SaveStack(false), m_ShowGrid(true),
          m_GridSpacing(20), m_LineWidthPixels(0.0),
          m_ShowBSplineMesh(false), m_IterCount(0),
          m_CSVHeaderWritten(false) {}

private:
    typename TImage::ConstPointer      m_FixedImage;
    typename TImage::ConstPointer      m_MovingImage;
    typename TTransform::Pointer       m_Transform;
    std::string                        m_OutputDir;
    unsigned int                       m_Every;
    bool                               m_SaveStack;
    bool                               m_ShowGrid;
    unsigned int                       m_GridSpacing;
    double                             m_LineWidthPixels;
    bool                               m_ShowBSplineMesh;
    unsigned long                      m_IterCount;

    // ── per-metric value logging ──────────────────────────────────────────
    std::function<std::map<std::string,double>()> m_MetricGetter;
    std::vector<unsigned long>                    m_IterHistory;
    std::vector<std::string>                      m_SeriesNames;
    std::vector<std::vector<double>>              m_SeriesValues;  // [series][snapshot]
    bool                                          m_CSVHeaderWritten;

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

    /** Rescale a slice to [0,255] using percentile-based windowing.
     *  Clips at [loPercentile, hiPercentile] to handle outliers, then maps
     *  each image independently — preserving good contrast regardless of the
     *  absolute intensity range (e.g. different scanners or modalities). */
    typename UCharSliceType::Pointer
    ToUCharPercentile(const SliceType* slice,
                      double loPercentile = 1.0,
                      double hiPercentile = 99.0) const
    {
        // Collect all pixel values into a sorted vector for percentile computation
        std::vector<PixelType> vals;
        vals.reserve(slice->GetLargestPossibleRegion().GetNumberOfPixels());
        itk::ImageRegionConstIterator<SliceType> it(slice, slice->GetLargestPossibleRegion());
        for (it.GoToBegin(); !it.IsAtEnd(); ++it)
            vals.push_back(it.Get());

        std::sort(vals.begin(), vals.end());

        const size_t N = vals.size();
        const PixelType wMin = vals[static_cast<size_t>(loPercentile / 100.0 * (N - 1))];
        const PixelType wMax = vals[static_cast<size_t>(hiPercentile / 100.0 * (N - 1))];

        using W = itk::IntensityWindowingImageFilter<SliceType, UCharSliceType>;
        auto w = W::New();
        w->SetInput(slice);
        w->SetWindowMinimum(wMin == wMax ? wMin - 1 : wMin);
        w->SetWindowMaximum(wMin == wMax ? wMax + 1 : wMax);
        w->SetOutputMinimum(0);
        w->SetOutputMaximum(255);
        w->Update();
        typename UCharSliceType::Pointer out = w->GetOutput();
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
        bool haveNode = false;

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

            if (!haveNode)
            {
                bb = { ci[0], ci[1], ci[0], ci[1] };
                haveNode = true;
                continue;
            }

            bb.minX = std::min(bb.minX, ci[0]);
            bb.minY = std::min(bb.minY, ci[1]);
            bb.maxX = std::max(bb.maxX, ci[0]);
            bb.maxY = std::max(bb.maxY, ci[1]);
        }

        if (!haveNode)
        {
            bb = {0.0, 0.0,
                  static_cast<double>(imgSz[0] - 1),
                  static_cast<double>(imgSz[1] - 1)};
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
     *  Grid intensity is used as alpha so resampled lines stay smooth. */
    typename RGBSliceType::Pointer
    OverlayGreenGrid(const RGBSliceType* anatomy,
                     const UCharSliceType* gridMask) const
    {
        auto out = RGBSliceType::New();
        out->CopyInformation(anatomy);
        out->SetRegions(anatomy->GetLargestPossibleRegion());
        out->Allocate();

        auto effectiveMask = ExpandMaskForLineWidth(
            gridMask,
            ResolveOverlayLineWidthPixels(static_cast<double>(m_GridSpacing)));

        itk::ImageRegionConstIterator<RGBSliceType>   aIt(anatomy,  anatomy->GetLargestPossibleRegion());
        itk::ImageRegionConstIterator<UCharSliceType>  gIt(effectiveMask, effectiveMask->GetLargestPossibleRegion());
        itk::ImageRegionIterator<RGBSliceType>         oIt(out,      out->GetLargestPossibleRegion());

        for (aIt.GoToBegin(), gIt.GoToBegin(), oIt.GoToBegin();
             !oIt.IsAtEnd(); ++aIt, ++gIt, ++oIt)
        {
            const double alpha = 0.85 * (static_cast<double>(gIt.Get()) / 255.0);
            if (alpha > 0.0)
            {
                RGBPixelType px = aIt.Get();
                px.SetRed(static_cast<UCharPixelType>(
                    std::round(px.GetRed() * (1.0 - alpha))));
                px.SetGreen(static_cast<UCharPixelType>(
                    std::round(px.GetGreen() * (1.0 - alpha) + 255.0 * alpha)));
                px.SetBlue(static_cast<UCharPixelType>(
                    std::round(px.GetBlue() * (1.0 - alpha))));
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

    double ResolveOverlayLineWidthPixels(double nominalSpacingPixels) const
    {
        if (m_LineWidthPixels > 0.0)
            return std::max(0.8, m_LineWidthPixels);
        if (nominalSpacingPixels > 0.0)
            return std::min(3.5, std::max(1.15, nominalSpacingPixels * 0.08));
        return 1.25;
    }

    typename UCharSliceType::Pointer
    ExpandMaskForLineWidth(const UCharSliceType* mask, double widthPixels) const
    {
        auto out = UCharSliceType::New();
        out->CopyInformation(mask);
        out->SetRegions(mask->GetLargestPossibleRegion());
        out->Allocate();

        const int radius = static_cast<int>(
            std::floor(std::max(0.0, widthPixels - 1.25) * 0.5));

        const auto region = mask->GetLargestPossibleRegion();
        const auto size = region.GetSize();
        const auto index = region.GetIndex();

        if (radius <= 0)
        {
            itk::ImageRegionConstIterator<UCharSliceType> src(mask, region);
            itk::ImageRegionIterator<UCharSliceType> dst(out, region);
            for (src.GoToBegin(), dst.GoToBegin(); !dst.IsAtEnd(); ++src, ++dst)
                dst.Set(src.Get());
            out->DisconnectPipeline();
            return out;
        }

        for (int y = 0; y < static_cast<int>(size[1]); ++y)
        {
            for (int x = 0; x < static_cast<int>(size[0]); ++x)
            {
                unsigned char maxValue = 0;
                for (int oy = -radius; oy <= radius; ++oy)
                {
                    for (int ox = -radius; ox <= radius; ++ox)
                    {
                        if (ox * ox + oy * oy > radius * radius)
                            continue;

                        const int sx = x + ox;
                        const int sy = y + oy;
                        if (sx < 0 || sy < 0 ||
                            sx >= static_cast<int>(size[0]) ||
                            sy >= static_cast<int>(size[1]))
                            continue;

                        typename UCharSliceType::IndexType sampleIdx;
                        sampleIdx[0] = index[0] + sx;
                        sampleIdx[1] = index[1] + sy;
                        maxValue = std::max(maxValue, mask->GetPixel(sampleIdx));
                    }
                }

                typename UCharSliceType::IndexType outIdx;
                outIdx[0] = index[0] + x;
                outIdx[1] = index[1] + y;
                out->SetPixel(outIdx, maxValue);
            }
        }

        out->DisconnectPipeline();
        return out;
    }

    /** Resample a 2D slice to isotropic pixel spacing (smallest spacing wins).
     *  If spacing is already isotropic (within 1% tolerance), returns the
     *  input unchanged.  This fixes the squashed-looking debug PNGs that
     *  occur with anisotropic voxels (e.g. sagittal acquisitions). */
    typename SliceType::Pointer
    ResampleSliceIsotropic(typename SliceType::Pointer slice) const
    {
        auto sp = slice->GetSpacing();
        double minSp = sp[0];
        for (unsigned d = 1; d < Dim - 1; ++d)
            if (sp[d] < minSp) minSp = sp[d];

        bool needResample = false;
        for (unsigned d = 0; d < Dim - 1; ++d)
        {
            if (std::abs(sp[d] - minSp) / minSp > 0.01)
            { needResample = true; break; }
        }
        if (!needResample) return slice;

        // Compute new size to cover the same physical extent
        auto oldSize = slice->GetLargestPossibleRegion().GetSize();
        typename SliceType::SizeType    newSize;
        typename SliceType::SpacingType newSpacing;
        for (unsigned d = 0; d < Dim - 1; ++d)
        {
            newSpacing[d] = minSp;
            newSize[d] = static_cast<typename SliceType::SizeType::SizeValueType>(
                std::ceil(oldSize[d] * sp[d] / minSp));
        }

        using ResampleSlice = itk::ResampleImageFilter<SliceType, SliceType>;
        auto rs = ResampleSlice::New();
        rs->SetInput(slice);
        rs->SetSize(newSize);
        rs->SetOutputSpacing(newSpacing);
        rs->SetOutputOrigin(slice->GetOrigin());
        rs->SetOutputDirection(slice->GetDirection());
        rs->SetDefaultPixelValue(0);
        using LinInterp = itk::LinearInterpolateImageFunction<SliceType, double>;
        rs->SetInterpolator(LinInterp::New());
        rs->Update();
        typename SliceType::Pointer out = rs->GetOutput();
        out->DisconnectPipeline();
        return out;
    }

    static double ClampUnit(double value)
    {
        return std::max(0.0, std::min(1.0, value));
    }

    static void BlendPixelRGB(RGBSliceType* img, int W, int H,
                              int x, int y,
                              unsigned char r, unsigned char g, unsigned char b,
                              double alpha)
    {
        if (x < 0 || x >= W || y < 0 || y >= H || alpha <= 0.0)
            return;

        typename RGBSliceType::IndexType idx;
        idx[0] = x; idx[1] = y;
        RGBPixelType px = img->GetPixel(idx);
        const double a = ClampUnit(alpha);
        px.SetRed(static_cast<UCharPixelType>(
            std::round(px.GetRed() * (1.0 - a) + r * a)));
        px.SetGreen(static_cast<UCharPixelType>(
            std::round(px.GetGreen() * (1.0 - a) + g * a)));
        px.SetBlue(static_cast<UCharPixelType>(
            std::round(px.GetBlue() * (1.0 - a) + b * a)));
        img->SetPixel(idx, px);
    }

    static double DistancePointToSegment(double px, double py,
                                         double x0, double y0,
                                         double x1, double y1)
    {
        const double dx = x1 - x0;
        const double dy = y1 - y0;
        const double len2 = dx * dx + dy * dy;
        if (len2 <= 1e-12)
            return std::hypot(px - x0, py - y0);

        const double t = ClampUnit(((px - x0) * dx + (py - y0) * dy) / len2);
        const double projX = x0 + t * dx;
        const double projY = y0 + t * dy;
        return std::hypot(px - projX, py - projY);
    }

    /** Anti-aliased alpha-blended line draw directly on an RGBSliceType. */
    static void DrawLineRGB(RGBSliceType* img, int W, int H,
                            double x0, double y0, double x1, double y1,
                            unsigned char r, unsigned char g, unsigned char b,
                            double widthPixels = 1.0,
                            double opacity = 1.0)
    {
        const double width = std::max(1.0, widthPixels);
        const double radius = 0.5 * width;

        const int minX = static_cast<int>(std::floor(std::min(x0, x1) - radius - 1.0));
        const int maxX = static_cast<int>(std::ceil (std::max(x0, x1) + radius + 1.0));
        const int minY = static_cast<int>(std::floor(std::min(y0, y1) - radius - 1.0));
        const int maxY = static_cast<int>(std::ceil (std::max(y0, y1) + radius + 1.0));

        for (int y = minY; y <= maxY; ++y)
        {
            for (int x = minX; x <= maxX; ++x)
            {
                const double dist = DistancePointToSegment(
                    x + 0.5, y + 0.5, x0, y0, x1, y1);
                const double coverage = ClampUnit(radius + 0.5 - dist);
                BlendPixelRGB(img, W, H, x, y, r, g, b, opacity * coverage);
            }
        }
    }

    /** Minimal 5×7 bitmap font for drawing text labels on RGB images.
     *  Covers A-Z, a-z (rendered as uppercase), 0-9, '.', '-', '+', ' '.
     *  Each glyph is 5 pixels wide, 7 pixels tall. */
    static void DrawText(RGBSliceType* img, int W, int H,
                         int startX, int startY, const std::string& text,
                         unsigned char r, unsigned char g, unsigned char b)
    {
        // 5×7 font bitmaps (each row is one column of 7 bits, LSB = top)
        // Index: 0-9 → digits, 10-35 → A-Z, 36 → '.', 37 → '-', 38 → '+', 39 → ' ', 40 → ':'
        static const unsigned char font[][5] = {
            {0x3E,0x51,0x49,0x45,0x3E}, // 0
            {0x00,0x42,0x7F,0x40,0x00}, // 1
            {0x42,0x61,0x51,0x49,0x46}, // 2
            {0x21,0x41,0x45,0x4B,0x31}, // 3
            {0x18,0x14,0x12,0x7F,0x10}, // 4
            {0x27,0x45,0x45,0x45,0x39}, // 5
            {0x3C,0x4A,0x49,0x49,0x30}, // 6
            {0x01,0x71,0x09,0x05,0x03}, // 7
            {0x36,0x49,0x49,0x49,0x36}, // 8
            {0x06,0x49,0x49,0x29,0x1E}, // 9
            {0x7E,0x11,0x11,0x11,0x7E}, // A 10
            {0x7F,0x49,0x49,0x49,0x36}, // B
            {0x3E,0x41,0x41,0x41,0x22}, // C
            {0x7F,0x41,0x41,0x22,0x1C}, // D
            {0x7F,0x49,0x49,0x49,0x41}, // E
            {0x7F,0x09,0x09,0x09,0x01}, // F
            {0x3E,0x41,0x49,0x49,0x7A}, // G
            {0x7F,0x08,0x08,0x08,0x7F}, // H
            {0x00,0x41,0x7F,0x41,0x00}, // I
            {0x20,0x40,0x41,0x3F,0x01}, // J
            {0x7F,0x08,0x14,0x22,0x41}, // K
            {0x7F,0x40,0x40,0x40,0x40}, // L
            {0x7F,0x02,0x0C,0x02,0x7F}, // M
            {0x7F,0x04,0x08,0x10,0x7F}, // N
            {0x3E,0x41,0x41,0x41,0x3E}, // O
            {0x7F,0x09,0x09,0x09,0x06}, // P
            {0x3E,0x41,0x51,0x21,0x5E}, // Q
            {0x7F,0x09,0x19,0x29,0x46}, // R
            {0x46,0x49,0x49,0x49,0x31}, // S
            {0x01,0x01,0x7F,0x01,0x01}, // T
            {0x3F,0x40,0x40,0x40,0x3F}, // U
            {0x1F,0x20,0x40,0x20,0x1F}, // V
            {0x3F,0x40,0x38,0x40,0x3F}, // W
            {0x63,0x14,0x08,0x14,0x63}, // X
            {0x07,0x08,0x70,0x08,0x07}, // Y
            {0x61,0x51,0x49,0x45,0x43}, // Z 35
            {0x00,0x60,0x60,0x00,0x00}, // . 36
            {0x08,0x08,0x08,0x08,0x08}, // - 37
            {0x08,0x08,0x3E,0x08,0x08}, // + 38
            {0x00,0x00,0x00,0x00,0x00}, // ' ' 39
            {0x00,0x36,0x36,0x00,0x00}, // : 40
        };

        int cx = startX;
        for (size_t ci = 0; ci < text.size(); ++ci)
        {
            char ch = text[ci];
            int gi = -1;
            if (ch >= '0' && ch <= '9') gi = ch - '0';
            else if (ch >= 'A' && ch <= 'Z') gi = ch - 'A' + 10;
            else if (ch >= 'a' && ch <= 'z') gi = ch - 'a' + 10; // lowercase → uppercase
            else if (ch == '.') gi = 36;
            else if (ch == '-') gi = 37;
            else if (ch == '+') gi = 38;
            else if (ch == ' ') gi = 39;
            else if (ch == ':') gi = 40;
            else if (ch == 'e' || ch == 'E') gi = 14; // E
            else { cx += 6; continue; } // unknown → skip

            for (int col = 0; col < 5; ++col)
            {
                unsigned char bits = font[gi][col];
                for (int row = 0; row < 7; ++row)
                {
                    if (bits & (1 << row))
                    {
                        int px = cx + col;
                        int py = startY + row;
                        if (px >= 0 && px < W && py >= 0 && py < H)
                        {
                            typename RGBSliceType::IndexType idx;
                            idx[0] = px; idx[1] = py;
                            RGBPixelType pixel;
                            pixel.SetRed(r); pixel.SetGreen(g); pixel.SetBlue(b);
                            img->SetPixel(idx, pixel);
                        }
                    }
                }
            }
            cx += 6; // 5px glyph + 1px spacing
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
                            int offsetX = 0, int offsetY = 0,
                            double scaleX = 1.0, double scaleY = 1.0) const
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
        const double nominalSpacingPx = std::min(
            std::abs(coeffImg->GetSpacing()[0] / fixedImg->GetSpacing()[0]) * scaleX,
            std::abs(coeffImg->GetSpacing()[1] / fixedImg->GetSpacing()[1]) * scaleY);
        const double lineWidthPx = ResolveOverlayLineWidthPixels(nominalSpacingPx);

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

            // Project to fixed image continuous index, scale for isotropic
            // resampling, then add canvas offset
            itk::ContinuousIndex<double, 3> ci;
            fixedImg->TransformPhysicalPointToContinuousIndex(displaced, ci);
            proj[(k*ny + j)*nx + i] = { ci[0] * scaleX + offsetX,
                                        ci[1] * scaleY + offsetY,
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
                        a.px, a.py, b.px, b.py,
                        0, 255, 0, lineWidthPx, 0.95);
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

                auto fixSliceRaw = ExtractAxialSlice(m_FixedImage.GetPointer(), midZ);
                auto movSliceRaw = ExtractAxialSlice(resampled.GetPointer(),    midZ);

                // Resample to isotropic pixels so anisotropic voxels
                // (e.g. sagittal acquisitions) render with correct aspect ratio
                auto fixSlice = ResampleSliceIsotropic(fixSliceRaw);
                auto movSlice = ResampleSliceIsotropic(movSliceRaw);

                // Normalise each image independently using percentile windowing
                // [1st, 99th] percentile clips outliers while preserving visible contrast.
                // Each image is treated independently so different intensity ranges
                // (e.g. different scanners/modalities) both look good.
                auto fUC = ToUCharPercentile(fixSlice.GetPointer());
                auto mUC = ToUCharPercentile(movSlice.GetPointer());

                // Checkerboard is built from the already-normalised UChar slices
                // so both patches have comparable contrast in the mosaic.
                using CBu = itk::CheckerBoardImageFilter<UCharSliceType>;
                auto cbu = CBu::New();
                cbu->SetInput1(fUC);
                cbu->SetInput2(mUC);
                typename CBu::PatternArrayType patu;
                patu.Fill(8);
                cbu->SetCheckerPattern(patu);
                cbu->Update();
                typename UCharSliceType::Pointer cUC = cbu->GetOutput();
                cUC->DisconnectPipeline();

                // Convert all 3 grayscale panels to RGB
                auto fRGB = GrayToRGB(fUC.GetPointer());
                auto mRGB = GrayToRGB(mUC.GetPointer());
                auto cRGB = GrayToRGB(cUC.GetPointer());

                // ── Compute padding to show full B-spline domain ────────
                // Use the (possibly resampled-to-isotropic) slice dimensions
                int padL = 0, padR = 0, padT = 0, padB = 0;
                auto sliceSz2 = fixSlice->GetLargestPossibleRegion().GetSize();
                const int imgW = static_cast<int>(sliceSz2[0]);
                const int imgH = static_cast<int>(sliceSz2[1]);

                // Compute scale factors from original fixed-image voxel
                // coords to isotropic-resampled pixel coords
                auto fixSp3D = m_FixedImage->GetSpacing();
                double minSp2D = std::min(fixSp3D[0], fixSp3D[1]);
                const double isoScaleX = fixSp3D[0] / minSp2D;
                const double isoScaleY = fixSp3D[1] / minSp2D;

                if (m_ShowBSplineMesh)
                {
                    auto bb = ComputeMeshBBox2D(m_FixedImage.GetPointer());
                    // Scale bbox to isotropic pixel space
                    bb.minX *= isoScaleX;  bb.maxX *= isoScaleX;
                    bb.minY *= isoScaleY;  bb.maxY *= isoScaleY;
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
                    OverlayBSplineMesh(gridRGB, m_FixedImage.GetPointer(), midZ,
                                       padL, padT, isoScaleX, isoScaleY);
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

                    auto warpedGridSliceRaw = ExtractAxialSlice(rsg->GetOutput(), midZ);
                    auto warpedGridSlice = ResampleSliceIsotropic(warpedGridSliceRaw);
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

            // ── log per-sub-metric values ─────────────────────────────────────
            if (m_MetricGetter)
            {
                auto vals = m_MetricGetter();

                // First call: build series list and write CSV header
                if (!m_CSVHeaderWritten)
                {
                    for (const auto& kv : vals) m_SeriesNames.push_back(kv.first);
                    m_SeriesValues.resize(m_SeriesNames.size());
                    std::ofstream hdr(m_OutputDir + "/metrics.csv");
                    hdr << "iteration";
                    for (const auto& n : m_SeriesNames) hdr << "," << n;
                    hdr << "\n";
                    m_CSVHeaderWritten = true;
                }

                // Append row to CSV
                m_IterHistory.push_back(m_IterCount);
                {
                    std::ofstream csv(m_OutputDir + "/metrics.csv", std::ios::app);
                    csv << m_IterCount;
                    for (size_t i = 0; i < m_SeriesNames.size(); ++i)
                    {
                        double v = 0.0;
                        auto it2 = vals.find(m_SeriesNames[i]);
                        if (it2 != vals.end()) v = it2->second;
                        csv << "," << std::setprecision(10) << v;
                        if (i < m_SeriesValues.size())
                            m_SeriesValues[i].push_back(v);
                    }
                    csv << "\n";
                }

                // Print to stdout
                std::cout << "[Metrics] iter " << m_IterCount;
                for (const auto& kv : vals)
                    std::cout << "  " << kv.first << "=" << std::setprecision(6) << kv.second;
                std::cout << std::endl;
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
