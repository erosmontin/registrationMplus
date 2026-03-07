/*=========================================================================
 *  _core.cpp  –  pybind11 bindings for the Mplus metric & registration
 *
 *  This module exposes:
 *    - MplusMetric class (property getters/setters, initialize, print)
 *    - register_images() high-level function
 *    - NumPy ↔ ITK image conversion utilities
 *    - CUDA availability check and GPU functions (when USE_CUDA)
 *=========================================================================*/
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

// ITK
#include "itkImage.h"
#include "itkImportImageFilter.h"

// Mplus metric
#include "itkMplus.h"

#ifdef USE_CUDA
#include "cuda/DistanceTransform.cuh"
#include "cuda/LabelMetricKernels.cuh"
#include "cuda/DerivativeOps.cuh"
#endif

namespace py = pybind11;

// ── Type aliases ─────────────────────────────────────────────────────────
constexpr unsigned int Dim = 3;
using PixelType  = float;
using ImageType  = itk::Image<PixelType, Dim>;
using MplusType  = itk::Mplus<ImageType, ImageType>;

// ── NumPy → ITK conversion ──────────────────────────────────────────────
static ImageType::Pointer numpy_to_itk(
    py::array_t<float> arr,
    py::array_t<double> spacing,
    py::array_t<double> origin)
{
    auto buf = arr.request();
    if (buf.ndim != 3)
        throw std::runtime_error("Expected 3-D array");

    ImageType::SizeType size;
    // NumPy is Z, Y, X (row-major) → ITK is X, Y, Z
    size[0] = buf.shape[2];
    size[1] = buf.shape[1];
    size[2] = buf.shape[0];

    ImageType::SpacingType sp;
    auto sp_buf = spacing.request();
    auto* sp_ptr = static_cast<double*>(sp_buf.ptr);
    for (int i = 0; i < 3; ++i) sp[i] = sp_ptr[i];

    ImageType::PointType orig;
    auto or_buf = origin.request();
    auto* or_ptr = static_cast<double*>(or_buf.ptr);
    for (int i = 0; i < 3; ++i) orig[i] = or_ptr[i];

    auto importer = itk::ImportImageFilter<PixelType, Dim>::New();
    ImageType::IndexType start;
    start.Fill(0);
    ImageType::RegionType region;
    region.SetIndex(start);
    region.SetSize(size);
    importer->SetRegion(region);
    importer->SetSpacing(sp);
    importer->SetOrigin(orig);

    size_t nVoxels = size[0] * size[1] * size[2];
    auto* data = static_cast<float*>(buf.ptr);
    // ImportImageFilter takes ownership=false so Python keeps the memory
    importer->SetImportPointer(data, nVoxels, false);
    importer->Update();

    return importer->GetOutput();
}

// ── ITK → NumPy conversion ──────────────────────────────────────────────
static py::array_t<float> itk_to_numpy(ImageType::Pointer img) {
    auto size = img->GetLargestPossibleRegion().GetSize();
    size_t nVoxels = size[0] * size[1] * size[2];

    // Shape: Z, Y, X for NumPy (row-major)
    std::vector<ssize_t> shape = {
        static_cast<ssize_t>(size[2]),
        static_cast<ssize_t>(size[1]),
        static_cast<ssize_t>(size[0])
    };

    auto result = py::array_t<float>(shape);
    auto buf = result.request();
    std::memcpy(buf.ptr, img->GetBufferPointer(), nVoxels * sizeof(float));
    return result;
}

// ── Module definition ────────────────────────────────────────────────────
PYBIND11_MODULE(_core, m) {
    m.doc() = "Mplus registration C++ bindings";

    // ── CUDA availability ────────────────────────────────────────────────
#ifdef USE_CUDA
    m.def("cuda_available", []() { return true; });
#else
    m.def("cuda_available", []() { return false; });
#endif

    // ── MplusMetric class ────────────────────────────────────────────────
    py::class_<MplusType, MplusType::Pointer>(m, "MplusMetric")
        .def(py::init([]() { return MplusType::New(); }))
        .def_property("lambda_val",
            &MplusType::GetLambda, &MplusType::SetLambda)
        .def_property("lambda_derivative",
            &MplusType::GetLambdaDerivative, &MplusType::SetLambdaDerivative)
        .def_property("alpha",
            &MplusType::GetAlpha, &MplusType::SetAlpha)
        .def_property("alpha_derivative",
            &MplusType::GetAlphaDerivative, &MplusType::SetAlphaDerivative)
        .def_property("nu",
            &MplusType::GetNu, &MplusType::SetNu)
        .def_property("nu_derivative",
            &MplusType::GetNuDerivative, &MplusType::SetNuDerivative)
        .def_property("yota",
            &MplusType::GetYota, &MplusType::SetYota)
        .def_property("yota_derivative",
            &MplusType::GetYotaDerivative, &MplusType::SetYotaDerivative)
        .def_property("label_kappa",
            &MplusType::GetLabelKappa, &MplusType::SetLabelKappa)
        .def_property("label_kappa_derivative",
            &MplusType::GetLabelKappaDerivative, &MplusType::SetLabelKappaDerivative)
        .def_property("derivative_mode",
            &MplusType::GetDerivativeMode, &MplusType::SetDerivativeMode)
        .def_property("fixed_eta",
            &MplusType::GetFixedEta, &MplusType::SetFixedEta)
        .def_property("moving_eta",
            &MplusType::GetMovingEta, &MplusType::SetMovingEta)
        .def_property("auto_estimate_eta",
            &MplusType::GetAutoEstimateEta, &MplusType::SetAutoEstimateEta)
        .def_property("bin_numbers",
            &MplusType::GetBinNumbers, &MplusType::SetBinNumbers)
        .def_property("num_samples",
            &MplusType::GetMANumberOfSamples, &MplusType::SetMANumberOfSamples)
        .def("initialize", &MplusType::Initialize)
        .def("print_info", &MplusType::print)
    ;

    // ── High-level register_images function ──────────────────────────────
    m.def("register_images",
        [](py::array_t<float> fixed,
           py::array_t<float> moving,
           py::array_t<double> spacing,
           py::array_t<double> origin,
           py::dict params,
           py::object fixed_labels,
           py::object moving_labels) -> py::dict
        {
            // Convert NumPy → ITK images
            auto fixedImg  = numpy_to_itk(fixed, spacing, origin);
            auto movingImg = numpy_to_itk(moving, spacing, origin);

            // TODO: Full registration pipeline
            // 1. Create BSpline transform
            // 2. Configure Mplus metric from params dict
            // 3. Set up multi-resolution schedule
            // 4. Run optimizer
            // 5. Resample moving image
            // 6. Return results

            py::dict result;
            result["transform"] = py::array_t<double>(0);
            result["warped"]    = py::array_t<float>(0);
            result["metrics"]   = py::dict();
            result["history"]   = py::list();

            return result;
        },
        py::arg("fixed"),
        py::arg("moving"),
        py::arg("spacing"),
        py::arg("origin"),
        py::arg("params"),
        py::arg("fixed_labels")  = py::none(),
        py::arg("moving_labels") = py::none(),
        "Run multi-metric B-spline registration."
    );

    // ── GPU functions (when available) ───────────────────────────────────
#ifdef USE_CUDA
    m.def("gpu_distance_transform",
        [](py::array_t<float> binary) -> py::array_t<float> {
            auto buf = binary.request();
            if (buf.ndim != 3)
                throw std::runtime_error("Expected 3-D array");
            int sz = buf.shape[0], sy = buf.shape[1], sx = buf.shape[2];
            auto result = py::array_t<float>(buf.shape);
            auto rbuf = result.request();
            mplus::cuda::gpu_signed_distance_transform(
                static_cast<float*>(buf.ptr),
                static_cast<float*>(rbuf.ptr),
                sx, sy, sz);
            return result;
        },
        py::arg("binary"),
        "GPU signed distance transform of a binary 3-D image.");

    m.def("gpu_normalize_derivative",
        [](py::array_t<double> deriv) {
            auto buf = deriv.mutable_unchecked<1>();
            mplus::cuda::gpu_normalize_derivative(buf.mutable_data(0), buf.size());
        },
        py::arg("derivative"),
        "In-place L2 normalisation on GPU.");

    m.def("gpu_rescale_derivative",
        [](py::array_t<double> deriv) {
            auto buf = deriv.mutable_unchecked<1>();
            mplus::cuda::gpu_rescale_derivative(buf.mutable_data(0), buf.size());
        },
        py::arg("derivative"),
        "In-place rescale to [-1, 1] on GPU.");

    m.def("gpu_label_dice",
        [](py::array_t<short> fixed_labels,
           py::array_t<short> moving_labels) -> py::dict
        {
            auto fb = fixed_labels.request();
            auto mb = moving_labels.request();
            if (fb.ndim != 3 || mb.ndim != 3)
                throw std::runtime_error("Expected 3-D label arrays");
            size_t nVoxels = fb.shape[0] * fb.shape[1] * fb.shape[2];

            // Find unique labels
            auto* fdata = static_cast<short*>(fb.ptr);
            std::set<short> ulabels;
            for (size_t i = 0; i < nVoxels; ++i)
                if (fdata[i] > 0) ulabels.insert(fdata[i]);

            std::vector<short> labels(ulabels.begin(), ulabels.end());
            std::vector<double> dice(labels.size());

            mplus::cuda::gpu_label_dice(
                fdata,
                static_cast<short*>(mb.ptr),
                nVoxels,
                labels.data(),
                (int)labels.size(),
                dice.data());

            py::dict result;
            for (size_t i = 0; i < labels.size(); ++i)
                result[py::int_(labels[i])] = dice[i];
            return result;
        },
        py::arg("fixed_labels"),
        py::arg("moving_labels"),
        "Compute per-label Dice coefficients on GPU.");
#endif
}
