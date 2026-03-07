/*=========================================================================
 *  itkMplusCompat.h  –  ITK 4 ↔ ITK 5 compatibility shims
 *
 *  Include this header instead of version-specific ITK headers when
 *  code must build against both ITK 4.x and ITK 5.x.
 *=========================================================================*/
#ifndef __itkMplusCompat_h
#define __itkMplusCompat_h

#include "itkVersion.h"

// ── Detect ITK major version ─────────────────────────────────────────────
#if ITK_VERSION_MAJOR >= 5
  #define MPLUS_ITK5 1
#else
  #define MPLUS_ITK5 0
#endif

// ── BSplineTransform name changed in ITK 5 ──────────────────────────────
//    ITK 4: itkBSplineDeformableTransform<>
//    ITK 5: itkBSplineTransform<>
//    Template parameters are unchanged.
#if MPLUS_ITK5
  #include "itkBSplineTransform.h"
  template <typename TScalar, unsigned int NDimensions, unsigned int VSplineOrder>
  using BSplineTransformCompat =
      itk::BSplineTransform<TScalar, NDimensions, VSplineOrder>;
#else
  #include "itkBSplineDeformableTransform.h"
  template <typename TScalar, unsigned int NDimensions, unsigned int VSplineOrder>
  using BSplineTransformCompat =
      itk::BSplineDeformableTransform<TScalar, NDimensions, VSplineOrder>;
#endif

// ── Threading model ──────────────────────────────────────────────────────
//    ITK 4: itkMultiThreader (explicit management)
//    ITK 5: itk::MultiThreaderBase (auto-managed TBB pool)
#if MPLUS_ITK5
  #include "itkMultiThreaderBase.h"
  using ThreaderType = itk::MultiThreaderBase;
#else
  #include "itkMultiThreader.h"
  using ThreaderType = itk::MultiThreader;
#endif

// ── itkTypeMacro → itkOverrideGetNameOfClassMacro ───────────────────────
//    In ITK 5.x the old itkTypeMacro(Self, Superclass) is deprecated.
//    Provide a single macro that resolves to the correct form.
#if MPLUS_ITK5
  #define MPLUS_TYPE_MACRO(Self, Superclass) \
    itkOverrideGetNameOfClassMacro(Self)
#else
  #define MPLUS_TYPE_MACRO(Self, Superclass) \
    itkTypeMacro(Self, Superclass)
#endif

// ── Smart pointer aliases (using vs typedef) ─────────────────────────────
//    Modern alias syntax for C++17 code (both ITK versions).
#define MPLUS_SMART_POINTER_ALIASES(Self) \
  using Pointer      = itk::SmartPointer<Self>;       \
  using ConstPointer = itk::SmartPointer<const Self>;

// ── Deprecation helper ──────────────────────────────────────────────────
#if defined(__cplusplus) && __cplusplus >= 201402L
  #define MPLUS_DEPRECATED [[deprecated]]
  #define MPLUS_DEPRECATED_MSG(msg) [[deprecated(msg)]]
#else
  #define MPLUS_DEPRECATED
  #define MPLUS_DEPRECATED_MSG(msg)
#endif

// ── Misc helpers ─────────────────────────────────────────────────────────
//    ITK 5 removed GetInput<N>() overloads; use ProcessObject::GetInput(n).
//    Wrap in a uniform accessor where needed.

#endif /* __itkMplusCompat_h */
