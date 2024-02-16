
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/Code")

MACRO(ITKEXE name sources)
  ADD_EXECUTABLE(${name} ${sources})
  TARGET_LINK_LIBRARIES(${name} ITKNumerics ITKAlgorithms ITKIO)
ENDMACRO(ITKEXE name sources)

ITKEXE(NGFGradientMap NGFGradientMap.cxx)
ITKEXE(ImageToNGF ImageToNGF.cxx)
ITKEXE(MultiResImageRegistration2DTranslateOnly MultiResImageRegistration2DTranslateOnly.cxx)
ITKEXE(MultiResImageRegistration2D MultiResImageRegistration2D.cxx)


