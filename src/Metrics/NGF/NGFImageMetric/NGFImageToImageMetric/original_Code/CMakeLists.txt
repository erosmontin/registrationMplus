
MACRO(NEW_TEST name)
  IF (Boost_FOUND) 
    ADD_EXECUTABLE(test${name} test${name}.cpp)
    TARGET_LINK_LIBRARIES(test${name} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ITKNumerics ITKAlgorithms )
    ADD_TEST(${name} test${name})
  ENDIF(Boost_FOUND) 
ENDMACRO(NEW_TEST src)


NEW_TEST(ImageToNGFFilter)
NEW_TEST(NGFMetricKernel)
NEW_TEST(NormalizedGradientFieldImageToImageMetric)


