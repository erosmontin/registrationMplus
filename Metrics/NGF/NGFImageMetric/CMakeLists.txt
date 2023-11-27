PROJECT("ITKNGF")

CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0 FATAL_ERROR)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

IF (NOT WIN32)
  SET(SKIP_ITK_BUILD_SETTINGS_FILE "ON")
ENDIF (NOT WIN32)

FIND_PACKAGE(ITK REQUIRED)

  INCLUDE(${ITK_USE_FILE})

FIND_PACKAGE( Boost 1.34.1 COMPONENTS unit_test_framework)
IF (Boost_FOUND) 
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
  ENABLE_TESTING()
  add_definitions(-DHAVE_BOOST)
ELSE (Boost_FOUND) 
  MESSAGE(WARN "Boost testing framework not found, test will not be compiled")
ENDIF(Boost_FOUND) 

ADD_SUBDIRECTORY(Code)
ADD_SUBDIRECTORY(Examples)
