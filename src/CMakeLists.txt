cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

project(RegistrationSuite LANGUAGES CXX)

# ── C++17 standard (required for ITK 5.x) ────────────────────────────────
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Options ───────────────────────────────────────────────────────────────
option(USE_CUDA    "Enable CUDA acceleration"      OFF)
option(BUILD_PYTHON "Build Python pybind11 module"  OFF)
option(BUILD_TESTS "Build test executables"         OFF)

# ── Find ITK (5.x preferred; 4.x still supported via itkMplusCompat.h) ──
# Pass -DITK_DIR=/path to cmake to select a specific installation.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
message(STATUS "ITK version: ${ITK_VERSION}")

# ── CUDA (optional) ──────────────────────────────────────────────────────
if(USE_CUDA)
  set(CMAKE_CUDA_FLAGS "-allow-unsupported-compiler --expt-extended-lambda -Xcompiler -fPIC" CACHE STRING "" FORCE)
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
  set(CMAKE_CUDA_ARCHITECTURES "50;60;70;75;80;86;89;90")
  add_compile_definitions(USE_CUDA)
  add_subdirectory(Metrics/Mplus/cuda)
endif()

# ── Platform-specific shared library settings ─────────────────────────────
if(WIN32)
  if(CMAKE_COMPILER_IS_GNUCXX)
    if(MINGW)
      set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS
        "-shared -Wl,--export-all-symbols -Wl,--enable-auto-import")
      set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS
        "-shared -Wl,--export-all-symbols -Wl,--enable-auto-import")
      set(CMAKE_EXE_LINKER_FLAGS "-Wl,--enable-auto-import")
    endif()
  else()
    if(BUILD_SHARED_LIBS)
      set(ITK_LIBRARY_BUILD_TYPE "SHARED")
    else()
      set(ITK_LIBRARY_BUILD_TYPE "STATIC")
    endif()
    set(BUILD_SHARED_LIBS ON)
  endif()
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# ── Sub-projects ──────────────────────────────────────────────────────────
add_subdirectory(3DRegistration)

if(BUILD_TESTS)
  add_subdirectory(test)
endif()

if(BUILD_PYTHON)
  add_subdirectory(${CMAKE_SOURCE_DIR}/../python ${CMAKE_BINARY_DIR}/python)
endif()


