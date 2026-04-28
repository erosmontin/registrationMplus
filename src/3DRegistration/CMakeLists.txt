# Output directories.

if(WIN32)
  # Some libraries (e.g. vxl libs) have no dllexport markup, so we can
  # build full shared libraries only with the GNU toolchain.
  if(CMAKE_COMPILER_IS_GNUCXX)
    # CMake adds --enable-all-exports on Cygwin (since Cygwin is
    # supposed to be UNIX-like), but we need to add it explicitly for
    # a native windows build with the MinGW tools.
    if(MINGW)
      set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS
        "-shared -Wl,--export-all-symbols -Wl,--enable-auto-import")
      set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS
        "-shared -Wl,--export-all-symbols -Wl,--enable-auto-import")
      set(CMAKE_EXE_LINKER_FLAGS "-Wl,--enable-auto-import")
    endif(MINGW)
  else(CMAKE_COMPILER_IS_GNUCXX)
   if(BUILD_SHARED_LIBS)
     set(ITK_LIBRARY_BUILD_TYPE "SHARED")
   else(BUILD_SHARED_LIBS)
     set(ITK_LIBRARY_BUILD_TYPE "STATIC")
   endif(BUILD_SHARED_LIBS)
   set(BUILD_SHARED_LIBS OFF)
  endif(CMAKE_COMPILER_IS_GNUCXX)
endif(WIN32)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# OpenMP (portable across MSVC, GCC, Clang). The `#pragma omp parallel for`
# directives in src/Metrics/Mplus are silently ignored without this — on
# Windows MSVC `-fopenmp` is not a valid flag and OpenMP never gets enabled.
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    if(MSVC)
        # MSVC's default `/openmp` is OpenMP 2.0 and does not understand
        # `#pragma omp simd` (used by NGF). `/openmp:experimental` enables the
        # newer SIMD/atomic constructs while remaining backward compatible
        # with the parallel-for loops used elsewhere.
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp:experimental")
        message(STATUS "OpenMP enabled (MSVC): /openmp:experimental")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
        message(STATUS "OpenMP enabled: ${OpenMP_CXX_FLAGS}")
    endif()
else()
    message(WARNING "OpenMP NOT found — derivative-combine loops will run single-threaded.")
endif()


ADD_SUBDIRECTORY(3DRegBsplines/src)
ADD_SUBDIRECTORY(3DRegAffine/src)

