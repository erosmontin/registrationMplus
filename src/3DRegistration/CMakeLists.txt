# ── 3DRegistration sub-project ─────────────────────────────────────────────

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
    set(BUILD_SHARED_LIBS OFF)
  endif()
endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_subdirectory(3DRegBsplines/src)
add_subdirectory(3DRegAffine/src)

