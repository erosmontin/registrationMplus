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

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -fpermissive")

ADD_SUBDIRECTORY(3DRecMultiTushort/src)
ADD_SUBDIRECTORY(3DRecMultiTshort/src)
ADD_SUBDIRECTORY(3DRecMultiTuchar/src)
ADD_SUBDIRECTORY(3DRecMultiRushort/src)
ADD_SUBDIRECTORY(3DRecMultiRshort/src)
ADD_SUBDIRECTORY(3DRecMultiSushort/src)
ADD_SUBDIRECTORY(3DRecMultiSshort/src)
ADD_SUBDIRECTORY(3DRecMultiAushort/src)
ADD_SUBDIRECTORY(3DRecMultiAshort/src)

#ADD_SUBDIRECTORY(3DTto3Dversortfm/src)
#ADD_SUBDIRECTORY(3DRecBspineMANGF2short/src)

ADD_SUBDIRECTORY(3DRecBsplinesMANGFushort/src)
ADD_SUBDIRECTORY(3DRecBsplinesushort/src)


#ADD_SUBDIRECTORY(2DRecMultiAchar/src) #version1
#ADD_SUBDIRECTORY(2DRecMultiAushort/src) #version1
#ADD_SUBDIRECTORY(2DRecMultiBChar/src)
#ADD_SUBDIRECTORY(2DRecMultiBuShort/src)
#ADD_SUBDIRECTORY(2DRecMultiRchar/src)
#ADD_SUBDIRECTORY(2DRecMultiRushort/src)
#ADD_SUBDIRECTORY(2DRecMultiSchar/src)
#ADD_SUBDIRECTORY(2DRecMultiSushort/src)
#ADD_SUBDIRECTORY(2DRecMultiTchar/src)

 #ADD_SUBDIRECTORY(2DRecDemonsChar/src)
#ADD_SUBDIRECTORY(UN/src)
