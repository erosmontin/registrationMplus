# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/cmake-gui

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /DATA/Dropbox/itkSRC/TRIAL/2DFFDMANGF

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /DATA/Dropbox/itkSRC/TRIAL/2DFFDMANGF/LBuild

# Include any dependencies generated for this target.
include CMakeFiles/2DFFDMANGF.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/2DFFDMANGF.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/2DFFDMANGF.dir/flags.make

CMakeFiles/2DFFDMANGF.dir/2DFFDMANGF.cxx.o: CMakeFiles/2DFFDMANGF.dir/flags.make
CMakeFiles/2DFFDMANGF.dir/2DFFDMANGF.cxx.o: ../2DFFDMANGF.cxx
	$(CMAKE_COMMAND) -E cmake_progress_report /DATA/Dropbox/itkSRC/TRIAL/2DFFDMANGF/LBuild/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/2DFFDMANGF.dir/2DFFDMANGF.cxx.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/2DFFDMANGF.dir/2DFFDMANGF.cxx.o -c /DATA/Dropbox/itkSRC/TRIAL/2DFFDMANGF/2DFFDMANGF.cxx

CMakeFiles/2DFFDMANGF.dir/2DFFDMANGF.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/2DFFDMANGF.dir/2DFFDMANGF.cxx.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /DATA/Dropbox/itkSRC/TRIAL/2DFFDMANGF/2DFFDMANGF.cxx > CMakeFiles/2DFFDMANGF.dir/2DFFDMANGF.cxx.i

CMakeFiles/2DFFDMANGF.dir/2DFFDMANGF.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/2DFFDMANGF.dir/2DFFDMANGF.cxx.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /DATA/Dropbox/itkSRC/TRIAL/2DFFDMANGF/2DFFDMANGF.cxx -o CMakeFiles/2DFFDMANGF.dir/2DFFDMANGF.cxx.s

CMakeFiles/2DFFDMANGF.dir/2DFFDMANGF.cxx.o.requires:
.PHONY : CMakeFiles/2DFFDMANGF.dir/2DFFDMANGF.cxx.o.requires

CMakeFiles/2DFFDMANGF.dir/2DFFDMANGF.cxx.o.provides: CMakeFiles/2DFFDMANGF.dir/2DFFDMANGF.cxx.o.requires
	$(MAKE) -f CMakeFiles/2DFFDMANGF.dir/build.make CMakeFiles/2DFFDMANGF.dir/2DFFDMANGF.cxx.o.provides.build
.PHONY : CMakeFiles/2DFFDMANGF.dir/2DFFDMANGF.cxx.o.provides

CMakeFiles/2DFFDMANGF.dir/2DFFDMANGF.cxx.o.provides.build: CMakeFiles/2DFFDMANGF.dir/2DFFDMANGF.cxx.o

# Object files for target 2DFFDMANGF
2DFFDMANGF_OBJECTS = \
"CMakeFiles/2DFFDMANGF.dir/2DFFDMANGF.cxx.o"

# External object files for target 2DFFDMANGF
2DFFDMANGF_EXTERNAL_OBJECTS =

2DFFDMANGF: CMakeFiles/2DFFDMANGF.dir/2DFFDMANGF.cxx.o
2DFFDMANGF: CMakeFiles/2DFFDMANGF.dir/build.make
2DFFDMANGF: /usr/lib/libitksys-4.3.so.1
2DFFDMANGF: /usr/lib/libitkvnl_algo-4.3.so.1
2DFFDMANGF: /usr/lib/libitkvnl-4.3.so.1
2DFFDMANGF: /usr/lib/libitkv3p_netlib-4.3.so.1
2DFFDMANGF: /usr/lib/libITKCommon-4.3.so.1
2DFFDMANGF: /usr/lib/libitkNetlibSlatec-4.3.so.1
2DFFDMANGF: /usr/lib/libITKStatistics-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOImageBase-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOBMP-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOBioRad-4.3.so.1
2DFFDMANGF: /usr/lib/libITKEXPAT-4.3.so.1
2DFFDMANGF: /usr/lib/libitkopenjpeg-4.3.so.1
2DFFDMANGF: /usr/lib/x86_64-linux-gnu/libz.so
2DFFDMANGF: /usr/lib/libITKGDCM-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOGDCM-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOGIPL-4.3.so.1
2DFFDMANGF: /usr/lib/x86_64-linux-gnu/libjpeg.so
2DFFDMANGF: /usr/lib/libITKIOJPEG-4.3.so.1
2DFFDMANGF: /usr/lib/libITKMetaIO-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOMeta-4.3.so.1
2DFFDMANGF: /usr/lib/libITKznz-4.3.so.1
2DFFDMANGF: /usr/lib/libITKniftiio-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIONIFTI-4.3.so.1
2DFFDMANGF: /usr/lib/libITKNrrdIO-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIONRRD-4.3.so.1
2DFFDMANGF: /usr/lib/x86_64-linux-gnu/libpng.so
2DFFDMANGF: /usr/lib/libITKIOPNG-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOStimulate-4.3.so.1
2DFFDMANGF: /usr/lib/x86_64-linux-gnu/libtiff.so
2DFFDMANGF: /usr/lib/libITKIOTIFF-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOVTK-4.3.so.1
2DFFDMANGF: /usr/lib/libITKMesh-4.3.so.1
2DFFDMANGF: /usr/lib/libITKSpatialObjects-4.3.so.1
2DFFDMANGF: /usr/lib/libITKPath-4.3.so.1
2DFFDMANGF: /usr/lib/libITKLabelMap-4.3.so.1
2DFFDMANGF: /usr/lib/libITKQuadEdgeMesh-4.3.so.1
2DFFDMANGF: /usr/lib/libITKOptimizers-4.3.so.1
2DFFDMANGF: /usr/lib/libITKPolynomials-4.3.so.1
2DFFDMANGF: /usr/lib/libITKBiasCorrection-4.3.so.1
2DFFDMANGF: /usr/lib/libITKBioCell-4.3.so.1
2DFFDMANGF: /usr/lib/libITKDICOMParser-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOXML-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOSpatialObjects-4.3.so.1
2DFFDMANGF: /usr/lib/libITKFEM-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOIPL-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOGE-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOSiemens-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOTransformBase-4.3.so.1
2DFFDMANGF: /usr/lib/libitkhdf5_cpp-4.3.so.1
2DFFDMANGF: /usr/lib/libitkhdf5-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOTransformHDF5-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOTransformInsightLegacy-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOTransformMatlab-4.3.so.1
2DFFDMANGF: /usr/lib/libITKKLMRegionGrowing-4.3.so.1
2DFFDMANGF: /usr/lib/libITKVTK-4.3.so.1
2DFFDMANGF: /usr/lib/libITKWatersheds-4.3.so.1
2DFFDMANGF: /usr/lib/libITKDeprecated-4.3.so.1
2DFFDMANGF: /usr/lib/libITKgiftiio-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOMesh-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOCSV-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOHDF5-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOPhilipsREC-4.3.so.1
2DFFDMANGF: /usr/lib/libITKOptimizersv4-4.3.so.1
2DFFDMANGF: /usr/lib/libITKReview-4.3.so.1
2DFFDMANGF: /usr/lib/libITKVideoCore-4.3.so.1
2DFFDMANGF: /usr/lib/libITKVideoIO-4.3.so.1
2DFFDMANGF: /usr/lib/libITKDICOMParser-4.3.so.1
2DFFDMANGF: /usr/lib/libITKgiftiio-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOBMP-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOBioRad-4.3.so.1
2DFFDMANGF: /usr/lib/libitkopenjpeg-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOGDCM-4.3.so.1
2DFFDMANGF: /usr/lib/libITKGDCM-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOGIPL-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOJPEG-4.3.so.1
2DFFDMANGF: /usr/lib/x86_64-linux-gnu/libjpeg.so
2DFFDMANGF: /usr/lib/libITKIOMeta-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIONIFTI-4.3.so.1
2DFFDMANGF: /usr/lib/libITKniftiio-4.3.so.1
2DFFDMANGF: /usr/lib/libITKznz-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIONRRD-4.3.so.1
2DFFDMANGF: /usr/lib/libITKNrrdIO-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOPNG-4.3.so.1
2DFFDMANGF: /usr/lib/x86_64-linux-gnu/libpng.so
2DFFDMANGF: /usr/lib/libITKIOStimulate-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOTIFF-4.3.so.1
2DFFDMANGF: /usr/lib/x86_64-linux-gnu/libtiff.so
2DFFDMANGF: /usr/lib/libITKIOVTK-4.3.so.1
2DFFDMANGF: /usr/lib/libITKLabelMap-4.3.so.1
2DFFDMANGF: /usr/lib/libITKQuadEdgeMesh-4.3.so.1
2DFFDMANGF: /usr/lib/libITKBiasCorrection-4.3.so.1
2DFFDMANGF: /usr/lib/libITKPolynomials-4.3.so.1
2DFFDMANGF: /usr/lib/libITKBioCell-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOSpatialObjects-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOXML-4.3.so.1
2DFFDMANGF: /usr/lib/libITKEXPAT-4.3.so.1
2DFFDMANGF: /usr/lib/libITKFEM-4.3.so.1
2DFFDMANGF: /usr/lib/libITKMetaIO-4.3.so.1
2DFFDMANGF: /usr/lib/libITKOptimizers-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOSiemens-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOGE-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOIPL-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOTransformHDF5-4.3.so.1
2DFFDMANGF: /usr/lib/libitkhdf5_cpp-4.3.so.1
2DFFDMANGF: /usr/lib/libitkhdf5-4.3.so.1
2DFFDMANGF: /usr/lib/x86_64-linux-gnu/libz.so
2DFFDMANGF: /usr/lib/libITKIOTransformInsightLegacy-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOTransformMatlab-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOTransformBase-4.3.so.1
2DFFDMANGF: /usr/lib/libITKKLMRegionGrowing-4.3.so.1
2DFFDMANGF: /usr/lib/libITKVTK-4.3.so.1
2DFFDMANGF: /usr/lib/libITKWatersheds-4.3.so.1
2DFFDMANGF: /usr/lib/libITKSpatialObjects-4.3.so.1
2DFFDMANGF: /usr/lib/libITKMesh-4.3.so.1
2DFFDMANGF: /usr/lib/libITKPath-4.3.so.1
2DFFDMANGF: /usr/lib/libITKStatistics-4.3.so.1
2DFFDMANGF: /usr/lib/libitkNetlibSlatec-4.3.so.1
2DFFDMANGF: /usr/lib/libITKIOImageBase-4.3.so.1
2DFFDMANGF: /usr/lib/libITKVideoCore-4.3.so.1
2DFFDMANGF: /usr/lib/libITKCommon-4.3.so.1
2DFFDMANGF: /usr/lib/libitksys-4.3.so.1
2DFFDMANGF: /usr/lib/libITKVNLInstantiation-4.3.so.1
2DFFDMANGF: /usr/lib/libitkvnl_algo-4.3.so.1
2DFFDMANGF: /usr/lib/libitkv3p_lsqr-4.3.so.1
2DFFDMANGF: /usr/lib/libitkvnl-4.3.so.1
2DFFDMANGF: /usr/lib/libitkvcl-4.3.so.1
2DFFDMANGF: /usr/lib/libitkv3p_netlib-4.3.so.1
2DFFDMANGF: CMakeFiles/2DFFDMANGF.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable 2DFFDMANGF"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/2DFFDMANGF.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/2DFFDMANGF.dir/build: 2DFFDMANGF
.PHONY : CMakeFiles/2DFFDMANGF.dir/build

CMakeFiles/2DFFDMANGF.dir/requires: CMakeFiles/2DFFDMANGF.dir/2DFFDMANGF.cxx.o.requires
.PHONY : CMakeFiles/2DFFDMANGF.dir/requires

CMakeFiles/2DFFDMANGF.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/2DFFDMANGF.dir/cmake_clean.cmake
.PHONY : CMakeFiles/2DFFDMANGF.dir/clean

CMakeFiles/2DFFDMANGF.dir/depend:
	cd /DATA/Dropbox/itkSRC/TRIAL/2DFFDMANGF/LBuild && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /DATA/Dropbox/itkSRC/TRIAL/2DFFDMANGF /DATA/Dropbox/itkSRC/TRIAL/2DFFDMANGF /DATA/Dropbox/itkSRC/TRIAL/2DFFDMANGF/LBuild /DATA/Dropbox/itkSRC/TRIAL/2DFFDMANGF/LBuild /DATA/Dropbox/itkSRC/TRIAL/2DFFDMANGF/LBuild/CMakeFiles/2DFFDMANGF.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/2DFFDMANGF.dir/depend
