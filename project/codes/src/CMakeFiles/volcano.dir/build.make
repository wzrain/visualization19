# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /opt/cmake-3.12.2/bin/cmake

# The command to remove a file.
RM = /opt/cmake-3.12.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/marshallee/Documents/ETHz/course/visualization/program/final/codes/visualization/project/codes

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/marshallee/Documents/ETHz/course/visualization/program/final/codes/visualization/project/codes/src

# Include any dependencies generated for this target.
include CMakeFiles/volcano.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/volcano.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/volcano.dir/flags.make

CMakeFiles/volcano.dir/volcano.cpp.o: CMakeFiles/volcano.dir/flags.make
CMakeFiles/volcano.dir/volcano.cpp.o: volcano.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/marshallee/Documents/ETHz/course/visualization/program/final/codes/visualization/project/codes/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/volcano.dir/volcano.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/volcano.dir/volcano.cpp.o -c /home/marshallee/Documents/ETHz/course/visualization/program/final/codes/visualization/project/codes/src/volcano.cpp

CMakeFiles/volcano.dir/volcano.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/volcano.dir/volcano.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/marshallee/Documents/ETHz/course/visualization/program/final/codes/visualization/project/codes/src/volcano.cpp > CMakeFiles/volcano.dir/volcano.cpp.i

CMakeFiles/volcano.dir/volcano.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/volcano.dir/volcano.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/marshallee/Documents/ETHz/course/visualization/program/final/codes/visualization/project/codes/src/volcano.cpp -o CMakeFiles/volcano.dir/volcano.cpp.s

# Object files for target volcano
volcano_OBJECTS = \
"CMakeFiles/volcano.dir/volcano.cpp.o"

# External object files for target volcano
volcano_EXTERNAL_OBJECTS =

volcano: CMakeFiles/volcano.dir/volcano.cpp.o
volcano: CMakeFiles/volcano.dir/build.make
volcano: /usr/local/lib/libvtkDomainsChemistryOpenGL2-8.2.so.1
volcano: /usr/local/lib/libvtkFiltersFlowPaths-8.2.so.1
volcano: /usr/local/lib/libvtkFiltersGeneric-8.2.so.1
volcano: /usr/local/lib/libvtkFiltersHyperTree-8.2.so.1
volcano: /usr/local/lib/libvtkFiltersParallelImaging-8.2.so.1
volcano: /usr/local/lib/libvtkFiltersPoints-8.2.so.1
volcano: /usr/local/lib/libvtkFiltersProgrammable-8.2.so.1
volcano: /usr/local/lib/libvtkFiltersSMP-8.2.so.1
volcano: /usr/local/lib/libvtkFiltersSelection-8.2.so.1
volcano: /usr/local/lib/libvtkFiltersTexture-8.2.so.1
volcano: /usr/local/lib/libvtkFiltersTopology-8.2.so.1
volcano: /usr/local/lib/libvtkFiltersVerdict-8.2.so.1
volcano: /usr/local/lib/libvtkGeovisCore-8.2.so.1
volcano: /usr/local/lib/libvtkIOAMR-8.2.so.1
volcano: /usr/local/lib/libvtkIOAsynchronous-8.2.so.1
volcano: /usr/local/lib/libvtkIOCityGML-8.2.so.1
volcano: /usr/local/lib/libvtkIOEnSight-8.2.so.1
volcano: /usr/local/lib/libvtkIOExodus-8.2.so.1
volcano: /usr/local/lib/libvtkIOExportOpenGL2-8.2.so.1
volcano: /usr/local/lib/libvtkIOExportPDF-8.2.so.1
volcano: /usr/local/lib/libvtkIOImport-8.2.so.1
volcano: /usr/local/lib/libvtkIOInfovis-8.2.so.1
volcano: /usr/local/lib/libvtkIOLSDyna-8.2.so.1
volcano: /usr/local/lib/libvtkIOMINC-8.2.so.1
volcano: /usr/local/lib/libvtkIOMovie-8.2.so.1
volcano: /usr/local/lib/libvtkIOPLY-8.2.so.1
volcano: /usr/local/lib/libvtkIOParallel-8.2.so.1
volcano: /usr/local/lib/libvtkIOParallelXML-8.2.so.1
volcano: /usr/local/lib/libvtkIOSQL-8.2.so.1
volcano: /usr/local/lib/libvtkIOSegY-8.2.so.1
volcano: /usr/local/lib/libvtkIOTecplotTable-8.2.so.1
volcano: /usr/local/lib/libvtkIOVeraOut-8.2.so.1
volcano: /usr/local/lib/libvtkIOVideo-8.2.so.1
volcano: /usr/local/lib/libvtkImagingMorphological-8.2.so.1
volcano: /usr/local/lib/libvtkImagingStatistics-8.2.so.1
volcano: /usr/local/lib/libvtkImagingStencil-8.2.so.1
volcano: /usr/local/lib/libvtkInteractionImage-8.2.so.1
volcano: /usr/local/lib/libvtkRenderingContextOpenGL2-8.2.so.1
volcano: /usr/local/lib/libvtkRenderingImage-8.2.so.1
volcano: /usr/local/lib/libvtkRenderingLOD-8.2.so.1
volcano: /usr/local/lib/libvtkRenderingVolumeOpenGL2-8.2.so.1
volcano: /usr/local/lib/libvtkViewsContext2D-8.2.so.1
volcano: /usr/local/lib/libvtkViewsInfovis-8.2.so.1
volcano: /usr/local/lib/libvtkDomainsChemistry-8.2.so.1
volcano: /usr/local/lib/libvtkverdict-8.2.so.1
volcano: /usr/local/lib/libvtkproj-8.2.so.1
volcano: /usr/local/lib/libvtkFiltersAMR-8.2.so.1
volcano: /usr/local/lib/libvtkpugixml-8.2.so.1
volcano: /usr/local/lib/libvtkIOExport-8.2.so.1
volcano: /usr/local/lib/libvtkRenderingGL2PSOpenGL2-8.2.so.1
volcano: /usr/local/lib/libvtkgl2ps-8.2.so.1
volcano: /usr/local/lib/libvtklibharu-8.2.so.1
volcano: /usr/local/lib/libvtklibxml2-8.2.so.1
volcano: /usr/local/lib/libvtktheora-8.2.so.1
volcano: /usr/local/lib/libvtkogg-8.2.so.1
volcano: /usr/local/lib/libvtkFiltersParallel-8.2.so.1
volcano: /usr/local/lib/libvtkexodusII-8.2.so.1
volcano: /usr/local/lib/libvtkIOGeometry-8.2.so.1
volcano: /usr/local/lib/libvtkIONetCDF-8.2.so.1
volcano: /usr/local/lib/libvtkNetCDF-8.2.so.1
volcano: /usr/local/lib/libvtkjsoncpp-8.2.so.1
volcano: /usr/local/lib/libvtkParallelCore-8.2.so.1
volcano: /usr/local/lib/libvtkIOLegacy-8.2.so.1
volcano: /usr/local/lib/libvtksqlite-8.2.so.1
volcano: /usr/local/lib/libvtkhdf5-8.2.so.1
volcano: /usr/local/lib/libvtkhdf5_hl-8.2.so.1
volcano: /usr/local/lib/libvtkRenderingOpenGL2-8.2.so.1
volcano: /usr/local/lib/libvtkglew-8.2.so.1
volcano: /usr/lib/x86_64-linux-gnu/libSM.so
volcano: /usr/lib/x86_64-linux-gnu/libICE.so
volcano: /usr/lib/x86_64-linux-gnu/libX11.so
volcano: /usr/lib/x86_64-linux-gnu/libXext.so
volcano: /usr/lib/x86_64-linux-gnu/libXt.so
volcano: /usr/local/lib/libvtkImagingMath-8.2.so.1
volcano: /usr/local/lib/libvtkChartsCore-8.2.so.1
volcano: /usr/local/lib/libvtkRenderingContext2D-8.2.so.1
volcano: /usr/local/lib/libvtkFiltersImaging-8.2.so.1
volcano: /usr/local/lib/libvtkInfovisLayout-8.2.so.1
volcano: /usr/local/lib/libvtkInfovisCore-8.2.so.1
volcano: /usr/local/lib/libvtkViewsCore-8.2.so.1
volcano: /usr/local/lib/libvtkInteractionWidgets-8.2.so.1
volcano: /usr/local/lib/libvtkFiltersHybrid-8.2.so.1
volcano: /usr/local/lib/libvtkImagingGeneral-8.2.so.1
volcano: /usr/local/lib/libvtkImagingSources-8.2.so.1
volcano: /usr/local/lib/libvtkFiltersModeling-8.2.so.1
volcano: /usr/local/lib/libvtkImagingHybrid-8.2.so.1
volcano: /usr/local/lib/libvtkIOImage-8.2.so.1
volcano: /usr/local/lib/libvtkDICOMParser-8.2.so.1
volcano: /usr/local/lib/libvtkmetaio-8.2.so.1
volcano: /usr/local/lib/libvtkjpeg-8.2.so.1
volcano: /usr/local/lib/libvtkpng-8.2.so.1
volcano: /usr/local/lib/libvtktiff-8.2.so.1
volcano: /usr/local/lib/libvtkInteractionStyle-8.2.so.1
volcano: /usr/local/lib/libvtkFiltersExtraction-8.2.so.1
volcano: /usr/local/lib/libvtkFiltersStatistics-8.2.so.1
volcano: /usr/local/lib/libvtkImagingFourier-8.2.so.1
volcano: /usr/local/lib/libvtkRenderingAnnotation-8.2.so.1
volcano: /usr/local/lib/libvtkImagingColor-8.2.so.1
volcano: /usr/local/lib/libvtkRenderingVolume-8.2.so.1
volcano: /usr/local/lib/libvtkImagingCore-8.2.so.1
volcano: /usr/local/lib/libvtkIOXML-8.2.so.1
volcano: /usr/local/lib/libvtkIOXMLParser-8.2.so.1
volcano: /usr/local/lib/libvtkIOCore-8.2.so.1
volcano: /usr/local/lib/libvtkdoubleconversion-8.2.so.1
volcano: /usr/local/lib/libvtklz4-8.2.so.1
volcano: /usr/local/lib/libvtklzma-8.2.so.1
volcano: /usr/local/lib/libvtkexpat-8.2.so.1
volcano: /usr/local/lib/libvtkRenderingLabel-8.2.so.1
volcano: /usr/local/lib/libvtkRenderingFreeType-8.2.so.1
volcano: /usr/local/lib/libvtkRenderingCore-8.2.so.1
volcano: /usr/local/lib/libvtkCommonColor-8.2.so.1
volcano: /usr/local/lib/libvtkFiltersGeometry-8.2.so.1
volcano: /usr/local/lib/libvtkFiltersSources-8.2.so.1
volcano: /usr/local/lib/libvtkFiltersGeneral-8.2.so.1
volcano: /usr/local/lib/libvtkCommonComputationalGeometry-8.2.so.1
volcano: /usr/local/lib/libvtkFiltersCore-8.2.so.1
volcano: /usr/local/lib/libvtkCommonExecutionModel-8.2.so.1
volcano: /usr/local/lib/libvtkCommonDataModel-8.2.so.1
volcano: /usr/local/lib/libvtkCommonMisc-8.2.so.1
volcano: /usr/local/lib/libvtkCommonSystem-8.2.so.1
volcano: /usr/local/lib/libvtksys-8.2.so.1
volcano: /usr/local/lib/libvtkCommonTransforms-8.2.so.1
volcano: /usr/local/lib/libvtkCommonMath-8.2.so.1
volcano: /usr/local/lib/libvtkCommonCore-8.2.so.1
volcano: /usr/local/lib/libvtkfreetype-8.2.so.1
volcano: /usr/local/lib/libvtkzlib-8.2.so.1
volcano: CMakeFiles/volcano.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/marshallee/Documents/ETHz/course/visualization/program/final/codes/visualization/project/codes/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable volcano"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/volcano.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/volcano.dir/build: volcano

.PHONY : CMakeFiles/volcano.dir/build

CMakeFiles/volcano.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/volcano.dir/cmake_clean.cmake
.PHONY : CMakeFiles/volcano.dir/clean

CMakeFiles/volcano.dir/depend:
	cd /home/marshallee/Documents/ETHz/course/visualization/program/final/codes/visualization/project/codes/src && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/marshallee/Documents/ETHz/course/visualization/program/final/codes/visualization/project/codes /home/marshallee/Documents/ETHz/course/visualization/program/final/codes/visualization/project/codes /home/marshallee/Documents/ETHz/course/visualization/program/final/codes/visualization/project/codes/src /home/marshallee/Documents/ETHz/course/visualization/program/final/codes/visualization/project/codes/src /home/marshallee/Documents/ETHz/course/visualization/program/final/codes/visualization/project/codes/src/CMakeFiles/volcano.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/volcano.dir/depend

