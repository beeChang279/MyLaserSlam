# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/beechang/下载/HW1/basicTransformStudy

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/beechang/下载/HW1/basicTransformStudy

# Include any dependencies generated for this target.
include CMakeFiles/basicTransformStudy.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/basicTransformStudy.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/basicTransformStudy.dir/flags.make

CMakeFiles/basicTransformStudy.dir/basic_transform_study.cpp.o: CMakeFiles/basicTransformStudy.dir/flags.make
CMakeFiles/basicTransformStudy.dir/basic_transform_study.cpp.o: basic_transform_study.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/beechang/下载/HW1/basicTransformStudy/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/basicTransformStudy.dir/basic_transform_study.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/basicTransformStudy.dir/basic_transform_study.cpp.o -c /home/beechang/下载/HW1/basicTransformStudy/basic_transform_study.cpp

CMakeFiles/basicTransformStudy.dir/basic_transform_study.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/basicTransformStudy.dir/basic_transform_study.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/beechang/下载/HW1/basicTransformStudy/basic_transform_study.cpp > CMakeFiles/basicTransformStudy.dir/basic_transform_study.cpp.i

CMakeFiles/basicTransformStudy.dir/basic_transform_study.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/basicTransformStudy.dir/basic_transform_study.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/beechang/下载/HW1/basicTransformStudy/basic_transform_study.cpp -o CMakeFiles/basicTransformStudy.dir/basic_transform_study.cpp.s

CMakeFiles/basicTransformStudy.dir/basic_transform_study.cpp.o.requires:

.PHONY : CMakeFiles/basicTransformStudy.dir/basic_transform_study.cpp.o.requires

CMakeFiles/basicTransformStudy.dir/basic_transform_study.cpp.o.provides: CMakeFiles/basicTransformStudy.dir/basic_transform_study.cpp.o.requires
	$(MAKE) -f CMakeFiles/basicTransformStudy.dir/build.make CMakeFiles/basicTransformStudy.dir/basic_transform_study.cpp.o.provides.build
.PHONY : CMakeFiles/basicTransformStudy.dir/basic_transform_study.cpp.o.provides

CMakeFiles/basicTransformStudy.dir/basic_transform_study.cpp.o.provides.build: CMakeFiles/basicTransformStudy.dir/basic_transform_study.cpp.o


# Object files for target basicTransformStudy
basicTransformStudy_OBJECTS = \
"CMakeFiles/basicTransformStudy.dir/basic_transform_study.cpp.o"

# External object files for target basicTransformStudy
basicTransformStudy_EXTERNAL_OBJECTS =

basicTransformStudy: CMakeFiles/basicTransformStudy.dir/basic_transform_study.cpp.o
basicTransformStudy: CMakeFiles/basicTransformStudy.dir/build.make
basicTransformStudy: CMakeFiles/basicTransformStudy.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/beechang/下载/HW1/basicTransformStudy/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable basicTransformStudy"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/basicTransformStudy.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/basicTransformStudy.dir/build: basicTransformStudy

.PHONY : CMakeFiles/basicTransformStudy.dir/build

CMakeFiles/basicTransformStudy.dir/requires: CMakeFiles/basicTransformStudy.dir/basic_transform_study.cpp.o.requires

.PHONY : CMakeFiles/basicTransformStudy.dir/requires

CMakeFiles/basicTransformStudy.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/basicTransformStudy.dir/cmake_clean.cmake
.PHONY : CMakeFiles/basicTransformStudy.dir/clean

CMakeFiles/basicTransformStudy.dir/depend:
	cd /home/beechang/下载/HW1/basicTransformStudy && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/beechang/下载/HW1/basicTransformStudy /home/beechang/下载/HW1/basicTransformStudy /home/beechang/下载/HW1/basicTransformStudy /home/beechang/下载/HW1/basicTransformStudy /home/beechang/下载/HW1/basicTransformStudy/CMakeFiles/basicTransformStudy.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/basicTransformStudy.dir/depend

