# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
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
CMAKE_COMMAND = /home/lijing/softwares/cmake/2.8.5/bin/cmake

# The command to remove a file.
RM = /home/lijing/softwares/cmake/2.8.5/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /home/lijing/softwares/cmake/2.8.5/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/lijing/cepc/v01-17-05/Arbor_KD_3.3

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/lijing/cepc/v01-17-05/Arbor_KD_3.3/build

# Utility rule file for ExperimentalMemCheck.

CMakeFiles/ExperimentalMemCheck:
	/home/lijing/softwares/cmake/2.8.5/bin/ctest -D ExperimentalMemCheck

ExperimentalMemCheck: CMakeFiles/ExperimentalMemCheck
ExperimentalMemCheck: CMakeFiles/ExperimentalMemCheck.dir/build.make
.PHONY : ExperimentalMemCheck

# Rule to build all files generated by this target.
CMakeFiles/ExperimentalMemCheck.dir/build: ExperimentalMemCheck
.PHONY : CMakeFiles/ExperimentalMemCheck.dir/build

CMakeFiles/ExperimentalMemCheck.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ExperimentalMemCheck.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ExperimentalMemCheck.dir/clean

CMakeFiles/ExperimentalMemCheck.dir/depend:
	cd /home/lijing/cepc/v01-17-05/Arbor_KD_3.3/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/lijing/cepc/v01-17-05/Arbor_KD_3.3 /home/lijing/cepc/v01-17-05/Arbor_KD_3.3 /home/lijing/cepc/v01-17-05/Arbor_KD_3.3/build /home/lijing/cepc/v01-17-05/Arbor_KD_3.3/build /home/lijing/cepc/v01-17-05/Arbor_KD_3.3/build/CMakeFiles/ExperimentalMemCheck.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ExperimentalMemCheck.dir/depend

