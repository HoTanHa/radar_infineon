# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/htha/Infineon/radar_sdk

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/htha/Infineon/radar_sdk/build

# Utility rule file for NightlyTest.

# Include the progress variables for this target.
include external/strata/contrib/pugixml/CMakeFiles/NightlyTest.dir/progress.make

external/strata/contrib/pugixml/CMakeFiles/NightlyTest:
	cd /home/htha/Infineon/radar_sdk/build/external/strata/contrib/pugixml && /usr/bin/ctest -D NightlyTest

NightlyTest: external/strata/contrib/pugixml/CMakeFiles/NightlyTest
NightlyTest: external/strata/contrib/pugixml/CMakeFiles/NightlyTest.dir/build.make

.PHONY : NightlyTest

# Rule to build all files generated by this target.
external/strata/contrib/pugixml/CMakeFiles/NightlyTest.dir/build: NightlyTest

.PHONY : external/strata/contrib/pugixml/CMakeFiles/NightlyTest.dir/build

external/strata/contrib/pugixml/CMakeFiles/NightlyTest.dir/clean:
	cd /home/htha/Infineon/radar_sdk/build/external/strata/contrib/pugixml && $(CMAKE_COMMAND) -P CMakeFiles/NightlyTest.dir/cmake_clean.cmake
.PHONY : external/strata/contrib/pugixml/CMakeFiles/NightlyTest.dir/clean

external/strata/contrib/pugixml/CMakeFiles/NightlyTest.dir/depend:
	cd /home/htha/Infineon/radar_sdk/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/htha/Infineon/radar_sdk /home/htha/Infineon/radar_sdk/external/strata/contrib/pugixml /home/htha/Infineon/radar_sdk/build /home/htha/Infineon/radar_sdk/build/external/strata/contrib/pugixml /home/htha/Infineon/radar_sdk/build/external/strata/contrib/pugixml/CMakeFiles/NightlyTest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/strata/contrib/pugixml/CMakeFiles/NightlyTest.dir/depend

