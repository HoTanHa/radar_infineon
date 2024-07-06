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

# Include any dependencies generated for this target.
include examples/c/BGT60LTR11AIP/advanced_motion_sensing/CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/depend.make

# Include the progress variables for this target.
include examples/c/BGT60LTR11AIP/advanced_motion_sensing/CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/progress.make

# Include the compile flags for this target's objects.
include examples/c/BGT60LTR11AIP/advanced_motion_sensing/CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/flags.make

examples/c/BGT60LTR11AIP/advanced_motion_sensing/CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/advanced_motion_sensing.c.o: examples/c/BGT60LTR11AIP/advanced_motion_sensing/CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/flags.make
examples/c/BGT60LTR11AIP/advanced_motion_sensing/CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/advanced_motion_sensing.c.o: ../examples/c/BGT60LTR11AIP/advanced_motion_sensing/advanced_motion_sensing.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/htha/Infineon/radar_sdk/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object examples/c/BGT60LTR11AIP/advanced_motion_sensing/CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/advanced_motion_sensing.c.o"
	cd /home/htha/Infineon/radar_sdk/build/examples/c/BGT60LTR11AIP/advanced_motion_sensing && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/advanced_motion_sensing.c.o   -c /home/htha/Infineon/radar_sdk/examples/c/BGT60LTR11AIP/advanced_motion_sensing/advanced_motion_sensing.c

examples/c/BGT60LTR11AIP/advanced_motion_sensing/CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/advanced_motion_sensing.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/advanced_motion_sensing.c.i"
	cd /home/htha/Infineon/radar_sdk/build/examples/c/BGT60LTR11AIP/advanced_motion_sensing && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/htha/Infineon/radar_sdk/examples/c/BGT60LTR11AIP/advanced_motion_sensing/advanced_motion_sensing.c > CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/advanced_motion_sensing.c.i

examples/c/BGT60LTR11AIP/advanced_motion_sensing/CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/advanced_motion_sensing.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/advanced_motion_sensing.c.s"
	cd /home/htha/Infineon/radar_sdk/build/examples/c/BGT60LTR11AIP/advanced_motion_sensing && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/htha/Infineon/radar_sdk/examples/c/BGT60LTR11AIP/advanced_motion_sensing/advanced_motion_sensing.c -o CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/advanced_motion_sensing.c.s

# Object files for target BGT60LTR11AIP_advanced_motion_sensing
BGT60LTR11AIP_advanced_motion_sensing_OBJECTS = \
"CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/advanced_motion_sensing.c.o"

# External object files for target BGT60LTR11AIP_advanced_motion_sensing
BGT60LTR11AIP_advanced_motion_sensing_EXTERNAL_OBJECTS =

bin/BGT60LTR11AIP_advanced_motion_sensing: examples/c/BGT60LTR11AIP/advanced_motion_sensing/CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/advanced_motion_sensing.c.o
bin/BGT60LTR11AIP_advanced_motion_sensing: examples/c/BGT60LTR11AIP/advanced_motion_sensing/CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/build.make
bin/BGT60LTR11AIP_advanced_motion_sensing: bin/libsdk_ltr11.so
bin/BGT60LTR11AIP_advanced_motion_sensing: 3rd_party/libs/argparse/libargparse.a
bin/BGT60LTR11AIP_advanced_motion_sensing: ../libs/linux_x64/libsdk_advanced_motion_sensing.so
bin/BGT60LTR11AIP_advanced_motion_sensing: bin/libsdk_radar_device_common.so
bin/BGT60LTR11AIP_advanced_motion_sensing: bin/libsdk_base.so
bin/BGT60LTR11AIP_advanced_motion_sensing: bin/libstrata_shared.so
bin/BGT60LTR11AIP_advanced_motion_sensing: external/strata/contrib/pugixml/libpugixml.a
bin/BGT60LTR11AIP_advanced_motion_sensing: examples/c/BGT60LTR11AIP/advanced_motion_sensing/CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/htha/Infineon/radar_sdk/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../../../bin/BGT60LTR11AIP_advanced_motion_sensing"
	cd /home/htha/Infineon/radar_sdk/build/examples/c/BGT60LTR11AIP/advanced_motion_sensing && /usr/bin/cmake -E copy /home/htha/Infineon/radar_sdk/libs/linux_x64/libsdk_advanced_motion_sensing.so /home/htha/Infineon/radar_sdk/build/bin
	cd /home/htha/Infineon/radar_sdk/build/examples/c/BGT60LTR11AIP/advanced_motion_sensing && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/c/BGT60LTR11AIP/advanced_motion_sensing/CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/build: bin/BGT60LTR11AIP_advanced_motion_sensing

.PHONY : examples/c/BGT60LTR11AIP/advanced_motion_sensing/CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/build

examples/c/BGT60LTR11AIP/advanced_motion_sensing/CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/clean:
	cd /home/htha/Infineon/radar_sdk/build/examples/c/BGT60LTR11AIP/advanced_motion_sensing && $(CMAKE_COMMAND) -P CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/cmake_clean.cmake
.PHONY : examples/c/BGT60LTR11AIP/advanced_motion_sensing/CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/clean

examples/c/BGT60LTR11AIP/advanced_motion_sensing/CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/depend:
	cd /home/htha/Infineon/radar_sdk/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/htha/Infineon/radar_sdk /home/htha/Infineon/radar_sdk/examples/c/BGT60LTR11AIP/advanced_motion_sensing /home/htha/Infineon/radar_sdk/build /home/htha/Infineon/radar_sdk/build/examples/c/BGT60LTR11AIP/advanced_motion_sensing /home/htha/Infineon/radar_sdk/build/examples/c/BGT60LTR11AIP/advanced_motion_sensing/CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/c/BGT60LTR11AIP/advanced_motion_sensing/CMakeFiles/BGT60LTR11AIP_advanced_motion_sensing.dir/depend
