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
include examples/c/BGT60ATR24C/range_doppler_map/CMakeFiles/BGT60ATR24C_range_doppler_map.dir/depend.make

# Include the progress variables for this target.
include examples/c/BGT60ATR24C/range_doppler_map/CMakeFiles/BGT60ATR24C_range_doppler_map.dir/progress.make

# Include the compile flags for this target's objects.
include examples/c/BGT60ATR24C/range_doppler_map/CMakeFiles/BGT60ATR24C_range_doppler_map.dir/flags.make

examples/c/BGT60ATR24C/range_doppler_map/CMakeFiles/BGT60ATR24C_range_doppler_map.dir/range_doppler_map.c.o: examples/c/BGT60ATR24C/range_doppler_map/CMakeFiles/BGT60ATR24C_range_doppler_map.dir/flags.make
examples/c/BGT60ATR24C/range_doppler_map/CMakeFiles/BGT60ATR24C_range_doppler_map.dir/range_doppler_map.c.o: ../examples/c/BGT60ATR24C/range_doppler_map/range_doppler_map.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/htha/Infineon/radar_sdk/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object examples/c/BGT60ATR24C/range_doppler_map/CMakeFiles/BGT60ATR24C_range_doppler_map.dir/range_doppler_map.c.o"
	cd /home/htha/Infineon/radar_sdk/build/examples/c/BGT60ATR24C/range_doppler_map && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/BGT60ATR24C_range_doppler_map.dir/range_doppler_map.c.o   -c /home/htha/Infineon/radar_sdk/examples/c/BGT60ATR24C/range_doppler_map/range_doppler_map.c

examples/c/BGT60ATR24C/range_doppler_map/CMakeFiles/BGT60ATR24C_range_doppler_map.dir/range_doppler_map.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/BGT60ATR24C_range_doppler_map.dir/range_doppler_map.c.i"
	cd /home/htha/Infineon/radar_sdk/build/examples/c/BGT60ATR24C/range_doppler_map && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/htha/Infineon/radar_sdk/examples/c/BGT60ATR24C/range_doppler_map/range_doppler_map.c > CMakeFiles/BGT60ATR24C_range_doppler_map.dir/range_doppler_map.c.i

examples/c/BGT60ATR24C/range_doppler_map/CMakeFiles/BGT60ATR24C_range_doppler_map.dir/range_doppler_map.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/BGT60ATR24C_range_doppler_map.dir/range_doppler_map.c.s"
	cd /home/htha/Infineon/radar_sdk/build/examples/c/BGT60ATR24C/range_doppler_map && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/htha/Infineon/radar_sdk/examples/c/BGT60ATR24C/range_doppler_map/range_doppler_map.c -o CMakeFiles/BGT60ATR24C_range_doppler_map.dir/range_doppler_map.c.s

# Object files for target BGT60ATR24C_range_doppler_map
BGT60ATR24C_range_doppler_map_OBJECTS = \
"CMakeFiles/BGT60ATR24C_range_doppler_map.dir/range_doppler_map.c.o"

# External object files for target BGT60ATR24C_range_doppler_map
BGT60ATR24C_range_doppler_map_EXTERNAL_OBJECTS =

bin/BGT60ATR24C_range_doppler_map: examples/c/BGT60ATR24C/range_doppler_map/CMakeFiles/BGT60ATR24C_range_doppler_map.dir/range_doppler_map.c.o
bin/BGT60ATR24C_range_doppler_map: examples/c/BGT60ATR24C/range_doppler_map/CMakeFiles/BGT60ATR24C_range_doppler_map.dir/build.make
bin/BGT60ATR24C_range_doppler_map: examples/c/BGT60ATR24C/common/libBGT60ATR24C_common.a
bin/BGT60ATR24C_range_doppler_map: bin/libsdk_radar.so
bin/BGT60ATR24C_range_doppler_map: bin/libsdk_avian.so
bin/BGT60ATR24C_range_doppler_map: bin/libsdk_fmcw.so
bin/BGT60ATR24C_range_doppler_map: bin/libsdk_radar_device_common.so
bin/BGT60ATR24C_range_doppler_map: bin/libsdk_algo.so
bin/BGT60ATR24C_range_doppler_map: bin/libsdk_base.so
bin/BGT60ATR24C_range_doppler_map: bin/libstrata_shared.so
bin/BGT60ATR24C_range_doppler_map: external/strata/contrib/pugixml/libpugixml.a
bin/BGT60ATR24C_range_doppler_map: 3rd_party/libs/argparse/libargparse.a
bin/BGT60ATR24C_range_doppler_map: examples/c/BGT60ATR24C/range_doppler_map/CMakeFiles/BGT60ATR24C_range_doppler_map.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/htha/Infineon/radar_sdk/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../../../bin/BGT60ATR24C_range_doppler_map"
	cd /home/htha/Infineon/radar_sdk/build/examples/c/BGT60ATR24C/range_doppler_map && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/BGT60ATR24C_range_doppler_map.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/c/BGT60ATR24C/range_doppler_map/CMakeFiles/BGT60ATR24C_range_doppler_map.dir/build: bin/BGT60ATR24C_range_doppler_map

.PHONY : examples/c/BGT60ATR24C/range_doppler_map/CMakeFiles/BGT60ATR24C_range_doppler_map.dir/build

examples/c/BGT60ATR24C/range_doppler_map/CMakeFiles/BGT60ATR24C_range_doppler_map.dir/clean:
	cd /home/htha/Infineon/radar_sdk/build/examples/c/BGT60ATR24C/range_doppler_map && $(CMAKE_COMMAND) -P CMakeFiles/BGT60ATR24C_range_doppler_map.dir/cmake_clean.cmake
.PHONY : examples/c/BGT60ATR24C/range_doppler_map/CMakeFiles/BGT60ATR24C_range_doppler_map.dir/clean

examples/c/BGT60ATR24C/range_doppler_map/CMakeFiles/BGT60ATR24C_range_doppler_map.dir/depend:
	cd /home/htha/Infineon/radar_sdk/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/htha/Infineon/radar_sdk /home/htha/Infineon/radar_sdk/examples/c/BGT60ATR24C/range_doppler_map /home/htha/Infineon/radar_sdk/build /home/htha/Infineon/radar_sdk/build/examples/c/BGT60ATR24C/range_doppler_map /home/htha/Infineon/radar_sdk/build/examples/c/BGT60ATR24C/range_doppler_map/CMakeFiles/BGT60ATR24C_range_doppler_map.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/c/BGT60ATR24C/range_doppler_map/CMakeFiles/BGT60ATR24C_range_doppler_map.dir/depend
