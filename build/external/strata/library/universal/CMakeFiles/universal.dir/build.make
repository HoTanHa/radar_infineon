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
include external/strata/library/universal/CMakeFiles/universal.dir/depend.make

# Include the progress variables for this target.
include external/strata/library/universal/CMakeFiles/universal.dir/progress.make

# Include the compile flags for this target's objects.
include external/strata/library/universal/CMakeFiles/universal.dir/flags.make

external/strata/library/universal/CMakeFiles/universal.dir/error_strings.c.o: external/strata/library/universal/CMakeFiles/universal.dir/flags.make
external/strata/library/universal/CMakeFiles/universal.dir/error_strings.c.o: ../external/strata/library/universal/error_strings.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/htha/Infineon/radar_sdk/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object external/strata/library/universal/CMakeFiles/universal.dir/error_strings.c.o"
	cd /home/htha/Infineon/radar_sdk/build/external/strata/library/universal && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/universal.dir/error_strings.c.o   -c /home/htha/Infineon/radar_sdk/external/strata/library/universal/error_strings.c

external/strata/library/universal/CMakeFiles/universal.dir/error_strings.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/universal.dir/error_strings.c.i"
	cd /home/htha/Infineon/radar_sdk/build/external/strata/library/universal && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/htha/Infineon/radar_sdk/external/strata/library/universal/error_strings.c > CMakeFiles/universal.dir/error_strings.c.i

external/strata/library/universal/CMakeFiles/universal.dir/error_strings.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/universal.dir/error_strings.c.s"
	cd /home/htha/Infineon/radar_sdk/build/external/strata/library/universal && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/htha/Infineon/radar_sdk/external/strata/library/universal/error_strings.c -o CMakeFiles/universal.dir/error_strings.c.s

universal: external/strata/library/universal/CMakeFiles/universal.dir/error_strings.c.o
universal: external/strata/library/universal/CMakeFiles/universal.dir/build.make

.PHONY : universal

# Rule to build all files generated by this target.
external/strata/library/universal/CMakeFiles/universal.dir/build: universal

.PHONY : external/strata/library/universal/CMakeFiles/universal.dir/build

external/strata/library/universal/CMakeFiles/universal.dir/clean:
	cd /home/htha/Infineon/radar_sdk/build/external/strata/library/universal && $(CMAKE_COMMAND) -P CMakeFiles/universal.dir/cmake_clean.cmake
.PHONY : external/strata/library/universal/CMakeFiles/universal.dir/clean

external/strata/library/universal/CMakeFiles/universal.dir/depend:
	cd /home/htha/Infineon/radar_sdk/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/htha/Infineon/radar_sdk /home/htha/Infineon/radar_sdk/external/strata/library/universal /home/htha/Infineon/radar_sdk/build /home/htha/Infineon/radar_sdk/build/external/strata/library/universal /home/htha/Infineon/radar_sdk/build/external/strata/library/universal/CMakeFiles/universal.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/strata/library/universal/CMakeFiles/universal.dir/depend

