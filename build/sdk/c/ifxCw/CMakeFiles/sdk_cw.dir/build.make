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
include sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/depend.make

# Include the progress variables for this target.
include sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/progress.make

# Include the compile flags for this target's objects.
include sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/flags.make

sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/DeviceCwTypes.c.o: sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/flags.make
sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/DeviceCwTypes.c.o: ../sdk/c/ifxCw/DeviceCwTypes.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/htha/Infineon/radar_sdk/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/DeviceCwTypes.c.o"
	cd /home/htha/Infineon/radar_sdk/build/sdk/c/ifxCw && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sdk_cw.dir/DeviceCwTypes.c.o   -c /home/htha/Infineon/radar_sdk/sdk/c/ifxCw/DeviceCwTypes.c

sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/DeviceCwTypes.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sdk_cw.dir/DeviceCwTypes.c.i"
	cd /home/htha/Infineon/radar_sdk/build/sdk/c/ifxCw && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/htha/Infineon/radar_sdk/sdk/c/ifxCw/DeviceCwTypes.c > CMakeFiles/sdk_cw.dir/DeviceCwTypes.c.i

sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/DeviceCwTypes.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sdk_cw.dir/DeviceCwTypes.c.s"
	cd /home/htha/Infineon/radar_sdk/build/sdk/c/ifxCw && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/htha/Infineon/radar_sdk/sdk/c/ifxCw/DeviceCwTypes.c -o CMakeFiles/sdk_cw.dir/DeviceCwTypes.c.s

sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/DeviceCw.cpp.o: sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/flags.make
sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/DeviceCw.cpp.o: ../sdk/c/ifxCw/DeviceCw.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/htha/Infineon/radar_sdk/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/DeviceCw.cpp.o"
	cd /home/htha/Infineon/radar_sdk/build/sdk/c/ifxCw && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/sdk_cw.dir/DeviceCw.cpp.o -c /home/htha/Infineon/radar_sdk/sdk/c/ifxCw/DeviceCw.cpp

sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/DeviceCw.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/sdk_cw.dir/DeviceCw.cpp.i"
	cd /home/htha/Infineon/radar_sdk/build/sdk/c/ifxCw && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/htha/Infineon/radar_sdk/sdk/c/ifxCw/DeviceCw.cpp > CMakeFiles/sdk_cw.dir/DeviceCw.cpp.i

sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/DeviceCw.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/sdk_cw.dir/DeviceCw.cpp.s"
	cd /home/htha/Infineon/radar_sdk/build/sdk/c/ifxCw && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/htha/Infineon/radar_sdk/sdk/c/ifxCw/DeviceCw.cpp -o CMakeFiles/sdk_cw.dir/DeviceCw.cpp.s

sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/DeviceCwBase.cpp.o: sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/flags.make
sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/DeviceCwBase.cpp.o: ../sdk/c/ifxCw/DeviceCwBase.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/htha/Infineon/radar_sdk/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/DeviceCwBase.cpp.o"
	cd /home/htha/Infineon/radar_sdk/build/sdk/c/ifxCw && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/sdk_cw.dir/DeviceCwBase.cpp.o -c /home/htha/Infineon/radar_sdk/sdk/c/ifxCw/DeviceCwBase.cpp

sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/DeviceCwBase.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/sdk_cw.dir/DeviceCwBase.cpp.i"
	cd /home/htha/Infineon/radar_sdk/build/sdk/c/ifxCw && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/htha/Infineon/radar_sdk/sdk/c/ifxCw/DeviceCwBase.cpp > CMakeFiles/sdk_cw.dir/DeviceCwBase.cpp.i

sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/DeviceCwBase.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/sdk_cw.dir/DeviceCwBase.cpp.s"
	cd /home/htha/Infineon/radar_sdk/build/sdk/c/ifxCw && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/htha/Infineon/radar_sdk/sdk/c/ifxCw/DeviceCwBase.cpp -o CMakeFiles/sdk_cw.dir/DeviceCwBase.cpp.s

sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/avian/DeviceCwAvian.cpp.o: sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/flags.make
sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/avian/DeviceCwAvian.cpp.o: ../sdk/c/ifxCw/avian/DeviceCwAvian.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/htha/Infineon/radar_sdk/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/avian/DeviceCwAvian.cpp.o"
	cd /home/htha/Infineon/radar_sdk/build/sdk/c/ifxCw && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/sdk_cw.dir/avian/DeviceCwAvian.cpp.o -c /home/htha/Infineon/radar_sdk/sdk/c/ifxCw/avian/DeviceCwAvian.cpp

sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/avian/DeviceCwAvian.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/sdk_cw.dir/avian/DeviceCwAvian.cpp.i"
	cd /home/htha/Infineon/radar_sdk/build/sdk/c/ifxCw && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/htha/Infineon/radar_sdk/sdk/c/ifxCw/avian/DeviceCwAvian.cpp > CMakeFiles/sdk_cw.dir/avian/DeviceCwAvian.cpp.i

sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/avian/DeviceCwAvian.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/sdk_cw.dir/avian/DeviceCwAvian.cpp.s"
	cd /home/htha/Infineon/radar_sdk/build/sdk/c/ifxCw && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/htha/Infineon/radar_sdk/sdk/c/ifxCw/avian/DeviceCwAvian.cpp -o CMakeFiles/sdk_cw.dir/avian/DeviceCwAvian.cpp.s

# Object files for target sdk_cw
sdk_cw_OBJECTS = \
"CMakeFiles/sdk_cw.dir/DeviceCwTypes.c.o" \
"CMakeFiles/sdk_cw.dir/DeviceCw.cpp.o" \
"CMakeFiles/sdk_cw.dir/DeviceCwBase.cpp.o" \
"CMakeFiles/sdk_cw.dir/avian/DeviceCwAvian.cpp.o"

# External object files for target sdk_cw
sdk_cw_EXTERNAL_OBJECTS =

bin/libsdk_cw.so: sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/DeviceCwTypes.c.o
bin/libsdk_cw.so: sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/DeviceCw.cpp.o
bin/libsdk_cw.so: sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/DeviceCwBase.cpp.o
bin/libsdk_cw.so: sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/avian/DeviceCwAvian.cpp.o
bin/libsdk_cw.so: sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/build.make
bin/libsdk_cw.so: bin/libsdk_radar_device_common.so
bin/libsdk_cw.so: bin/liblib_avian.so
bin/libsdk_cw.so: bin/libsdk_base.so
bin/libsdk_cw.so: bin/libstrata_shared.so
bin/libsdk_cw.so: external/strata/contrib/pugixml/libpugixml.a
bin/libsdk_cw.so: sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/htha/Infineon/radar_sdk/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX shared library ../../../bin/libsdk_cw.so"
	cd /home/htha/Infineon/radar_sdk/build/sdk/c/ifxCw && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sdk_cw.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/build: bin/libsdk_cw.so

.PHONY : sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/build

sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/clean:
	cd /home/htha/Infineon/radar_sdk/build/sdk/c/ifxCw && $(CMAKE_COMMAND) -P CMakeFiles/sdk_cw.dir/cmake_clean.cmake
.PHONY : sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/clean

sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/depend:
	cd /home/htha/Infineon/radar_sdk/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/htha/Infineon/radar_sdk /home/htha/Infineon/radar_sdk/sdk/c/ifxCw /home/htha/Infineon/radar_sdk/build /home/htha/Infineon/radar_sdk/build/sdk/c/ifxCw /home/htha/Infineon/radar_sdk/build/sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : sdk/c/ifxCw/CMakeFiles/sdk_cw.dir/depend

