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
include sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/depend.make

# Include the progress variables for this target.
include sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/progress.make

# Include the compile flags for this target's objects.
include sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/flags.make

sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/COMPort_Unix.c.o: sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/flags.make
sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/COMPort_Unix.c.o: ../sdk/c/ifxComPort/COMPort_Unix.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/htha/Infineon/radar_sdk/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/COMPort_Unix.c.o"
	cd /home/htha/Infineon/radar_sdk/build/sdk/c/ifxComPort && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sdk_comport_obj.dir/COMPort_Unix.c.o   -c /home/htha/Infineon/radar_sdk/sdk/c/ifxComPort/COMPort_Unix.c

sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/COMPort_Unix.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sdk_comport_obj.dir/COMPort_Unix.c.i"
	cd /home/htha/Infineon/radar_sdk/build/sdk/c/ifxComPort && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/htha/Infineon/radar_sdk/sdk/c/ifxComPort/COMPort_Unix.c > CMakeFiles/sdk_comport_obj.dir/COMPort_Unix.c.i

sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/COMPort_Unix.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sdk_comport_obj.dir/COMPort_Unix.c.s"
	cd /home/htha/Infineon/radar_sdk/build/sdk/c/ifxComPort && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/htha/Infineon/radar_sdk/sdk/c/ifxComPort/COMPort_Unix.c -o CMakeFiles/sdk_comport_obj.dir/COMPort_Unix.c.s

sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/COMPort_Windows.c.o: sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/flags.make
sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/COMPort_Windows.c.o: ../sdk/c/ifxComPort/COMPort_Windows.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/htha/Infineon/radar_sdk/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/COMPort_Windows.c.o"
	cd /home/htha/Infineon/radar_sdk/build/sdk/c/ifxComPort && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sdk_comport_obj.dir/COMPort_Windows.c.o   -c /home/htha/Infineon/radar_sdk/sdk/c/ifxComPort/COMPort_Windows.c

sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/COMPort_Windows.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sdk_comport_obj.dir/COMPort_Windows.c.i"
	cd /home/htha/Infineon/radar_sdk/build/sdk/c/ifxComPort && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/htha/Infineon/radar_sdk/sdk/c/ifxComPort/COMPort_Windows.c > CMakeFiles/sdk_comport_obj.dir/COMPort_Windows.c.i

sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/COMPort_Windows.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sdk_comport_obj.dir/COMPort_Windows.c.s"
	cd /home/htha/Infineon/radar_sdk/build/sdk/c/ifxComPort && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/htha/Infineon/radar_sdk/sdk/c/ifxComPort/COMPort_Windows.c -o CMakeFiles/sdk_comport_obj.dir/COMPort_Windows.c.s

sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/ComPort.cpp.o: sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/flags.make
sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/ComPort.cpp.o: ../sdk/c/ifxComPort/ComPort.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/htha/Infineon/radar_sdk/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/ComPort.cpp.o"
	cd /home/htha/Infineon/radar_sdk/build/sdk/c/ifxComPort && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/sdk_comport_obj.dir/ComPort.cpp.o -c /home/htha/Infineon/radar_sdk/sdk/c/ifxComPort/ComPort.cpp

sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/ComPort.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/sdk_comport_obj.dir/ComPort.cpp.i"
	cd /home/htha/Infineon/radar_sdk/build/sdk/c/ifxComPort && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/htha/Infineon/radar_sdk/sdk/c/ifxComPort/ComPort.cpp > CMakeFiles/sdk_comport_obj.dir/ComPort.cpp.i

sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/ComPort.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/sdk_comport_obj.dir/ComPort.cpp.s"
	cd /home/htha/Infineon/radar_sdk/build/sdk/c/ifxComPort && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/htha/Infineon/radar_sdk/sdk/c/ifxComPort/ComPort.cpp -o CMakeFiles/sdk_comport_obj.dir/ComPort.cpp.s

sdk_comport_obj: sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/COMPort_Unix.c.o
sdk_comport_obj: sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/COMPort_Windows.c.o
sdk_comport_obj: sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/ComPort.cpp.o
sdk_comport_obj: sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/build.make

.PHONY : sdk_comport_obj

# Rule to build all files generated by this target.
sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/build: sdk_comport_obj

.PHONY : sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/build

sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/clean:
	cd /home/htha/Infineon/radar_sdk/build/sdk/c/ifxComPort && $(CMAKE_COMMAND) -P CMakeFiles/sdk_comport_obj.dir/cmake_clean.cmake
.PHONY : sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/clean

sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/depend:
	cd /home/htha/Infineon/radar_sdk/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/htha/Infineon/radar_sdk /home/htha/Infineon/radar_sdk/sdk/c/ifxComPort /home/htha/Infineon/radar_sdk/build /home/htha/Infineon/radar_sdk/build/sdk/c/ifxComPort /home/htha/Infineon/radar_sdk/build/sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : sdk/c/ifxComPort/CMakeFiles/sdk_comport_obj.dir/depend

