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
include external/strata/library/CMakeFiles/strata_static.dir/depend.make

# Include the progress variables for this target.
include external/strata/library/CMakeFiles/strata_static.dir/progress.make

# Include the compile flags for this target's objects.
include external/strata/library/CMakeFiles/strata_static.dir/flags.make

external/strata/library/CMakeFiles/strata_static.dir/Dummy.cpp.o: external/strata/library/CMakeFiles/strata_static.dir/flags.make
external/strata/library/CMakeFiles/strata_static.dir/Dummy.cpp.o: ../external/strata/library/Dummy.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/htha/Infineon/radar_sdk/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object external/strata/library/CMakeFiles/strata_static.dir/Dummy.cpp.o"
	cd /home/htha/Infineon/radar_sdk/build/external/strata/library && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/strata_static.dir/Dummy.cpp.o -c /home/htha/Infineon/radar_sdk/external/strata/library/Dummy.cpp

external/strata/library/CMakeFiles/strata_static.dir/Dummy.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/strata_static.dir/Dummy.cpp.i"
	cd /home/htha/Infineon/radar_sdk/build/external/strata/library && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/htha/Infineon/radar_sdk/external/strata/library/Dummy.cpp > CMakeFiles/strata_static.dir/Dummy.cpp.i

external/strata/library/CMakeFiles/strata_static.dir/Dummy.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/strata_static.dir/Dummy.cpp.s"
	cd /home/htha/Infineon/radar_sdk/build/external/strata/library && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/htha/Infineon/radar_sdk/external/strata/library/Dummy.cpp -o CMakeFiles/strata_static.dir/Dummy.cpp.s

# Object files for target strata_static
strata_static_OBJECTS = \
"CMakeFiles/strata_static.dir/Dummy.cpp.o"

# External object files for target strata_static
strata_static_EXTERNAL_OBJECTS = \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/CMakeFiles/strata_boards.dir/BoardList.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/CMakeFiles/strata_version.dir/Version.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/CMakeFiles/strata_library.dir/Library.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/common/CMakeFiles/common.dir/crc/Crc6.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/common/CMakeFiles/common.dir/crc/Crc8.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/common/CMakeFiles/common.dir/crc/Crc16.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/common/CMakeFiles/common.dir/crc/Crc32.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/common/CMakeFiles/common.dir/endian/LittleEndianReader.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/common/CMakeFiles/common.dir/Logger.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/common/CMakeFiles/common.dir/Profiler.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/common/CMakeFiles/common.dir/ProductVersion.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/components/CMakeFiles/components.dir/FlashImage.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/components/CMakeFiles/components.dir/Registers8bitPec.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/components/CMakeFiles/components.dir/imager/ImagerIrs.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/components/CMakeFiles/components.dir/imager/ImagerIrs11x5.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/components/CMakeFiles/components.dir/imager/ImagerIrs16x5.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/components/CMakeFiles/components.dir/imager/PinsIrs.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/components/CMakeFiles/components.dir/imager/RegistersIrs.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/components/CMakeFiles/components.dir/powerSupply/PowerSupplyMax2043xPec.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/components/CMakeFiles/components.dir/powerSupply/SupplyMonitorIna231.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/components/CMakeFiles/components.dir/radar/TypeSerialization.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/components/CMakeFiles/components.dir/nonvolatileMemory/NonvolatileMemory.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/components/CMakeFiles/components.dir/nonvolatileMemory/NonvolatileMemoryEepromI2c.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/components/CMakeFiles/components.dir/nonvolatileMemory/NonvolatileMemoryFlash.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/components/CMakeFiles/components.dir/nonvolatileMemory/NonvolatileMemoryFlashSpi.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/components/CMakeFiles/components.dir/temperature/TemperatureSensorTMP102.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/components/CMakeFiles/components.dir/temperature/TemperatureSensorMCP98x43.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/components/CMakeFiles/components.dir/processing/ProcessingRadar.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/BoardAny.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/BoardDescriptor.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/BoardInstance.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/BoardManager.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/NamedMemory.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/Memory.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/RegisterGenerator.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/bridge/BridgeControl.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/bridge/BridgeData.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/bridge/BridgeProtocol.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/bridge/BridgeProtocolI2c.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/bridge/BridgeProtocolGpio.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/bridge/BridgeProtocolSpi.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/bridge/BridgeProtocolMemory.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/bridge/BridgeProtocolFlash.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/bridge/BridgeProtocolData.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/bridge/BridgeWrapperBase.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/bridge/VendorCommandsImpl.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/ethernet/BridgeEthernet.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/ethernet/BridgeEthernetControl.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/ethernet/BridgeEthernetData.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/ethernet/BridgeEthernetTcp.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/ethernet/BridgeEthernetUdp.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/ethernet/EnumeratorEthernet.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/ethernet/BoardEthernetTcp.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/ethernet/BoardEthernetUdp.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/ethernet/BoardDescriptorEthernet.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/ethernet/SocketTcp.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/ethernet/SocketUdp.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/fpga/BridgeFpgaIrpli.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/frames/DebugFrame.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/frames/ErrorFrame.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/frames/Frame.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/frames/FrameBase.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/frames/FrameForwarder.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/frames/FramePool.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/frames/FrameQueue.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/frames/FrameHelper.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/macro/BoardInstanceMacro.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/macro/BridgeMacro.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/serial/SerialPort.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/serial/BridgeSerial.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/serial/BoardSerial.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/uvc/BoardUvc.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/uvc/VendorExtensionCypress.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/uvc/VendorExtensionRealtek.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/uvc/VendorExtensionRealtekFlash.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/uvc/VendorExtensionRealtekI2c.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/boards/Board.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/boards/BoardGeneric.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/boards/BoardRemote.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/impl/Linux/serial/SerialPortImpl.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/impl/Linux/serial/EnumeratorSerialImpl.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/impl/Linux/video/EnumeratorV4l2.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/impl/Linux/video/BridgeV4l2.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/impl/Linux/video/FrameV4l2.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/impl/Linux/video/FramePoolV4l2.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/impl/unix/ethernet/SocketImpl.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/impl/unix/ethernet/SocketTcpImpl.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/impl/unix/ethernet/SocketUdpImpl.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/impl/unix/serial/SerialPortImplBase.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/platform/CMakeFiles/platform.dir/impl/unix/serial/EnumeratorSerialImplBase.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/remote/CMakeFiles/remote.dir/RemoteGasBoyle.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/remote/CMakeFiles/remote.dir/RemoteProtocolAtr22.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/remote/CMakeFiles/remote.dir/RemoteProtocolAvian.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/remote/CMakeFiles/remote.dir/RemoteProtocolLtr11.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/remote/CMakeFiles/remote.dir/RemoteProtocolSmartar.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/remote/CMakeFiles/remote.dir/RemotePinsAvian.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/remote/CMakeFiles/remote.dir/RemotePinsLtr11.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/remote/CMakeFiles/remote.dir/RemotePinsSmartar.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/remote/CMakeFiles/remote.dir/RemoteProcessingRadar.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/remote/CMakeFiles/remote.dir/RemoteRegisters.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/remote/CMakeFiles/remote.dir/RemoteRadarAvian.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/remote/CMakeFiles/remote.dir/RemoteRadarAtr22.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/remote/CMakeFiles/remote.dir/RemoteRadarLtr11.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/remote/CMakeFiles/remote.dir/RemoteRadarSmartar.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/remote/CMakeFiles/remote.dir/RemoteVendorCommands.cpp.o" \
"/home/htha/Infineon/radar_sdk/build/external/strata/library/universal/CMakeFiles/universal.dir/error_strings.c.o"

external/strata/library/libstrata_static.a: external/strata/library/CMakeFiles/strata_static.dir/Dummy.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/CMakeFiles/strata_boards.dir/BoardList.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/CMakeFiles/strata_version.dir/Version.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/CMakeFiles/strata_library.dir/Library.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/common/CMakeFiles/common.dir/crc/Crc6.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/common/CMakeFiles/common.dir/crc/Crc8.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/common/CMakeFiles/common.dir/crc/Crc16.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/common/CMakeFiles/common.dir/crc/Crc32.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/common/CMakeFiles/common.dir/endian/LittleEndianReader.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/common/CMakeFiles/common.dir/Logger.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/common/CMakeFiles/common.dir/Profiler.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/common/CMakeFiles/common.dir/ProductVersion.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/components/CMakeFiles/components.dir/FlashImage.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/components/CMakeFiles/components.dir/Registers8bitPec.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/components/CMakeFiles/components.dir/imager/ImagerIrs.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/components/CMakeFiles/components.dir/imager/ImagerIrs11x5.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/components/CMakeFiles/components.dir/imager/ImagerIrs16x5.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/components/CMakeFiles/components.dir/imager/PinsIrs.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/components/CMakeFiles/components.dir/imager/RegistersIrs.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/components/CMakeFiles/components.dir/powerSupply/PowerSupplyMax2043xPec.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/components/CMakeFiles/components.dir/powerSupply/SupplyMonitorIna231.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/components/CMakeFiles/components.dir/radar/TypeSerialization.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/components/CMakeFiles/components.dir/nonvolatileMemory/NonvolatileMemory.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/components/CMakeFiles/components.dir/nonvolatileMemory/NonvolatileMemoryEepromI2c.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/components/CMakeFiles/components.dir/nonvolatileMemory/NonvolatileMemoryFlash.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/components/CMakeFiles/components.dir/nonvolatileMemory/NonvolatileMemoryFlashSpi.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/components/CMakeFiles/components.dir/temperature/TemperatureSensorTMP102.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/components/CMakeFiles/components.dir/temperature/TemperatureSensorMCP98x43.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/components/CMakeFiles/components.dir/processing/ProcessingRadar.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/BoardAny.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/BoardDescriptor.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/BoardInstance.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/BoardManager.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/NamedMemory.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/Memory.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/RegisterGenerator.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/bridge/BridgeControl.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/bridge/BridgeData.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/bridge/BridgeProtocol.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/bridge/BridgeProtocolI2c.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/bridge/BridgeProtocolGpio.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/bridge/BridgeProtocolSpi.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/bridge/BridgeProtocolMemory.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/bridge/BridgeProtocolFlash.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/bridge/BridgeProtocolData.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/bridge/BridgeWrapperBase.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/bridge/VendorCommandsImpl.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/ethernet/BridgeEthernet.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/ethernet/BridgeEthernetControl.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/ethernet/BridgeEthernetData.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/ethernet/BridgeEthernetTcp.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/ethernet/BridgeEthernetUdp.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/ethernet/EnumeratorEthernet.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/ethernet/BoardEthernetTcp.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/ethernet/BoardEthernetUdp.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/ethernet/BoardDescriptorEthernet.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/ethernet/SocketTcp.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/ethernet/SocketUdp.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/fpga/BridgeFpgaIrpli.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/frames/DebugFrame.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/frames/ErrorFrame.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/frames/Frame.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/frames/FrameBase.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/frames/FrameForwarder.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/frames/FramePool.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/frames/FrameQueue.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/frames/FrameHelper.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/macro/BoardInstanceMacro.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/macro/BridgeMacro.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/serial/SerialPort.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/serial/BridgeSerial.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/serial/BoardSerial.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/uvc/BoardUvc.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/uvc/VendorExtensionCypress.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/uvc/VendorExtensionRealtek.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/uvc/VendorExtensionRealtekFlash.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/uvc/VendorExtensionRealtekI2c.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/boards/Board.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/boards/BoardGeneric.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/boards/BoardRemote.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/impl/Linux/serial/SerialPortImpl.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/impl/Linux/serial/EnumeratorSerialImpl.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/impl/Linux/video/EnumeratorV4l2.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/impl/Linux/video/BridgeV4l2.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/impl/Linux/video/FrameV4l2.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/impl/Linux/video/FramePoolV4l2.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/impl/unix/ethernet/SocketImpl.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/impl/unix/ethernet/SocketTcpImpl.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/impl/unix/ethernet/SocketUdpImpl.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/impl/unix/serial/SerialPortImplBase.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/platform/CMakeFiles/platform.dir/impl/unix/serial/EnumeratorSerialImplBase.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/remote/CMakeFiles/remote.dir/RemoteGasBoyle.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/remote/CMakeFiles/remote.dir/RemoteProtocolAtr22.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/remote/CMakeFiles/remote.dir/RemoteProtocolAvian.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/remote/CMakeFiles/remote.dir/RemoteProtocolLtr11.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/remote/CMakeFiles/remote.dir/RemoteProtocolSmartar.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/remote/CMakeFiles/remote.dir/RemotePinsAvian.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/remote/CMakeFiles/remote.dir/RemotePinsLtr11.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/remote/CMakeFiles/remote.dir/RemotePinsSmartar.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/remote/CMakeFiles/remote.dir/RemoteProcessingRadar.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/remote/CMakeFiles/remote.dir/RemoteRegisters.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/remote/CMakeFiles/remote.dir/RemoteRadarAvian.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/remote/CMakeFiles/remote.dir/RemoteRadarAtr22.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/remote/CMakeFiles/remote.dir/RemoteRadarLtr11.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/remote/CMakeFiles/remote.dir/RemoteRadarSmartar.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/remote/CMakeFiles/remote.dir/RemoteVendorCommands.cpp.o
external/strata/library/libstrata_static.a: external/strata/library/universal/CMakeFiles/universal.dir/error_strings.c.o
external/strata/library/libstrata_static.a: external/strata/library/CMakeFiles/strata_static.dir/build.make
external/strata/library/libstrata_static.a: external/strata/library/CMakeFiles/strata_static.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/htha/Infineon/radar_sdk/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libstrata_static.a"
	cd /home/htha/Infineon/radar_sdk/build/external/strata/library && $(CMAKE_COMMAND) -P CMakeFiles/strata_static.dir/cmake_clean_target.cmake
	cd /home/htha/Infineon/radar_sdk/build/external/strata/library && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/strata_static.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
external/strata/library/CMakeFiles/strata_static.dir/build: external/strata/library/libstrata_static.a

.PHONY : external/strata/library/CMakeFiles/strata_static.dir/build

external/strata/library/CMakeFiles/strata_static.dir/clean:
	cd /home/htha/Infineon/radar_sdk/build/external/strata/library && $(CMAKE_COMMAND) -P CMakeFiles/strata_static.dir/cmake_clean.cmake
.PHONY : external/strata/library/CMakeFiles/strata_static.dir/clean

external/strata/library/CMakeFiles/strata_static.dir/depend:
	cd /home/htha/Infineon/radar_sdk/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/htha/Infineon/radar_sdk /home/htha/Infineon/radar_sdk/external/strata/library /home/htha/Infineon/radar_sdk/build /home/htha/Infineon/radar_sdk/build/external/strata/library /home/htha/Infineon/radar_sdk/build/external/strata/library/CMakeFiles/strata_static.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/strata/library/CMakeFiles/strata_static.dir/depend

