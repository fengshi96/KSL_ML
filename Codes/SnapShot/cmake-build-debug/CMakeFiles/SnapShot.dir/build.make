# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/shifeng/IDE/clion-2020.2.5/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/shifeng/IDE/clion-2020.2.5/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Barn/Projects/7.Kitaev_ML/Codes/SnapShot

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Barn/Projects/7.Kitaev_ML/Codes/SnapShot/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/SnapShot.dir/depend.make
# Include the progress variables for this target.
include CMakeFiles/SnapShot.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/SnapShot.dir/flags.make

CMakeFiles/SnapShot.dir/main.cpp.o: CMakeFiles/SnapShot.dir/flags.make
CMakeFiles/SnapShot.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Barn/Projects/7.Kitaev_ML/Codes/SnapShot/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/SnapShot.dir/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/SnapShot.dir/main.cpp.o -c /Barn/Projects/7.Kitaev_ML/Codes/SnapShot/main.cpp

CMakeFiles/SnapShot.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SnapShot.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Barn/Projects/7.Kitaev_ML/Codes/SnapShot/main.cpp > CMakeFiles/SnapShot.dir/main.cpp.i

CMakeFiles/SnapShot.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SnapShot.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Barn/Projects/7.Kitaev_ML/Codes/SnapShot/main.cpp -o CMakeFiles/SnapShot.dir/main.cpp.s

# Object files for target SnapShot
SnapShot_OBJECTS = \
"CMakeFiles/SnapShot.dir/main.cpp.o"

# External object files for target SnapShot
SnapShot_EXTERNAL_OBJECTS =

SnapShot: CMakeFiles/SnapShot.dir/main.cpp.o
SnapShot: CMakeFiles/SnapShot.dir/build.make
SnapShot: CMakeFiles/SnapShot.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Barn/Projects/7.Kitaev_ML/Codes/SnapShot/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable SnapShot"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/SnapShot.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/SnapShot.dir/build: SnapShot
.PHONY : CMakeFiles/SnapShot.dir/build

CMakeFiles/SnapShot.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/SnapShot.dir/cmake_clean.cmake
.PHONY : CMakeFiles/SnapShot.dir/clean

CMakeFiles/SnapShot.dir/depend:
	cd /Barn/Projects/7.Kitaev_ML/Codes/SnapShot/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Barn/Projects/7.Kitaev_ML/Codes/SnapShot /Barn/Projects/7.Kitaev_ML/Codes/SnapShot /Barn/Projects/7.Kitaev_ML/Codes/SnapShot/cmake-build-debug /Barn/Projects/7.Kitaev_ML/Codes/SnapShot/cmake-build-debug /Barn/Projects/7.Kitaev_ML/Codes/SnapShot/cmake-build-debug/CMakeFiles/SnapShot.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/SnapShot.dir/depend

