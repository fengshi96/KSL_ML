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
CMAKE_SOURCE_DIR = /Barn/Projects/7.Kitaev_ML/Codes/Observe

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Barn/Projects/7.Kitaev_ML/Codes/Observe

# Include any dependencies generated for this target.
include CMakeFiles/Observe.dir/depend.make
# Include the progress variables for this target.
include CMakeFiles/Observe.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Observe.dir/flags.make

CMakeFiles/Observe.dir/main.cpp.o: CMakeFiles/Observe.dir/flags.make
CMakeFiles/Observe.dir/main.cpp.o: main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Barn/Projects/7.Kitaev_ML/Codes/Observe/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Observe.dir/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Observe.dir/main.cpp.o -c /Barn/Projects/7.Kitaev_ML/Codes/Observe/main.cpp

CMakeFiles/Observe.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Observe.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Barn/Projects/7.Kitaev_ML/Codes/Observe/main.cpp > CMakeFiles/Observe.dir/main.cpp.i

CMakeFiles/Observe.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Observe.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Barn/Projects/7.Kitaev_ML/Codes/Observe/main.cpp -o CMakeFiles/Observe.dir/main.cpp.s

# Object files for target Observe
Observe_OBJECTS = \
"CMakeFiles/Observe.dir/main.cpp.o"

# External object files for target Observe
Observe_EXTERNAL_OBJECTS =

Observe: CMakeFiles/Observe.dir/main.cpp.o
Observe: CMakeFiles/Observe.dir/build.make
Observe: CMakeFiles/Observe.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Barn/Projects/7.Kitaev_ML/Codes/Observe/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Observe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Observe.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Observe.dir/build: Observe
.PHONY : CMakeFiles/Observe.dir/build

CMakeFiles/Observe.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Observe.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Observe.dir/clean

CMakeFiles/Observe.dir/depend:
	cd /Barn/Projects/7.Kitaev_ML/Codes/Observe && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Barn/Projects/7.Kitaev_ML/Codes/Observe /Barn/Projects/7.Kitaev_ML/Codes/Observe /Barn/Projects/7.Kitaev_ML/Codes/Observe /Barn/Projects/7.Kitaev_ML/Codes/Observe /Barn/Projects/7.Kitaev_ML/Codes/Observe/CMakeFiles/Observe.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Observe.dir/depend

