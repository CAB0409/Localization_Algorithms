# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_SOURCE_DIR = /home/chris/software/project_x/workspace

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/chris/software/project_x/workspace/build

# Include any dependencies generated for this target.
include CMakeFiles/tdoa.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/tdoa.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/tdoa.dir/flags.make

CMakeFiles/tdoa.dir/lib/main.cpp.o: CMakeFiles/tdoa.dir/flags.make
CMakeFiles/tdoa.dir/lib/main.cpp.o: ../lib/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chris/software/project_x/workspace/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/tdoa.dir/lib/main.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/tdoa.dir/lib/main.cpp.o -c /home/chris/software/project_x/workspace/lib/main.cpp

CMakeFiles/tdoa.dir/lib/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tdoa.dir/lib/main.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chris/software/project_x/workspace/lib/main.cpp > CMakeFiles/tdoa.dir/lib/main.cpp.i

CMakeFiles/tdoa.dir/lib/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tdoa.dir/lib/main.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chris/software/project_x/workspace/lib/main.cpp -o CMakeFiles/tdoa.dir/lib/main.cpp.s

CMakeFiles/tdoa.dir/lib/main.cpp.o.requires:

.PHONY : CMakeFiles/tdoa.dir/lib/main.cpp.o.requires

CMakeFiles/tdoa.dir/lib/main.cpp.o.provides: CMakeFiles/tdoa.dir/lib/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/tdoa.dir/build.make CMakeFiles/tdoa.dir/lib/main.cpp.o.provides.build
.PHONY : CMakeFiles/tdoa.dir/lib/main.cpp.o.provides

CMakeFiles/tdoa.dir/lib/main.cpp.o.provides.build: CMakeFiles/tdoa.dir/lib/main.cpp.o


# Object files for target tdoa
tdoa_OBJECTS = \
"CMakeFiles/tdoa.dir/lib/main.cpp.o"

# External object files for target tdoa
tdoa_EXTERNAL_OBJECTS =

tdoa: CMakeFiles/tdoa.dir/lib/main.cpp.o
tdoa: CMakeFiles/tdoa.dir/build.make
tdoa: /usr/lib/x86_64-linux-gnu/libarmadillo.so
tdoa: /home/chris/software/2018matlab/bin/glnxa64/libmex.so
tdoa: CMakeFiles/tdoa.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/chris/software/project_x/workspace/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable tdoa"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tdoa.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/tdoa.dir/build: tdoa

.PHONY : CMakeFiles/tdoa.dir/build

CMakeFiles/tdoa.dir/requires: CMakeFiles/tdoa.dir/lib/main.cpp.o.requires

.PHONY : CMakeFiles/tdoa.dir/requires

CMakeFiles/tdoa.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/tdoa.dir/cmake_clean.cmake
.PHONY : CMakeFiles/tdoa.dir/clean

CMakeFiles/tdoa.dir/depend:
	cd /home/chris/software/project_x/workspace/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chris/software/project_x/workspace /home/chris/software/project_x/workspace /home/chris/software/project_x/workspace/build /home/chris/software/project_x/workspace/build /home/chris/software/project_x/workspace/build/CMakeFiles/tdoa.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/tdoa.dir/depend
