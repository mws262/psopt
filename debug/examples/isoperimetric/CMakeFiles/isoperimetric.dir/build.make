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
CMAKE_SOURCE_DIR = /home/matt/git/psopt

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/matt/git/psopt/debug

# Include any dependencies generated for this target.
include examples/isoperimetric/CMakeFiles/isoperimetric.dir/depend.make

# Include the progress variables for this target.
include examples/isoperimetric/CMakeFiles/isoperimetric.dir/progress.make

# Include the compile flags for this target's objects.
include examples/isoperimetric/CMakeFiles/isoperimetric.dir/flags.make

examples/isoperimetric/CMakeFiles/isoperimetric.dir/isoperimetric.cxx.o: examples/isoperimetric/CMakeFiles/isoperimetric.dir/flags.make
examples/isoperimetric/CMakeFiles/isoperimetric.dir/isoperimetric.cxx.o: ../examples/isoperimetric/isoperimetric.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/matt/git/psopt/debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/isoperimetric/CMakeFiles/isoperimetric.dir/isoperimetric.cxx.o"
	cd /home/matt/git/psopt/debug/examples/isoperimetric && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/isoperimetric.dir/isoperimetric.cxx.o -c /home/matt/git/psopt/examples/isoperimetric/isoperimetric.cxx

examples/isoperimetric/CMakeFiles/isoperimetric.dir/isoperimetric.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/isoperimetric.dir/isoperimetric.cxx.i"
	cd /home/matt/git/psopt/debug/examples/isoperimetric && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/matt/git/psopt/examples/isoperimetric/isoperimetric.cxx > CMakeFiles/isoperimetric.dir/isoperimetric.cxx.i

examples/isoperimetric/CMakeFiles/isoperimetric.dir/isoperimetric.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/isoperimetric.dir/isoperimetric.cxx.s"
	cd /home/matt/git/psopt/debug/examples/isoperimetric && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/matt/git/psopt/examples/isoperimetric/isoperimetric.cxx -o CMakeFiles/isoperimetric.dir/isoperimetric.cxx.s

# Object files for target isoperimetric
isoperimetric_OBJECTS = \
"CMakeFiles/isoperimetric.dir/isoperimetric.cxx.o"

# External object files for target isoperimetric
isoperimetric_EXTERNAL_OBJECTS =

examples/isoperimetric/isoperimetric: examples/isoperimetric/CMakeFiles/isoperimetric.dir/isoperimetric.cxx.o
examples/isoperimetric/isoperimetric: examples/isoperimetric/CMakeFiles/isoperimetric.dir/build.make
examples/isoperimetric/isoperimetric: libPSOPT.a
examples/isoperimetric/isoperimetric: /usr/lib/libipopt.so
examples/isoperimetric/isoperimetric: /usr/lib/x86_64-linux-gnu/libdmumps_seq.so
examples/isoperimetric/isoperimetric: /usr/lib/x86_64-linux-gnu/libblas.so
examples/isoperimetric/isoperimetric: /usr/lib/x86_64-linux-gnu/liblapack.so
examples/isoperimetric/isoperimetric: /usr/lib/x86_64-linux-gnu/libdl.so
examples/isoperimetric/isoperimetric: /usr/lib/gcc/x86_64-linux-gnu/9/libgfortran.so
examples/isoperimetric/isoperimetric: /usr/lib/gcc/x86_64-linux-gnu/9/libquadmath.so
examples/isoperimetric/isoperimetric: /usr/lib/x86_64-linux-gnu/libdmumps_seq.so
examples/isoperimetric/isoperimetric: /usr/lib/x86_64-linux-gnu/libblas.so
examples/isoperimetric/isoperimetric: /usr/lib/x86_64-linux-gnu/liblapack.so
examples/isoperimetric/isoperimetric: /usr/lib/x86_64-linux-gnu/libdl.so
examples/isoperimetric/isoperimetric: /usr/lib/gcc/x86_64-linux-gnu/9/libgfortran.so
examples/isoperimetric/isoperimetric: /usr/lib/gcc/x86_64-linux-gnu/9/libquadmath.so
examples/isoperimetric/isoperimetric: /usr/lib/x86_64-linux-gnu/libadolc.so
examples/isoperimetric/isoperimetric: /usr/lib/x86_64-linux-gnu/libboost_system.so
examples/isoperimetric/isoperimetric: /usr/lib/x86_64-linux-gnu/libColPack.so
examples/isoperimetric/isoperimetric: /usr/lib/x86_64-linux-gnu/libm.so
examples/isoperimetric/isoperimetric: examples/isoperimetric/CMakeFiles/isoperimetric.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/matt/git/psopt/debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable isoperimetric"
	cd /home/matt/git/psopt/debug/examples/isoperimetric && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/isoperimetric.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/isoperimetric/CMakeFiles/isoperimetric.dir/build: examples/isoperimetric/isoperimetric

.PHONY : examples/isoperimetric/CMakeFiles/isoperimetric.dir/build

examples/isoperimetric/CMakeFiles/isoperimetric.dir/clean:
	cd /home/matt/git/psopt/debug/examples/isoperimetric && $(CMAKE_COMMAND) -P CMakeFiles/isoperimetric.dir/cmake_clean.cmake
.PHONY : examples/isoperimetric/CMakeFiles/isoperimetric.dir/clean

examples/isoperimetric/CMakeFiles/isoperimetric.dir/depend:
	cd /home/matt/git/psopt/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/matt/git/psopt /home/matt/git/psopt/examples/isoperimetric /home/matt/git/psopt/debug /home/matt/git/psopt/debug/examples/isoperimetric /home/matt/git/psopt/debug/examples/isoperimetric/CMakeFiles/isoperimetric.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/isoperimetric/CMakeFiles/isoperimetric.dir/depend

