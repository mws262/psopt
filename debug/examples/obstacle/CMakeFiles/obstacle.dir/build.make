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
include examples/obstacle/CMakeFiles/obstacle.dir/depend.make

# Include the progress variables for this target.
include examples/obstacle/CMakeFiles/obstacle.dir/progress.make

# Include the compile flags for this target's objects.
include examples/obstacle/CMakeFiles/obstacle.dir/flags.make

examples/obstacle/CMakeFiles/obstacle.dir/obstacle.cxx.o: examples/obstacle/CMakeFiles/obstacle.dir/flags.make
examples/obstacle/CMakeFiles/obstacle.dir/obstacle.cxx.o: ../examples/obstacle/obstacle.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/matt/git/psopt/debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/obstacle/CMakeFiles/obstacle.dir/obstacle.cxx.o"
	cd /home/matt/git/psopt/debug/examples/obstacle && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/obstacle.dir/obstacle.cxx.o -c /home/matt/git/psopt/examples/obstacle/obstacle.cxx

examples/obstacle/CMakeFiles/obstacle.dir/obstacle.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/obstacle.dir/obstacle.cxx.i"
	cd /home/matt/git/psopt/debug/examples/obstacle && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/matt/git/psopt/examples/obstacle/obstacle.cxx > CMakeFiles/obstacle.dir/obstacle.cxx.i

examples/obstacle/CMakeFiles/obstacle.dir/obstacle.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/obstacle.dir/obstacle.cxx.s"
	cd /home/matt/git/psopt/debug/examples/obstacle && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/matt/git/psopt/examples/obstacle/obstacle.cxx -o CMakeFiles/obstacle.dir/obstacle.cxx.s

# Object files for target obstacle
obstacle_OBJECTS = \
"CMakeFiles/obstacle.dir/obstacle.cxx.o"

# External object files for target obstacle
obstacle_EXTERNAL_OBJECTS =

examples/obstacle/obstacle: examples/obstacle/CMakeFiles/obstacle.dir/obstacle.cxx.o
examples/obstacle/obstacle: examples/obstacle/CMakeFiles/obstacle.dir/build.make
examples/obstacle/obstacle: libPSOPT.a
examples/obstacle/obstacle: /usr/lib/libipopt.so
examples/obstacle/obstacle: /usr/lib/x86_64-linux-gnu/libdmumps_seq.so
examples/obstacle/obstacle: /usr/lib/x86_64-linux-gnu/libblas.so
examples/obstacle/obstacle: /usr/lib/x86_64-linux-gnu/liblapack.so
examples/obstacle/obstacle: /usr/lib/x86_64-linux-gnu/libdl.so
examples/obstacle/obstacle: /usr/lib/gcc/x86_64-linux-gnu/9/libgfortran.so
examples/obstacle/obstacle: /usr/lib/gcc/x86_64-linux-gnu/9/libquadmath.so
examples/obstacle/obstacle: /usr/lib/x86_64-linux-gnu/libdmumps_seq.so
examples/obstacle/obstacle: /usr/lib/x86_64-linux-gnu/libblas.so
examples/obstacle/obstacle: /usr/lib/x86_64-linux-gnu/liblapack.so
examples/obstacle/obstacle: /usr/lib/x86_64-linux-gnu/libdl.so
examples/obstacle/obstacle: /usr/lib/gcc/x86_64-linux-gnu/9/libgfortran.so
examples/obstacle/obstacle: /usr/lib/gcc/x86_64-linux-gnu/9/libquadmath.so
examples/obstacle/obstacle: /usr/lib/x86_64-linux-gnu/libadolc.so
examples/obstacle/obstacle: /usr/lib/x86_64-linux-gnu/libboost_system.so
examples/obstacle/obstacle: /usr/lib/x86_64-linux-gnu/libColPack.so
examples/obstacle/obstacle: /usr/lib/x86_64-linux-gnu/libm.so
examples/obstacle/obstacle: examples/obstacle/CMakeFiles/obstacle.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/matt/git/psopt/debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable obstacle"
	cd /home/matt/git/psopt/debug/examples/obstacle && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/obstacle.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/obstacle/CMakeFiles/obstacle.dir/build: examples/obstacle/obstacle

.PHONY : examples/obstacle/CMakeFiles/obstacle.dir/build

examples/obstacle/CMakeFiles/obstacle.dir/clean:
	cd /home/matt/git/psopt/debug/examples/obstacle && $(CMAKE_COMMAND) -P CMakeFiles/obstacle.dir/cmake_clean.cmake
.PHONY : examples/obstacle/CMakeFiles/obstacle.dir/clean

examples/obstacle/CMakeFiles/obstacle.dir/depend:
	cd /home/matt/git/psopt/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/matt/git/psopt /home/matt/git/psopt/examples/obstacle /home/matt/git/psopt/debug /home/matt/git/psopt/debug/examples/obstacle /home/matt/git/psopt/debug/examples/obstacle/CMakeFiles/obstacle.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/obstacle/CMakeFiles/obstacle.dir/depend

