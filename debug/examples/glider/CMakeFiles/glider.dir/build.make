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
include examples/glider/CMakeFiles/glider.dir/depend.make

# Include the progress variables for this target.
include examples/glider/CMakeFiles/glider.dir/progress.make

# Include the compile flags for this target's objects.
include examples/glider/CMakeFiles/glider.dir/flags.make

examples/glider/CMakeFiles/glider.dir/glider.cxx.o: examples/glider/CMakeFiles/glider.dir/flags.make
examples/glider/CMakeFiles/glider.dir/glider.cxx.o: ../examples/glider/glider.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/matt/git/psopt/debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/glider/CMakeFiles/glider.dir/glider.cxx.o"
	cd /home/matt/git/psopt/debug/examples/glider && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/glider.dir/glider.cxx.o -c /home/matt/git/psopt/examples/glider/glider.cxx

examples/glider/CMakeFiles/glider.dir/glider.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/glider.dir/glider.cxx.i"
	cd /home/matt/git/psopt/debug/examples/glider && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/matt/git/psopt/examples/glider/glider.cxx > CMakeFiles/glider.dir/glider.cxx.i

examples/glider/CMakeFiles/glider.dir/glider.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/glider.dir/glider.cxx.s"
	cd /home/matt/git/psopt/debug/examples/glider && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/matt/git/psopt/examples/glider/glider.cxx -o CMakeFiles/glider.dir/glider.cxx.s

# Object files for target glider
glider_OBJECTS = \
"CMakeFiles/glider.dir/glider.cxx.o"

# External object files for target glider
glider_EXTERNAL_OBJECTS =

examples/glider/glider: examples/glider/CMakeFiles/glider.dir/glider.cxx.o
examples/glider/glider: examples/glider/CMakeFiles/glider.dir/build.make
examples/glider/glider: libPSOPT.a
examples/glider/glider: /usr/lib/libipopt.so
examples/glider/glider: /usr/lib/x86_64-linux-gnu/libdmumps_seq.so
examples/glider/glider: /usr/lib/x86_64-linux-gnu/libblas.so
examples/glider/glider: /usr/lib/x86_64-linux-gnu/liblapack.so
examples/glider/glider: /usr/lib/x86_64-linux-gnu/libdl.so
examples/glider/glider: /usr/lib/gcc/x86_64-linux-gnu/9/libgfortran.so
examples/glider/glider: /usr/lib/gcc/x86_64-linux-gnu/9/libquadmath.so
examples/glider/glider: /usr/lib/x86_64-linux-gnu/libdmumps_seq.so
examples/glider/glider: /usr/lib/x86_64-linux-gnu/libblas.so
examples/glider/glider: /usr/lib/x86_64-linux-gnu/liblapack.so
examples/glider/glider: /usr/lib/x86_64-linux-gnu/libdl.so
examples/glider/glider: /usr/lib/gcc/x86_64-linux-gnu/9/libgfortran.so
examples/glider/glider: /usr/lib/gcc/x86_64-linux-gnu/9/libquadmath.so
examples/glider/glider: /usr/lib/x86_64-linux-gnu/libadolc.so
examples/glider/glider: /usr/lib/x86_64-linux-gnu/libboost_system.so
examples/glider/glider: /usr/lib/x86_64-linux-gnu/libColPack.so
examples/glider/glider: /usr/lib/x86_64-linux-gnu/libm.so
examples/glider/glider: examples/glider/CMakeFiles/glider.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/matt/git/psopt/debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable glider"
	cd /home/matt/git/psopt/debug/examples/glider && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/glider.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/glider/CMakeFiles/glider.dir/build: examples/glider/glider

.PHONY : examples/glider/CMakeFiles/glider.dir/build

examples/glider/CMakeFiles/glider.dir/clean:
	cd /home/matt/git/psopt/debug/examples/glider && $(CMAKE_COMMAND) -P CMakeFiles/glider.dir/cmake_clean.cmake
.PHONY : examples/glider/CMakeFiles/glider.dir/clean

examples/glider/CMakeFiles/glider.dir/depend:
	cd /home/matt/git/psopt/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/matt/git/psopt /home/matt/git/psopt/examples/glider /home/matt/git/psopt/debug /home/matt/git/psopt/debug/examples/glider /home/matt/git/psopt/debug/examples/glider/CMakeFiles/glider.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/glider/CMakeFiles/glider.dir/depend

