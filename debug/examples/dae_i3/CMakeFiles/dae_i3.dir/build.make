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
include examples/dae_i3/CMakeFiles/dae_i3.dir/depend.make

# Include the progress variables for this target.
include examples/dae_i3/CMakeFiles/dae_i3.dir/progress.make

# Include the compile flags for this target's objects.
include examples/dae_i3/CMakeFiles/dae_i3.dir/flags.make

examples/dae_i3/CMakeFiles/dae_i3.dir/dae_i3.cxx.o: examples/dae_i3/CMakeFiles/dae_i3.dir/flags.make
examples/dae_i3/CMakeFiles/dae_i3.dir/dae_i3.cxx.o: ../examples/dae_i3/dae_i3.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/matt/git/psopt/debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/dae_i3/CMakeFiles/dae_i3.dir/dae_i3.cxx.o"
	cd /home/matt/git/psopt/debug/examples/dae_i3 && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/dae_i3.dir/dae_i3.cxx.o -c /home/matt/git/psopt/examples/dae_i3/dae_i3.cxx

examples/dae_i3/CMakeFiles/dae_i3.dir/dae_i3.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dae_i3.dir/dae_i3.cxx.i"
	cd /home/matt/git/psopt/debug/examples/dae_i3 && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/matt/git/psopt/examples/dae_i3/dae_i3.cxx > CMakeFiles/dae_i3.dir/dae_i3.cxx.i

examples/dae_i3/CMakeFiles/dae_i3.dir/dae_i3.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dae_i3.dir/dae_i3.cxx.s"
	cd /home/matt/git/psopt/debug/examples/dae_i3 && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/matt/git/psopt/examples/dae_i3/dae_i3.cxx -o CMakeFiles/dae_i3.dir/dae_i3.cxx.s

# Object files for target dae_i3
dae_i3_OBJECTS = \
"CMakeFiles/dae_i3.dir/dae_i3.cxx.o"

# External object files for target dae_i3
dae_i3_EXTERNAL_OBJECTS =

examples/dae_i3/dae_i3: examples/dae_i3/CMakeFiles/dae_i3.dir/dae_i3.cxx.o
examples/dae_i3/dae_i3: examples/dae_i3/CMakeFiles/dae_i3.dir/build.make
examples/dae_i3/dae_i3: libPSOPT.a
examples/dae_i3/dae_i3: /usr/lib/libipopt.so
examples/dae_i3/dae_i3: /usr/lib/x86_64-linux-gnu/libdmumps_seq.so
examples/dae_i3/dae_i3: /usr/lib/x86_64-linux-gnu/libblas.so
examples/dae_i3/dae_i3: /usr/lib/x86_64-linux-gnu/liblapack.so
examples/dae_i3/dae_i3: /usr/lib/x86_64-linux-gnu/libdl.so
examples/dae_i3/dae_i3: /usr/lib/gcc/x86_64-linux-gnu/9/libgfortran.so
examples/dae_i3/dae_i3: /usr/lib/gcc/x86_64-linux-gnu/9/libquadmath.so
examples/dae_i3/dae_i3: /usr/lib/x86_64-linux-gnu/libdmumps_seq.so
examples/dae_i3/dae_i3: /usr/lib/x86_64-linux-gnu/libblas.so
examples/dae_i3/dae_i3: /usr/lib/x86_64-linux-gnu/liblapack.so
examples/dae_i3/dae_i3: /usr/lib/x86_64-linux-gnu/libdl.so
examples/dae_i3/dae_i3: /usr/lib/gcc/x86_64-linux-gnu/9/libgfortran.so
examples/dae_i3/dae_i3: /usr/lib/gcc/x86_64-linux-gnu/9/libquadmath.so
examples/dae_i3/dae_i3: /usr/lib/x86_64-linux-gnu/libadolc.so
examples/dae_i3/dae_i3: /usr/lib/x86_64-linux-gnu/libboost_system.so
examples/dae_i3/dae_i3: /usr/lib/x86_64-linux-gnu/libColPack.so
examples/dae_i3/dae_i3: /usr/lib/x86_64-linux-gnu/libm.so
examples/dae_i3/dae_i3: examples/dae_i3/CMakeFiles/dae_i3.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/matt/git/psopt/debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable dae_i3"
	cd /home/matt/git/psopt/debug/examples/dae_i3 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/dae_i3.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/dae_i3/CMakeFiles/dae_i3.dir/build: examples/dae_i3/dae_i3

.PHONY : examples/dae_i3/CMakeFiles/dae_i3.dir/build

examples/dae_i3/CMakeFiles/dae_i3.dir/clean:
	cd /home/matt/git/psopt/debug/examples/dae_i3 && $(CMAKE_COMMAND) -P CMakeFiles/dae_i3.dir/cmake_clean.cmake
.PHONY : examples/dae_i3/CMakeFiles/dae_i3.dir/clean

examples/dae_i3/CMakeFiles/dae_i3.dir/depend:
	cd /home/matt/git/psopt/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/matt/git/psopt /home/matt/git/psopt/examples/dae_i3 /home/matt/git/psopt/debug /home/matt/git/psopt/debug/examples/dae_i3 /home/matt/git/psopt/debug/examples/dae_i3/CMakeFiles/dae_i3.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/dae_i3/CMakeFiles/dae_i3.dir/depend

