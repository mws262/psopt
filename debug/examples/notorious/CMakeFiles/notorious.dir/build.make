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
include examples/notorious/CMakeFiles/notorious.dir/depend.make

# Include the progress variables for this target.
include examples/notorious/CMakeFiles/notorious.dir/progress.make

# Include the compile flags for this target's objects.
include examples/notorious/CMakeFiles/notorious.dir/flags.make

examples/notorious/CMakeFiles/notorious.dir/notorious.cxx.o: examples/notorious/CMakeFiles/notorious.dir/flags.make
examples/notorious/CMakeFiles/notorious.dir/notorious.cxx.o: ../examples/notorious/notorious.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/matt/git/psopt/debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/notorious/CMakeFiles/notorious.dir/notorious.cxx.o"
	cd /home/matt/git/psopt/debug/examples/notorious && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/notorious.dir/notorious.cxx.o -c /home/matt/git/psopt/examples/notorious/notorious.cxx

examples/notorious/CMakeFiles/notorious.dir/notorious.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/notorious.dir/notorious.cxx.i"
	cd /home/matt/git/psopt/debug/examples/notorious && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/matt/git/psopt/examples/notorious/notorious.cxx > CMakeFiles/notorious.dir/notorious.cxx.i

examples/notorious/CMakeFiles/notorious.dir/notorious.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/notorious.dir/notorious.cxx.s"
	cd /home/matt/git/psopt/debug/examples/notorious && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/matt/git/psopt/examples/notorious/notorious.cxx -o CMakeFiles/notorious.dir/notorious.cxx.s

# Object files for target notorious
notorious_OBJECTS = \
"CMakeFiles/notorious.dir/notorious.cxx.o"

# External object files for target notorious
notorious_EXTERNAL_OBJECTS =

examples/notorious/notorious: examples/notorious/CMakeFiles/notorious.dir/notorious.cxx.o
examples/notorious/notorious: examples/notorious/CMakeFiles/notorious.dir/build.make
examples/notorious/notorious: libPSOPT.a
examples/notorious/notorious: /usr/lib/libipopt.so
examples/notorious/notorious: /usr/lib/x86_64-linux-gnu/libdmumps_seq.so
examples/notorious/notorious: /usr/lib/x86_64-linux-gnu/libblas.so
examples/notorious/notorious: /usr/lib/x86_64-linux-gnu/liblapack.so
examples/notorious/notorious: /usr/lib/x86_64-linux-gnu/libdl.so
examples/notorious/notorious: /usr/lib/gcc/x86_64-linux-gnu/9/libgfortran.so
examples/notorious/notorious: /usr/lib/gcc/x86_64-linux-gnu/9/libquadmath.so
examples/notorious/notorious: /usr/lib/x86_64-linux-gnu/libdmumps_seq.so
examples/notorious/notorious: /usr/lib/x86_64-linux-gnu/libblas.so
examples/notorious/notorious: /usr/lib/x86_64-linux-gnu/liblapack.so
examples/notorious/notorious: /usr/lib/x86_64-linux-gnu/libdl.so
examples/notorious/notorious: /usr/lib/gcc/x86_64-linux-gnu/9/libgfortran.so
examples/notorious/notorious: /usr/lib/gcc/x86_64-linux-gnu/9/libquadmath.so
examples/notorious/notorious: /usr/lib/x86_64-linux-gnu/libadolc.so
examples/notorious/notorious: /usr/lib/x86_64-linux-gnu/libboost_system.so
examples/notorious/notorious: /usr/lib/x86_64-linux-gnu/libColPack.so
examples/notorious/notorious: /usr/lib/x86_64-linux-gnu/libm.so
examples/notorious/notorious: examples/notorious/CMakeFiles/notorious.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/matt/git/psopt/debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable notorious"
	cd /home/matt/git/psopt/debug/examples/notorious && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/notorious.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/notorious/CMakeFiles/notorious.dir/build: examples/notorious/notorious

.PHONY : examples/notorious/CMakeFiles/notorious.dir/build

examples/notorious/CMakeFiles/notorious.dir/clean:
	cd /home/matt/git/psopt/debug/examples/notorious && $(CMAKE_COMMAND) -P CMakeFiles/notorious.dir/cmake_clean.cmake
.PHONY : examples/notorious/CMakeFiles/notorious.dir/clean

examples/notorious/CMakeFiles/notorious.dir/depend:
	cd /home/matt/git/psopt/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/matt/git/psopt /home/matt/git/psopt/examples/notorious /home/matt/git/psopt/debug /home/matt/git/psopt/debug/examples/notorious /home/matt/git/psopt/debug/examples/notorious/CMakeFiles/notorious.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/notorious/CMakeFiles/notorious.dir/depend

