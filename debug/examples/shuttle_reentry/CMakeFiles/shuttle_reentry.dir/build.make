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
include examples/shuttle_reentry/CMakeFiles/shuttle_reentry.dir/depend.make

# Include the progress variables for this target.
include examples/shuttle_reentry/CMakeFiles/shuttle_reentry.dir/progress.make

# Include the compile flags for this target's objects.
include examples/shuttle_reentry/CMakeFiles/shuttle_reentry.dir/flags.make

examples/shuttle_reentry/CMakeFiles/shuttle_reentry.dir/shuttle_reentry.cxx.o: examples/shuttle_reentry/CMakeFiles/shuttle_reentry.dir/flags.make
examples/shuttle_reentry/CMakeFiles/shuttle_reentry.dir/shuttle_reentry.cxx.o: ../examples/shuttle_reentry/shuttle_reentry.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/matt/git/psopt/debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/shuttle_reentry/CMakeFiles/shuttle_reentry.dir/shuttle_reentry.cxx.o"
	cd /home/matt/git/psopt/debug/examples/shuttle_reentry && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/shuttle_reentry.dir/shuttle_reentry.cxx.o -c /home/matt/git/psopt/examples/shuttle_reentry/shuttle_reentry.cxx

examples/shuttle_reentry/CMakeFiles/shuttle_reentry.dir/shuttle_reentry.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/shuttle_reentry.dir/shuttle_reentry.cxx.i"
	cd /home/matt/git/psopt/debug/examples/shuttle_reentry && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/matt/git/psopt/examples/shuttle_reentry/shuttle_reentry.cxx > CMakeFiles/shuttle_reentry.dir/shuttle_reentry.cxx.i

examples/shuttle_reentry/CMakeFiles/shuttle_reentry.dir/shuttle_reentry.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/shuttle_reentry.dir/shuttle_reentry.cxx.s"
	cd /home/matt/git/psopt/debug/examples/shuttle_reentry && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/matt/git/psopt/examples/shuttle_reentry/shuttle_reentry.cxx -o CMakeFiles/shuttle_reentry.dir/shuttle_reentry.cxx.s

# Object files for target shuttle_reentry
shuttle_reentry_OBJECTS = \
"CMakeFiles/shuttle_reentry.dir/shuttle_reentry.cxx.o"

# External object files for target shuttle_reentry
shuttle_reentry_EXTERNAL_OBJECTS =

examples/shuttle_reentry/shuttle_reentry: examples/shuttle_reentry/CMakeFiles/shuttle_reentry.dir/shuttle_reentry.cxx.o
examples/shuttle_reentry/shuttle_reentry: examples/shuttle_reentry/CMakeFiles/shuttle_reentry.dir/build.make
examples/shuttle_reentry/shuttle_reentry: libPSOPT.a
examples/shuttle_reentry/shuttle_reentry: /usr/lib/libipopt.so
examples/shuttle_reentry/shuttle_reentry: /usr/lib/x86_64-linux-gnu/libdmumps_seq.so
examples/shuttle_reentry/shuttle_reentry: /usr/lib/x86_64-linux-gnu/libblas.so
examples/shuttle_reentry/shuttle_reentry: /usr/lib/x86_64-linux-gnu/liblapack.so
examples/shuttle_reentry/shuttle_reentry: /usr/lib/x86_64-linux-gnu/libdl.so
examples/shuttle_reentry/shuttle_reentry: /usr/lib/gcc/x86_64-linux-gnu/9/libgfortran.so
examples/shuttle_reentry/shuttle_reentry: /usr/lib/gcc/x86_64-linux-gnu/9/libquadmath.so
examples/shuttle_reentry/shuttle_reentry: /usr/lib/x86_64-linux-gnu/libdmumps_seq.so
examples/shuttle_reentry/shuttle_reentry: /usr/lib/x86_64-linux-gnu/libblas.so
examples/shuttle_reentry/shuttle_reentry: /usr/lib/x86_64-linux-gnu/liblapack.so
examples/shuttle_reentry/shuttle_reentry: /usr/lib/x86_64-linux-gnu/libdl.so
examples/shuttle_reentry/shuttle_reentry: /usr/lib/gcc/x86_64-linux-gnu/9/libgfortran.so
examples/shuttle_reentry/shuttle_reentry: /usr/lib/gcc/x86_64-linux-gnu/9/libquadmath.so
examples/shuttle_reentry/shuttle_reentry: /usr/lib/x86_64-linux-gnu/libadolc.so
examples/shuttle_reentry/shuttle_reentry: /usr/lib/x86_64-linux-gnu/libboost_system.so
examples/shuttle_reentry/shuttle_reentry: /usr/lib/x86_64-linux-gnu/libColPack.so
examples/shuttle_reentry/shuttle_reentry: /usr/lib/x86_64-linux-gnu/libm.so
examples/shuttle_reentry/shuttle_reentry: examples/shuttle_reentry/CMakeFiles/shuttle_reentry.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/matt/git/psopt/debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable shuttle_reentry"
	cd /home/matt/git/psopt/debug/examples/shuttle_reentry && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/shuttle_reentry.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/shuttle_reentry/CMakeFiles/shuttle_reentry.dir/build: examples/shuttle_reentry/shuttle_reentry

.PHONY : examples/shuttle_reentry/CMakeFiles/shuttle_reentry.dir/build

examples/shuttle_reentry/CMakeFiles/shuttle_reentry.dir/clean:
	cd /home/matt/git/psopt/debug/examples/shuttle_reentry && $(CMAKE_COMMAND) -P CMakeFiles/shuttle_reentry.dir/cmake_clean.cmake
.PHONY : examples/shuttle_reentry/CMakeFiles/shuttle_reentry.dir/clean

examples/shuttle_reentry/CMakeFiles/shuttle_reentry.dir/depend:
	cd /home/matt/git/psopt/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/matt/git/psopt /home/matt/git/psopt/examples/shuttle_reentry /home/matt/git/psopt/debug /home/matt/git/psopt/debug/examples/shuttle_reentry /home/matt/git/psopt/debug/examples/shuttle_reentry/CMakeFiles/shuttle_reentry.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/shuttle_reentry/CMakeFiles/shuttle_reentry.dir/depend

