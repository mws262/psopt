#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "PSOPT" for configuration "Debug"
set_property(TARGET PSOPT APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(PSOPT PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/PSOPT/libPSOPT.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS PSOPT )
list(APPEND _IMPORT_CHECK_FILES_FOR_PSOPT "${_IMPORT_PREFIX}/lib/PSOPT/libPSOPT.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
