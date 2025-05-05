#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "umf::umf" for configuration "Release"
set_property(TARGET umf::umf APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(umf::umf PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libumf.so.0.10.0"
  IMPORTED_SONAME_RELEASE "libumf.so.0"
  )

list(APPEND _cmake_import_check_targets umf::umf )
list(APPEND _cmake_import_check_files_for_umf::umf "${_IMPORT_PREFIX}/lib/libumf.so.0.10.0" )

# Import target "umf::disjoint_pool" for configuration "Release"
set_property(TARGET umf::disjoint_pool APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(umf::disjoint_pool PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C;CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libdisjoint_pool.a"
  )

list(APPEND _cmake_import_check_targets umf::disjoint_pool )
list(APPEND _cmake_import_check_files_for_umf::disjoint_pool "${_IMPORT_PREFIX}/lib/libdisjoint_pool.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
