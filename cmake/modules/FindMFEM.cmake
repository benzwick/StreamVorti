# FindMFEM.cmake
#
# Finds the MFEM library, supporting both CMake-based and simple library installations.
#
# This module will first try to use MFEM's CMake config files (MFEMConfig.cmake).
# If those are not found, it falls back to direct library and header detection
# (which handles Spack installations and other simple builds).
#
# Input variables:
#   MFEM_DIR - Path to MFEM installation prefix
#
# Output variables:
#   MFEM_FOUND - True if MFEM was found
#   MFEM_INCLUDE_DIRS - Include directories for MFEM
#   MFEM_LIBRARIES - Libraries to link against
#   MFEM_VERSION - Version of MFEM (extracted from config.hpp if available)
#
# Imported targets:
#   mfem - The MFEM library

# First, try to find MFEM using CMake config files
# From https://github.com/GLVis/glvis/blob/9f2df220342d13ef42c2985257ed96ba43f3bb70/CMakeLists.txt
# The following variables can be used to help CMake find MFEM:
#  * MFEM_DIR - absolute path to the MFEM build or install prefix.
#  * mfem_DIR - absolute path to where MFEMConfig.cmake is.

message(STATUS "FindMFEM: Starting MFEM detection")
message(STATUS "FindMFEM: MFEM_DIR = ${MFEM_DIR}")

if(MFEM_DIR)
    message(STATUS "FindMFEM: Searching for CMake config in:")
    message(STATUS "  - ${MFEM_DIR}")
    message(STATUS "  - ${MFEM_DIR}/lib/cmake/mfem")
    find_package(mfem QUIET NAMES MFEM HINTS "${MFEM_DIR}"
                 "${MFEM_DIR}/lib/cmake/mfem" NO_DEFAULT_PATH)
else()
    message(STATUS "FindMFEM: MFEM_DIR not set, searching default paths")
    find_package(mfem QUIET NAMES MFEM)
endif()

if(mfem_FOUND)
    message(STATUS "FindMFEM: CMake config found - mfem_FOUND = TRUE")
else()
    message(STATUS "FindMFEM: CMake config NOT found - falling back to direct detection")
endif()

if(mfem_FOUND)
    # CMake-based MFEM found
    message(STATUS "Found MFEM via CMake config: ${mfem_DIR}")

    # The config file should set these, but ensure they're available
    if(NOT MFEM_VERSION AND DEFINED mfem_VERSION)
        set(MFEM_VERSION ${mfem_VERSION})
    endif()

else()
    # CMake config not found - try direct library and header detection
    # This handles Spack installations and other simple builds
    message(STATUS "FindMFEM: CMake config not found - trying direct library/header detection")

    find_library(MFEM_LIBRARY
        NAMES mfem libmfem
        HINTS ${MFEM_DIR}
        PATH_SUFFIXES lib lib64
        DOC "MFEM library"
    )

    find_path(MFEM_INCLUDE_DIR
        NAMES mfem.hpp
        HINTS ${MFEM_DIR}
        PATH_SUFFIXES include
        DOC "MFEM include directory"
    )

    if(MFEM_LIBRARY AND MFEM_INCLUDE_DIR)
        message(STATUS "FindMFEM: SUCCESS - Found MFEM library and headers")
        message(STATUS "FindMFEM: Library: ${MFEM_LIBRARY}")
        message(STATUS "FindMFEM: Include: ${MFEM_INCLUDE_DIR}")

        set(MFEM_INCLUDE_DIRS ${MFEM_INCLUDE_DIR})
        set(MFEM_LIBRARIES ${MFEM_LIBRARY})

        # Try to extract version from config header
        if(EXISTS "${MFEM_INCLUDE_DIR}/mfem/config/config.hpp")
            file(READ "${MFEM_INCLUDE_DIR}/mfem/config/config.hpp" MFEM_CONFIG_CONTENT)
            string(REGEX MATCH "#define MFEM_VERSION_STRING \"([0-9.]+)\"" _ "${MFEM_CONFIG_CONTENT}")
            if(CMAKE_MATCH_1)
                set(MFEM_VERSION ${CMAKE_MATCH_1})
                message(STATUS "FindMFEM: Extracted version: ${MFEM_VERSION}")
            endif()
        endif()

        # Create imported target
        if(NOT TARGET mfem)
            add_library(mfem UNKNOWN IMPORTED)
            set_target_properties(mfem PROPERTIES
                IMPORTED_LOCATION "${MFEM_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${MFEM_INCLUDE_DIRS}"
            )
        endif()

        message(STATUS "FindMFEM: Created mfem imported target")
    else()
        message(STATUS "FindMFEM: FAILED - Could not find MFEM library or headers")
        if(NOT MFEM_LIBRARY)
            message(STATUS "  MFEM_LIBRARY not found in ${MFEM_DIR}/lib")
        endif()
        if(NOT MFEM_INCLUDE_DIR)
            message(STATUS "  MFEM_INCLUDE_DIR not found in ${MFEM_DIR}/include")
        endif()
    endif()
endif()

# Handle standard find_package arguments
include(FindPackageHandleStandardArgs)

# Check which method succeeded and validate accordingly
if(mfem_FOUND)
    # CMake config was found - check for MFEM_INCLUDE_DIRS set by MFEMConfig.cmake
    find_package_handle_standard_args(MFEM
        REQUIRED_VARS MFEM_INCLUDE_DIRS
        VERSION_VAR MFEM_VERSION
    )
else()
    # Direct detection - check for library and headers
    find_package_handle_standard_args(MFEM
        REQUIRED_VARS MFEM_LIBRARY MFEM_INCLUDE_DIR
        VERSION_VAR MFEM_VERSION
    )
endif()

if(MFEM_FOUND)
    message(STATUS "MFEM version: ${MFEM_VERSION}")
    message(STATUS "MFEM include dirs: ${MFEM_INCLUDE_DIRS}")
    message(STATUS "MFEM libraries: ${MFEM_LIBRARIES}")
endif()

mark_as_advanced(MFEM_LIBRARY MFEM_INCLUDE_DIR)
