# FindMFEM.cmake
#
# Finds the MFEM library, supporting both CMake-based and Make-based installations.
#
# This module will first try to use MFEM's CMake config files (MFEMConfig.cmake).
# If those are not found, it parses config.mk from Make-based MFEM installations
# (Spack and traditional Make builds).
#
# Input variables:
#   MFEM_DIR - Path to MFEM installation prefix
#
# Output variables:
#   MFEM_FOUND - True if MFEM was found
#   MFEM_INCLUDE_DIRS - Include directories for MFEM
#   MFEM_LIBRARIES - Libraries to link against
#   MFEM_VERSION - Version of MFEM
#   MFEM_CXX_COMPILER - C++ compiler used to build MFEM (from config.mk only)
#   MFEM_CXX_FLAGS - C++ flags used to build MFEM (from config.mk only)
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
    # CMake config not found - parse config.mk from Make-based MFEM installation
    message(STATUS "FindMFEM: CMake config not found - looking for config.mk")

    find_file(MFEM_CONFIG_MK
        NAMES config.mk
        PATHS ${MFEM_DIR}/share/mfem
        NO_DEFAULT_PATH
        DOC "MFEM config.mk file"
    )

    if(NOT MFEM_CONFIG_MK)
        message(FATAL_ERROR "FindMFEM: Could not find config.mk at ${MFEM_DIR}/share/mfem/config.mk")
    endif()

    message(STATUS "FindMFEM: Found config.mk: ${MFEM_CONFIG_MK}")

    # Read and parse config.mk
    file(STRINGS "${MFEM_CONFIG_MK}" MFEM_CONFIG_LINES)

    # Extract variables from config.mk
    foreach(line ${MFEM_CONFIG_LINES})
        if(line MATCHES "^MFEM_VERSION_STRING[ ]*=[ ]*(.+)$")
            set(MFEM_VERSION "${CMAKE_MATCH_1}")
            string(STRIP "${MFEM_VERSION}" MFEM_VERSION)
        elseif(line MATCHES "^MFEM_CXX[ ]*=[ ]*(.+)$")
            set(MFEM_CXX_COMPILER "${CMAKE_MATCH_1}")
            string(STRIP "${MFEM_CXX_COMPILER}" MFEM_CXX_COMPILER)
        elseif(line MATCHES "^MFEM_CXXFLAGS[ ]*=[ ]*(.+)$")
            set(MFEM_CXX_FLAGS "${CMAKE_MATCH_1}")
            string(STRIP "${MFEM_CXX_FLAGS}" MFEM_CXX_FLAGS)
        elseif(line MATCHES "^MFEM_INC_DIR[ ]*=[ ]*(.+)$")
            set(MFEM_INCLUDE_DIRS "${CMAKE_MATCH_1}")
            string(STRIP "${MFEM_INCLUDE_DIRS}" MFEM_INCLUDE_DIRS)
        elseif(line MATCHES "^MFEM_LIB_DIR[ ]*=[ ]*(.+)$")
            set(MFEM_LIB_DIR "${CMAKE_MATCH_1}")
            string(STRIP "${MFEM_LIB_DIR}" MFEM_LIB_DIR)
        elseif(line MATCHES "^MFEM_TPLFLAGS[ ]*=[ ]*(.+)$")
            set(MFEM_TPLFLAGS "${CMAKE_MATCH_1}")
            string(STRIP "${MFEM_TPLFLAGS}" MFEM_TPLFLAGS)
        elseif(line MATCHES "^MFEM_EXT_LIBS[ ]*=[ ]*(.+)$")
            set(MFEM_EXT_LIBS "${CMAKE_MATCH_1}")
            string(STRIP "${MFEM_EXT_LIBS}" MFEM_EXT_LIBS)
        endif()
    endforeach()

    # Parse MFEM_TPLFLAGS to extract third-party library include directories
    # MFEM_TPLFLAGS contains flags like: -I/path/to/hypre/include -I/path/to/metis/include
    if(MFEM_TPLFLAGS)
        string(REGEX MATCHALL "-I([^ ]+)" TPL_INCLUDES "${MFEM_TPLFLAGS}")
        foreach(inc ${TPL_INCLUDES})
            string(REGEX REPLACE "^-I" "" inc_path "${inc}")
            list(APPEND MFEM_INCLUDE_DIRS "${inc_path}")
        endforeach()
    endif()

    # Construct library path from MFEM_LIB_DIR
    set(MFEM_LIBRARY "${MFEM_LIB_DIR}/libmfem.a")

    # Include external libraries that MFEM depends on
    if(MFEM_EXT_LIBS)
        # Convert the space-separated string to a list and append to MFEM_LIBRARIES
        separate_arguments(MFEM_EXT_LIBS_LIST UNIX_COMMAND "${MFEM_EXT_LIBS}")
        set(MFEM_LIBRARIES ${MFEM_LIBRARY} ${MFEM_EXT_LIBS_LIST})
    else()
        set(MFEM_LIBRARIES ${MFEM_LIBRARY})
    endif()

    message(STATUS "FindMFEM: Parsed config.mk successfully")
    message(STATUS "  Version: ${MFEM_VERSION}")
    message(STATUS "  Compiler: ${MFEM_CXX_COMPILER}")
    message(STATUS "  Include: ${MFEM_INCLUDE_DIRS}")
    message(STATUS "  Library: ${MFEM_LIBRARY}")

    # Create imported target
    if(NOT TARGET mfem)
        add_library(mfem UNKNOWN IMPORTED)
        set_target_properties(mfem PROPERTIES
            IMPORTED_LOCATION "${MFEM_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${MFEM_INCLUDE_DIRS}"
        )
    endif()
endif()

# Handle standard find_package arguments
include(FindPackageHandleStandardArgs)

# Check which method succeeded and validate accordingly
if(mfem_FOUND)
    # CMake config was found
    find_package_handle_standard_args(MFEM
        REQUIRED_VARS MFEM_INCLUDE_DIRS
        VERSION_VAR MFEM_VERSION
    )
else()
    # config.mk was parsed
    find_package_handle_standard_args(MFEM
        REQUIRED_VARS MFEM_LIBRARY MFEM_INCLUDE_DIRS MFEM_CXX_COMPILER
        VERSION_VAR MFEM_VERSION
    )
endif()

if(MFEM_FOUND)
    message(STATUS "MFEM version: ${MFEM_VERSION}")
    message(STATUS "MFEM include dirs: ${MFEM_INCLUDE_DIRS}")
    message(STATUS "MFEM libraries: ${MFEM_LIBRARIES}")
endif()

mark_as_advanced(MFEM_CONFIG_MK MFEM_LIBRARY)
