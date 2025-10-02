# FindMFEM.cmake
#
# Finds the MFEM library, supporting both CMake-based and Makefile-based installations.
#
# This module will first try to use MFEM's CMake config files (MFEMConfig.cmake).
# If those are not found, it falls back to using mfem-config script (from Makefile builds).
#
# Input variables:
#   MFEM_DIR - Path to MFEM installation prefix
#
# Output variables:
#   MFEM_FOUND - True if MFEM was found
#   MFEM_INCLUDE_DIRS - Include directories for MFEM
#   MFEM_LIBRARIES - Libraries to link against
#   MFEM_VERSION - Version of MFEM
#   MFEM_CXX_COMPILER - C++ compiler used to build MFEM
#
# Imported targets:
#   mfem - The MFEM library

# First, try to find MFEM using CMake config files
find_package(mfem CONFIG QUIET)

if(mfem_FOUND)
    # CMake-based MFEM found
    message(STATUS "Found MFEM via CMake config: ${mfem_DIR}")

    # The config file should set these, but ensure they're available
    if(NOT MFEM_VERSION AND DEFINED mfem_VERSION)
        set(MFEM_VERSION ${mfem_VERSION})
    endif()

else()
    # Try to find mfem-config script (Makefile-based installation)
    find_program(MFEM_CONFIG
        NAMES mfem-config
        HINTS ${MFEM_DIR}
        PATH_SUFFIXES bin
        DOC "MFEM configuration script"
    )

    if(MFEM_CONFIG)
        message(STATUS "Found mfem-config: ${MFEM_CONFIG}")

        # Execute mfem-config to get compilation and link flags
        execute_process(
            COMMAND ${MFEM_CONFIG} --cxx
            OUTPUT_VARIABLE MFEM_CXX_COMPILER
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )

        execute_process(
            COMMAND ${MFEM_CONFIG} --cppflags
            OUTPUT_VARIABLE MFEM_CXX_FLAGS
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )

        execute_process(
            COMMAND ${MFEM_CONFIG} --incflags
            OUTPUT_VARIABLE MFEM_INCLUDE_FLAGS
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )

        execute_process(
            COMMAND ${MFEM_CONFIG} --libs
            OUTPUT_VARIABLE MFEM_LINK_FLAGS
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )

        execute_process(
            COMMAND ${MFEM_CONFIG} --version
            OUTPUT_VARIABLE MFEM_VERSION
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )

        execute_process(
            COMMAND ${MFEM_CONFIG} --prefix
            OUTPUT_VARIABLE MFEM_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )

        # Parse include directories from flags
        string(REGEX MATCHALL "-I[^ ]+" MFEM_INCLUDE_FLAGS_LIST "${MFEM_INCLUDE_FLAGS}")
        set(MFEM_INCLUDE_DIRS "")
        foreach(flag ${MFEM_INCLUDE_FLAGS_LIST})
            string(SUBSTRING ${flag} 2 -1 include_dir)
            list(APPEND MFEM_INCLUDE_DIRS ${include_dir})
        endforeach()

        # Parse libraries from link flags
        string(REGEX MATCHALL "-L[^ ]+" MFEM_LIBRARY_DIRS_LIST "${MFEM_LINK_FLAGS}")
        string(REGEX MATCHALL "-l[^ ]+" MFEM_LIBRARIES_LIST "${MFEM_LINK_FLAGS}")

        set(MFEM_LIBRARY_DIRS "")
        foreach(flag ${MFEM_LIBRARY_DIRS_LIST})
            string(SUBSTRING ${flag} 2 -1 lib_dir)
            list(APPEND MFEM_LIBRARY_DIRS ${lib_dir})
        endforeach()

        set(MFEM_LIBRARIES "")
        foreach(flag ${MFEM_LIBRARIES_LIST})
            string(SUBSTRING ${flag} 2 -1 lib_name)
            list(APPEND MFEM_LIBRARIES ${lib_name})
        endforeach()

        # Add full path to MFEM library
        if(MFEM_PREFIX)
            list(INSERT MFEM_LIBRARIES 0 "${MFEM_PREFIX}/lib/libmfem.a")
        endif()

        # Create imported target
        if(NOT TARGET mfem)
            add_library(mfem INTERFACE IMPORTED)
            set_target_properties(mfem PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES "${MFEM_INCLUDE_DIRS}"
                INTERFACE_LINK_LIBRARIES "${MFEM_LIBRARIES}"
            )

            # Add any additional compile flags
            if(MFEM_CXX_FLAGS)
                set_target_properties(mfem PROPERTIES
                    INTERFACE_COMPILE_OPTIONS "${MFEM_CXX_FLAGS}"
                )
            endif()
        endif()
    endif()
endif()

# Handle standard find_package arguments
include(FindPackageHandleStandardArgs)

# For CMake-based MFEM, check MFEM_LIBRARY_DIR (set by MFEMConfig.cmake)
# For Makefile-based MFEM, check MFEM_CONFIG (path to mfem-config script)
if(mfem_FOUND)
    find_package_handle_standard_args(MFEM
        REQUIRED_VARS MFEM_LIBRARY_DIR
        VERSION_VAR MFEM_VERSION
    )
else()
    find_package_handle_standard_args(MFEM
        REQUIRED_VARS MFEM_CONFIG
        VERSION_VAR MFEM_VERSION
    )
endif()

if(MFEM_FOUND)
    message(STATUS "MFEM version: ${MFEM_VERSION}")
    message(STATUS "MFEM include dirs: ${MFEM_INCLUDE_DIRS}")
    message(STATUS "MFEM libraries: ${MFEM_LIBRARIES}")
    if(MFEM_CXX_COMPILER)
        message(STATUS "MFEM C++ compiler: ${MFEM_CXX_COMPILER}")
    endif()
endif()

mark_as_advanced(MFEM_CONFIG MFEM_CXX_COMPILER)
