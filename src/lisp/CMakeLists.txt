# src/lisp/CMakeLists.txt
#
# StreamVorti - Lisp integration module
# Copyright (C) 2020-2025 Benjamin F. Zwick

# Source files for the Lisp integration module
set(LISP_SOURCES
    ecl_runtime.cpp
    lisp_bridge.cpp
    lisp_mesh.cpp
    lisp_dcpse.cpp
    lisp_loader.cpp
)

# Create object library for Lisp integration
add_library(LispIntegration OBJECT ${LISP_SOURCES})

# Add include directories
target_include_directories(LispIntegration
    PUBLIC
        ${CMAKE_SOURCE_DIR}/include
        ${ECL_INCLUDE_DIRS}
        ${MFEM_INCLUDE_DIRS}
        ${EIGEN3_INCLUDE_DIRS}
)

# Set compile definitions
if(ECL_FOUND)
    target_compile_definitions(LispIntegration PUBLIC STREAMVORTI_WITH_ECL)
endif()

# Position independent code for shared library
set_property(TARGET LispIntegration PROPERTY POSITION_INDEPENDENT_CODE ON)
