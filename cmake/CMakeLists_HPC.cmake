# This CMake is for StreamVorti simulation running in a Spack environment on HPC system

project(StreamVorti VERSION 2.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build type configuration
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Compiler flags for optimisation (GCC 14.3.0 specific)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -mtune=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wpedantic")


# =============================================================================
# Profiling configuration (gprof)
# =============================================================================
option(ENABLE_PROFILING "Enable gprof profiling" OFF)

if(ENABLE_PROFILING)
    message(STATUS "Profiling enabled with gprof (-pg flag)")
    # Add -pg flag for both compilation and linking
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")

    # Use -O2 instead of -O3 for profiling to get more accurate results
    # -O3 can inline functions which makes profiling less meaningful
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=native -mtune=native -DNDEBUG -pg")

    # Note: -pg adds overhead (~10-30%), but provides detailed profiling data
endif()

# Export compile commands for IDE support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# RPATH Configuration
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Add GCC library path for OpenMP
set(CMAKE_INSTALL_RPATH
    "/uwahpc/rocky9/devel/compiler/gcc/14.3.0/lib64"
    "${SPACK_VIEW_PATH}/lib"
    "${SPACK_VIEW_PATH}/lib64"
)

# =============================================================================
# Detect Spack environment
# =============================================================================

if(DEFINED ENV{SPACK_ENV})
  message(STATUS "Detected Spack environment: $ENV{SPACK_ENV}")
  set(SPACK_ENV_PATH $ENV{SPACK_ENV})

  # Set Spack view path
  set(SPACK_VIEW_PATH "${SPACK_ENV_PATH}/.spack-env/view")
  if(NOT EXISTS ${SPACK_VIEW_PATH})
    set(SPACK_VIEW_PATH "${SPACK_ENV_PATH}/view")
  endif()

  message(STATUS "Using Spack view: ${SPACK_VIEW_PATH}")
  list(APPEND CMAKE_PREFIX_PATH ${SPACK_VIEW_PATH})
endif()

# =============================================================================
# MPI Configuration
# =============================================================================

option(ENABLE_MPI "Enable MPI parallel support" ON)

if(ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS CXX)

  if(MPI_FOUND)
    message(STATUS "MPI Configuration:")
    message(STATUS "  MPI C++ compiler:     ${MPI_CXX_COMPILER}")
    message(STATUS "  MPI include dirs:     ${MPI_CXX_INCLUDE_DIRS}")
    message(STATUS "  MPI libraries:        ${MPI_CXX_LIBRARIES}")

    # Use MPI compiler wrapper
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

    add_compile_definitions(MFEM_USE_MPI)
    add_compile_definitions(STREAMVORTI_USE_MPI)
  endif()
endif()

# =============================================================================
# Find MFEM (Make-based installation from Spack)
# =============================================================================

message(STATUS "Searching for MFEM...")

find_path(MFEM_INCLUDE_DIR
  NAMES mfem.hpp
  HINTS ${SPACK_VIEW_PATH}
  PATH_SUFFIXES include
)

find_library(MFEM_LIBRARY
  NAMES mfem
  HINTS ${SPACK_VIEW_PATH}
  PATH_SUFFIXES lib lib64
)

if(MFEM_INCLUDE_DIR AND MFEM_LIBRARY)
  set(MFEM_FOUND TRUE)
  set(MFEM_INCLUDE_DIRS ${MFEM_INCLUDE_DIR})
  set(MFEM_LIBRARIES ${MFEM_LIBRARY})

  message(STATUS "MFEM found:")
  message(STATUS "  Include dir: ${MFEM_INCLUDE_DIR}")
  message(STATUS "  Library:     ${MFEM_LIBRARY}")

  # Check for MFEM features
  find_file(MFEM_CONFIG_HPP
    NAMES config.hpp
    HINTS ${MFEM_INCLUDE_DIR}
    PATH_SUFFIXES mfem config
    NO_DEFAULT_PATH
  )

  if(MFEM_CONFIG_HPP)
    file(READ ${MFEM_CONFIG_HPP} MFEM_CONFIG_CONTENT)

    if(MFEM_CONFIG_CONTENT MATCHES "MFEM_USE_MPI")
      message(STATUS "  MPI support:     Enabled")
    endif()
    if(MFEM_CONFIG_CONTENT MATCHES "MFEM_USE_METIS")
      message(STATUS "  METIS support:   Enabled")
    endif()
    if(MFEM_CONFIG_CONTENT MATCHES "MFEM_USE_HYPRE")
      message(STATUS "  HYPRE support:   Enabled")
    endif()
    if(MFEM_CONFIG_CONTENT MATCHES "MFEM_USE_SUITESPARSE")
      message(STATUS "  SuiteSparse:     Enabled")
      add_compile_definitions(MFEM_USE_SUITESPARSE)
    endif()
  endif()
else()
  message(FATAL_ERROR "MFEM not found. Ensure Spack environment is activated.")
endif()

# =============================================================================
# Find Eigen3
# =============================================================================

find_path(EIGEN3_INCLUDE_DIR
  NAMES Eigen/Core
  HINTS ${SPACK_VIEW_PATH}
  PATH_SUFFIXES include include/eigen3
)

if(EIGEN3_INCLUDE_DIR)
  message(STATUS "Eigen3 found:")
  message(STATUS "  Include dir: ${EIGEN3_INCLUDE_DIR}")
else()
  message(FATAL_ERROR "Eigen3 not found.")
endif()

# =============================================================================
# Find CGAL - Enhanced detection for Spack installation
# =============================================================================

message(STATUS "Searching for CGAL...")

find_path(CGAL_INCLUDE_DIR
  NAMES CGAL/basic.h CGAL/version.h
  HINTS ${SPACK_VIEW_PATH}
  PATH_SUFFIXES include
)

# Try different possible library names
set(CGAL_LIBRARY "CGAL_LIBRARY-NOTFOUND")
foreach(lib_name CGAL CGAL_Core cgal)
  find_library(CGAL_LIB_TEMP
    NAMES ${lib_name}
    HINTS ${SPACK_VIEW_PATH}
    PATH_SUFFIXES lib lib64
  )
  if(CGAL_LIB_TEMP)
    set(CGAL_LIBRARY ${CGAL_LIB_TEMP})
    break()
  endif()
  unset(CGAL_LIB_TEMP CACHE)
endforeach()

if(CGAL_INCLUDE_DIR)
  message(STATUS "CGAL found:")
  message(STATUS "  Include dir: ${CGAL_INCLUDE_DIR}")

  if(CGAL_LIBRARY)
    message(STATUS "  Library:     ${CGAL_LIBRARY}")
    set(CGAL_LIBRARIES ${CGAL_LIBRARY})
  else()
    message(STATUS "  Library:     Header-only")
    set(CGAL_LIBRARIES "")
  endif()

  set(CGAL_FOUND TRUE)

  # CGAL dependencies
  find_library(GMP_LIBRARY NAMES gmp HINTS ${SPACK_VIEW_PATH} PATH_SUFFIXES lib lib64)
  find_library(MPFR_LIBRARY NAMES mpfr HINTS ${SPACK_VIEW_PATH} PATH_SUFFIXES lib lib64)

  if(GMP_LIBRARY)
    list(APPEND CGAL_LIBRARIES ${GMP_LIBRARY})
    message(STATUS "  GMP library: ${GMP_LIBRARY}")
  endif()

  if(MPFR_LIBRARY)
    list(APPEND CGAL_LIBRARIES ${MPFR_LIBRARY})
    message(STATUS "  MPFR library: ${MPFR_LIBRARY}")
  endif()

  # Boost libraries for CGAL
  foreach(boost_lib thread system)
    find_library(BOOST_${boost_lib}_LIB
      NAMES boost_${boost_lib}
      HINTS ${SPACK_VIEW_PATH}
      PATH_SUFFIXES lib lib64
    )
    if(BOOST_${boost_lib}_LIB)
      list(APPEND CGAL_LIBRARIES ${BOOST_${boost_lib}_LIB})
      message(STATUS "  Boost.${boost_lib}: ${BOOST_${boost_lib}_LIB}")
    endif()
  endforeach()
else()
  message(FATAL_ERROR "CGAL not found in ${SPACK_VIEW_PATH}")
endif()

# =============================================================================
# Find additional dependencies for MFEM
# =============================================================================

# zlib (required by MFEM for mesh compression)
find_library(ZLIB_LIBRARY
  NAMES z zlib
  HINTS ${SPACK_VIEW_PATH}
  PATH_SUFFIXES lib lib64
)

if(ZLIB_LIBRARY)
  message(STATUS "zlib found: ${ZLIB_LIBRARY}")
  list(APPEND MFEM_LIBRARIES ${ZLIB_LIBRARY})
else()
  message(WARNING "zlib not found - MFEM may have linking issues")
endif()

# METIS
find_library(METIS_LIBRARY NAMES metis HINTS ${SPACK_VIEW_PATH} PATH_SUFFIXES lib lib64)
if(METIS_LIBRARY)
  message(STATUS "METIS found: ${METIS_LIBRARY}")
  list(APPEND MFEM_LIBRARIES ${METIS_LIBRARY})
endif()

# ParMETIS (for MPI)
if(ENABLE_MPI)
  find_library(PARMETIS_LIBRARY NAMES parmetis HINTS ${SPACK_VIEW_PATH} PATH_SUFFIXES lib lib64)
  if(PARMETIS_LIBRARY)
    message(STATUS "ParMETIS found: ${PARMETIS_LIBRARY}")
    list(APPEND MFEM_LIBRARIES ${PARMETIS_LIBRARY})
  endif()
endif()

# HYPRE (for MPI)
if(ENABLE_MPI)
  find_library(HYPRE_LIBRARY NAMES HYPRE HINTS ${SPACK_VIEW_PATH} PATH_SUFFIXES lib lib64)
  if(HYPRE_LIBRARY)
    message(STATUS "HYPRE found: ${HYPRE_LIBRARY}")
    list(APPEND MFEM_LIBRARIES ${HYPRE_LIBRARY})
    add_compile_definitions(MFEM_USE_HYPRE)
  endif()
endif()

# SuiteSparse/UMFPACK
find_library(UMFPACK_LIBRARY NAMES umfpack HINTS ${SPACK_VIEW_PATH} PATH_SUFFIXES lib lib64)
if(UMFPACK_LIBRARY)
  message(STATUS "UMFPACK found: ${UMFPACK_LIBRARY}")
  list(APPEND MFEM_LIBRARIES ${UMFPACK_LIBRARY})

  # SuiteSparse components
  foreach(ss_lib amd cholmod suitesparseconfig)
    find_library(${ss_lib}_LIBRARY
      NAMES ${ss_lib}
      HINTS ${SPACK_VIEW_PATH}
      PATH_SUFFIXES lib lib64
    )
    if(${ss_lib}_LIBRARY)
      list(APPEND MFEM_LIBRARIES ${${ss_lib}_LIBRARY})
      message(STATUS "  ${ss_lib}: ${${ss_lib}_LIBRARY}")
    endif()
  endforeach()
endif()

# Comprehensive SuiteSparse detection
message(STATUS "Searching for SuiteSparse components...")

set(SUITESPARSE_COMPONENTS
  umfpack
  klu
  btf
  colamd
  amd
  cholmod
  ccolamd
  camd
  suitesparseconfig
)

foreach(component ${SUITESPARSE_COMPONENTS})
  string(TOUPPER ${component} component_upper)
  find_library(${component_upper}_LIBRARY
    NAMES ${component}
    HINTS ${SPACK_VIEW_PATH}
    PATH_SUFFIXES lib lib64
  )

  if(${component_upper}_LIBRARY)
    message(STATUS "  ${component}: ${${component_upper}_LIBRARY}")
    list(APPEND MFEM_LIBRARIES ${${component_upper}_LIBRARY})
  else()
    message(STATUS "  ${component}: NOT FOUND")
  endif()
endforeach()

# BLAS/LAPACK
find_library(BLAS_LIBRARY NAMES blas openblas HINTS ${SPACK_VIEW_PATH} PATH_SUFFIXES lib lib64)
if(BLAS_LIBRARY)
  message(STATUS "BLAS found: ${BLAS_LIBRARY}")
  list(APPEND MFEM_LIBRARIES ${BLAS_LIBRARY})
endif()

find_library(LAPACK_LIBRARY NAMES lapack HINTS ${SPACK_VIEW_PATH} PATH_SUFFIXES lib lib64)
if(LAPACK_LIBRARY)
  message(STATUS "LAPACK found: ${LAPACK_LIBRARY}")
  list(APPEND MFEM_LIBRARIES ${LAPACK_LIBRARY})
endif()

# Standard math and threading libraries
list(APPEND MFEM_LIBRARIES m pthread)

# gfortran (often needed for LAPACK/BLAS)
find_library(GFORTRAN_LIBRARY NAMES gfortran HINTS ${SPACK_VIEW_PATH} PATH_SUFFIXES lib lib64)
if(GFORTRAN_LIBRARY)
  message(STATUS "gfortran found: ${GFORTRAN_LIBRARY}")
  list(APPEND MFEM_LIBRARIES ${GFORTRAN_LIBRARY})
endif()

# =============================================================================
# Collect source files
# =============================================================================

file(GLOB_RECURSE HEADER_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h
)

file(GLOB APPROXIMANTS_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/approximants/*.cpp
)

file(GLOB SUPPORT_DOMAIN_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/support_domain/*.cpp
)

set(LIBRARY_SOURCES
    ${APPROXIMANTS_SOURCES}
    ${SUPPORT_DOMAIN_SOURCES}
)

message(STATUS "Library source files:")
foreach(src ${LIBRARY_SOURCES})
  get_filename_component(src_name ${src} NAME)
  message(STATUS "  ${src_name}")
endforeach()

# =============================================================================
# Build StreamVorti library
# =============================================================================

add_library(streamvorti ${LIBRARY_SOURCES} ${HEADER_FILES})

target_include_directories(streamvorti
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
        ${MFEM_INCLUDE_DIRS}
        ${EIGEN3_INCLUDE_DIR}
        ${CGAL_INCLUDE_DIR}
)

target_link_libraries(streamvorti
    PUBLIC
        ${MFEM_LIBRARIES}
        ${CGAL_LIBRARIES}
)

# Link MPI
if(ENABLE_MPI AND MPI_FOUND)
    target_link_libraries(streamvorti PUBLIC MPI::MPI_CXX)
endif()

# Compiler warnings
target_compile_options(streamvorti PRIVATE
    -Wall -Wextra -Wpedantic
    -Wno-unused-parameter
    -Wno-sign-compare
)

# =============================================================================
# Build main executable
# =============================================================================
add_executable(StreamVorti ${CMAKE_CURRENT_SOURCE_DIR}/streamvorti.cpp)

target_include_directories(StreamVorti
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(StreamVorti
    PRIVATE
        streamvorti
)

target_compile_options(StreamVorti PRIVATE
    -Wall -Wextra -Wpedantic
    -Wno-unused-parameter
)

# =============================================================================
# OpenMP support
# =============================================================================

option(ENABLE_OPENMP "Enable OpenMP for thread parallelism" ON)

if(ENABLE_OPENMP)
    find_package(OpenMP REQUIRED)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP Configuration:")
        message(STATUS "OpenMP version: ${OpenMP_CXX_VERSION}")
        message(STATUS "  OpenMP flags: ${OpenMP_CXX_FLAGS}")

        target_link_libraries(streamvorti PUBLIC OpenMP::OpenMP_CXX)
        target_link_libraries(StreamVorti PRIVATE OpenMP::OpenMP_CXX)

        add_compile_definitions(STREAMVORTI_USE_OPENMP)

        # Set default number of threads at compile time if desired
        target_compile_definitions(streamvorti PUBLIC OMP_NUM_THREADS=1) #serial
    else()
        message(WARNING "OpenMP requested but not found")
    endif()
endif()

# =============================================================================
# Kaya HPC optimisations
# =============================================================================

option(ENABLE_AMD_OPTIMISATION "Enable AMD-specific optimisations" OFF)

if(ENABLE_AMD_OPTIMISATION)
    message(STATUS "Enabling AMD (AMD EPYC Zen4) optimisations")
    target_compile_options(streamvorti PUBLIC -march=znver4 -mtune=znver4)
    target_compile_options(StreamVorti PRIVATE -march=znver4 -mtune=znver4)
endif()

# =============================================================================
# Installation
# =============================================================================

include(GNUInstallDirs)

install(TARGETS streamvorti StreamVorti
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h"
)

# =============================================================================
# Configuration Summary
# =============================================================================

message(STATUS "")
message(STATUS "========================================")
message(STATUS "StreamVorti Build Configuration")
message(STATUS "========================================")
message(STATUS "")
message(STATUS "General:")
message(STATUS "  Version:              ${PROJECT_VERSION}")
message(STATUS "  Build type:           ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ compiler:         ${CMAKE_CXX_COMPILER}")
message(STATUS "  Install prefix:       ${CMAKE_INSTALL_PREFIX}")
message(STATUS "  Profiling:            ${ENABLE_PROFILING}")
if(DEFINED SPACK_VIEW_PATH)
message(STATUS "  Spack view:           ${SPACK_VIEW_PATH}")
endif()
message(STATUS "")
message(STATUS "Dependencies:")
message(STATUS "  MFEM:                 Found")
message(STATUS "  Eigen3:               Found")
message(STATUS "  CGAL:                 Found")
message(STATUS "  zlib:                 ${ZLIB_LIBRARY}")
message(STATUS "")
message(STATUS "Parallel Computing:")
message(STATUS "  MPI:                  ${ENABLE_MPI}")
message(STATUS "  OpenMP:               ${ENABLE_OPENMP}")
message(STATUS "")
message(STATUS "========================================")
message(STATUS "")