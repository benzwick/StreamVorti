#--------------------------------------------------------------
CMAKE_MINIMUM_REQUIRED(VERSION 3.30)

IF(POLICY CMP0181)
  CMAKE_POLICY(SET CMP0181 NEW) # CMake 3.30
ENDIF()

# Extract MFEM compiler before PROJECT() to enable using the same compiler as MFEM
# This is needed when MFEM was built with MPI using an MPI compiler wrapper
if(DEFINED MFEM_DIR AND EXISTS "${MFEM_DIR}/share/mfem/config.mk")
  file(STRINGS "${MFEM_DIR}/share/mfem/config.mk" MFEM_CONFIG_LINES)
  foreach(line ${MFEM_CONFIG_LINES})
    if(line MATCHES "^MFEM_CXX[ ]*=[ ]*(.+)$")
      set(MFEM_CXX_COMPILER "${CMAKE_MATCH_1}")
      string(STRIP "${MFEM_CXX_COMPILER}" MFEM_CXX_COMPILER)
      if(NOT CMAKE_CXX_COMPILER)
        set(CMAKE_CXX_COMPILER "${MFEM_CXX_COMPILER}")
        message(STATUS "Using MFEM compiler: ${CMAKE_CXX_COMPILER}")
      endif()
      break()
    endif()
  endforeach()
endif()

PROJECT(StreamVorti VERSION 2.0.0)


# Set project's CMake modules path
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules")

# Library compilation options.
OPTION(STREAMVORTI_BUILD_SHARED_LIBS "Build as shared library" OFF )
OPTION(STREAMVORTI_BUILD_STATIC_LIBS "Build as static library" ON )
OPTION(ENABLE_PROFILING "Enable gprof profiling" OFF )
OPTION(STREAMVORTI_WITH_ECL "Enable ECL (Embeddable Common Lisp) for SDL support" OFF )
OPTION(STREAMVORTI_BUILD_TESTS "Build unit tests" OFF )

# Set the version of C++.
SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_EXTENSIONS OFF)

# Add profiling flags if enabled
IF(ENABLE_PROFILING AND "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg -g")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg -g")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
    MESSAGE(STATUS "Profiling enabled with gprof")
ENDIF()

# Set additional C++ flags for supported compilers.
IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Ob2")
ELSEIF (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mfpmath=sse -O3 -march=native -msse2 -ffast-math -fPIC -Wall")
    #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfpmath=sse -O3 -march=native -msse2 -ffast-math -fPIC -Wall")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfpmath=sse -O3 -msse2 -fPIC -Wall")
ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -frounding-math -fsignaling-nans -mfpmath=sse -march=native -msse2 -O3 -ffast-math -fPIC -Wall")
    #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -frounding-math -fsignaling-nans -mfpmath=sse -march=native -msse2 -O3 -ffast-math -fPIC -Wall")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -frounding-math -fsignaling-nans -mfpmath=sse -msse2 -O3 -fPIC -Wall")
ENDIF ()

# Include path to the headers of StreamVorti.
INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/include )


# Set the headers of the StreamVorti library.
SET(HEADERS  ${CMAKE_SOURCE_DIR}/include/StreamVorti/mfem_main.hpp )


# Find dependencies of StreamVorti.
# -------------------------

# Find CGAL
# See https://doc.cgal.org/latest/Manual/devman_create_and_use_a_cmakelist.html
FIND_PACKAGE(CGAL  REQUIRED COMPONENTS Core)
MESSAGE(STATUS "Found CGAL version ${CGAL_VERSION}")
IF (CGAL_VERSION VERSION_LESS "5.0")
    MESSAGE(FATAL_ERROR "CGAL version 5.0 or higher is required")
ENDIF()

# Find Eigen
FIND_PACKAGE(Eigen3  REQUIRED)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIRS})

# Find OpenMP
if(APPLE)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(OpenMP_CXX_FLAGS "-fopenmp=libomp")
        set(OpenMP_CXX_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY omp)
    endif()
endif()

FIND_PACKAGE(OpenMP COMPONENTS CXX REQUIRED)

# Import MFEM. The following variables can be used to help CMake find MFEM:
#  * MFEM_DIR - absolute path to the MFEM build or install prefix.
# FindMFEM.cmake will automatically detect whether MFEM was built with CMake or Make
message(STATUS "Looking for mfem ...")
set(MFEM_DIR "" CACHE PATH "Path to the MFEM build or install prefix.")
find_package(MFEM REQUIRED)

# Include paths and libraries needed by MFEM
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MFEM_CXX_FLAGS}")
INCLUDE_DIRECTORIES(${MFEM_INCLUDE_DIRS})

# Find ECL (Embeddable Common Lisp) for SDL support
IF(STREAMVORTI_WITH_ECL)
    # Try pkg-config first, fall back to manual detection
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
        pkg_check_modules(ECL QUIET ecl)
    endif()

    if(NOT ECL_FOUND)
        # Manual detection for systems without pkg-config file (e.g., Ubuntu)
        find_path(ECL_INCLUDE_DIRS ecl/ecl.h
            PATHS /usr/include /usr/local/include
        )
        find_library(ECL_LIBRARIES ecl
            PATHS /usr/lib /usr/lib/x86_64-linux-gnu /usr/local/lib
        )
        if(ECL_INCLUDE_DIRS AND ECL_LIBRARIES)
            set(ECL_FOUND TRUE)
            # Try to get version from ecl-config or header
            execute_process(
                COMMAND ecl-config --version
                OUTPUT_VARIABLE ECL_VERSION
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
            )
            if(NOT ECL_VERSION)
                set(ECL_VERSION "unknown")
            endif()
        endif()
    endif()

    if(NOT ECL_FOUND)
        message(FATAL_ERROR "ECL (Embeddable Common Lisp) not found. "
            "Install with: apt install ecl libecl-dev")
    endif()

    message(STATUS "Found ECL version ${ECL_VERSION}")
    message(STATUS "  ECL include dirs: ${ECL_INCLUDE_DIRS}")
    message(STATUS "  ECL libraries: ${ECL_LIBRARIES}")
    INCLUDE_DIRECTORIES(${ECL_INCLUDE_DIRS})
    add_compile_definitions(STREAMVORTI_WITH_ECL)
ENDIF()

# Print compiler information for all components
message(STATUS "========================================")
message(STATUS "Compiler Configuration:")
message(STATUS "  StreamVorti C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "  MFEM was built with:      ${MFEM_CXX_COMPILER}")
message(STATUS "  CGAL version:             ${CGAL_VERSION}")
message(STATUS "  Eigen3 include dir:       ${EIGEN3_INCLUDE_DIRS}")
IF(STREAMVORTI_WITH_ECL)
message(STATUS "  ECL version:              ${ECL_VERSION}")
message(STATUS "  SDL support:              ENABLED")
ELSE()
message(STATUS "  SDL support:              DISABLED (use -DSTREAMVORTI_WITH_ECL=ON)")
ENDIF()
message(STATUS "========================================")

# -------------------------

# Subdirectories of StreamVorti.
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/doc)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/src/approximants)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/src/support_domain)

# Add Lisp integration module if ECL is enabled
IF(STREAMVORTI_WITH_ECL)
    ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/src/lisp)
ENDIF()


# Compilation of the shared library.
IF( STREAMVORTI_BUILD_SHARED_LIBS )
    IF(STREAMVORTI_WITH_ECL)
        ADD_LIBRARY( ${PROJECT_NAME} SHARED  ${HEADERS}
            $<TARGET_OBJECTS:Approximants>
            $<TARGET_OBJECTS:SupportDomain>
            $<TARGET_OBJECTS:LispIntegration>)
    ELSE()
        ADD_LIBRARY( ${PROJECT_NAME} SHARED  ${HEADERS}
            $<TARGET_OBJECTS:Approximants>
            $<TARGET_OBJECTS:SupportDomain>)
    ENDIF()

    # Link with dependencies
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${CGAL_LIBRARIES} ${MFEM_LIBRARIES})
    IF(STREAMVORTI_WITH_ECL)
        TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${ECL_LIBRARIES})
    ENDIF()

    # Set build folders.
    SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES
                                          LINKER_LANGUAGE CXX
                                          ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib/StreamVorti"
                                          LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib/StreamVorti"
                                          RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/StreamVorti")

    # Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX.
    INSTALL(TARGETS ${PROJECT_NAME} ARCHIVE DESTINATION lib/StreamVorti
                                    LIBRARY DESTINATION lib/StreamVorti
                                    RUNTIME DESTINATION bin/StreamVorti)
ENDIF()


# Compilation of the static library.
IF( STREAMVORTI_BUILD_STATIC_LIBS )
    IF(STREAMVORTI_WITH_ECL)
        ADD_LIBRARY( ${PROJECT_NAME}_static STATIC  ${HEADERS}
            $<TARGET_OBJECTS:Approximants>
            $<TARGET_OBJECTS:SupportDomain>
            $<TARGET_OBJECTS:LispIntegration>)
    ELSE()
        ADD_LIBRARY( ${PROJECT_NAME}_static STATIC  ${HEADERS}
            $<TARGET_OBJECTS:Approximants>
            $<TARGET_OBJECTS:SupportDomain>)
    ENDIF()

    # Link with dependencies
    TARGET_LINK_LIBRARIES(${PROJECT_NAME}_static ${CGAL_LIBRARIES} ${MFEM_LIBRARIES})
    IF(STREAMVORTI_WITH_ECL)
        TARGET_LINK_LIBRARIES(${PROJECT_NAME}_static ${ECL_LIBRARIES})
    ENDIF()

    # Set build folders
    SET_TARGET_PROPERTIES(${PROJECT_NAME}_static PROPERTIES
                                                 LINKER_LANGUAGE CXX
                                                 ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib/StreamVorti"
                                                 LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib/StreamVorti"
                                                 RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/StreamVorti" )

    # Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
    INSTALL(TARGETS ${PROJECT_NAME}_static ARCHIVE DESTINATION lib/StreamVorti
                                           LIBRARY DESTINATION lib/StreamVorti
                                           RUNTIME DESTINATION bin/StreamVorti)

ENDIF()

# MFEM simulation executable.
ADD_EXECUTABLE(MfemRun  mfem_main.cpp)

target_link_libraries(MfemRun CGAL::CGAL CGAL::CGAL_Core)  # Link with CGAL targets

IF(OpenMP_CXX_FOUND) # Link with OPENMP
    message(STATUS "Found OpenMP")
    #TARGET_LINK_LIBRARIES(MfemRun OpenMP::OpenMP_CXX)
ENDIF()

IF(STREAMVORTI_BUILD_SHARED_LIBS)
    TARGET_LINK_LIBRARIES(MfemRun  ${PROJECT_NAME})
ELSE()
    TARGET_LINK_LIBRARIES(MfemRun  ${PROJECT_NAME}_static)
ENDIF()
INSTALL(TARGETS MfemRun RUNTIME DESTINATION bin/StreamVorti)

# StreamVorti executable.
ADD_EXECUTABLE(StreamVorti  streamvorti.cpp)

target_link_libraries(StreamVorti CGAL::CGAL CGAL::CGAL_Core)  # Link with CGAL targets

IF(OpenMP_CXX_FOUND) # Link with OPENMP
    message(STATUS "Found OpenMP")
    #TARGET_LINK_LIBRARIES(StreamVorti OpenMP::OpenMP_CXX)
ENDIF()

IF(STREAMVORTI_BUILD_SHARED_LIBS)
    TARGET_LINK_LIBRARIES(StreamVorti  ${PROJECT_NAME})
ELSE()
    TARGET_LINK_LIBRARIES(StreamVorti  ${PROJECT_NAME}_static)
ENDIF()
INSTALL(TARGETS StreamVorti RUNTIME DESTINATION bin/StreamVorti)

# Install the header files to include/StreamVorti directory under CMAKE_INSTALL_PREFIX.
INSTALL(FILES ${HEADERS} DESTINATION include/StreamVorti)

# Install Lisp headers if ECL is enabled
IF(STREAMVORTI_WITH_ECL)
    INSTALL(FILES
        ${CMAKE_SOURCE_DIR}/include/StreamVorti/lisp/ecl_runtime.hpp
        ${CMAKE_SOURCE_DIR}/include/StreamVorti/lisp/lisp_bridge.hpp
        ${CMAKE_SOURCE_DIR}/include/StreamVorti/lisp/lisp_mesh.hpp
        ${CMAKE_SOURCE_DIR}/include/StreamVorti/lisp/lisp_dcpse.hpp
        ${CMAKE_SOURCE_DIR}/include/StreamVorti/lisp/lisp_loader.hpp
        DESTINATION include/StreamVorti/lisp)

    # Install SDL Lisp source files
    INSTALL(FILES
        ${CMAKE_SOURCE_DIR}/lisp/packages.lisp
        ${CMAKE_SOURCE_DIR}/lisp/geometry.lisp
        ${CMAKE_SOURCE_DIR}/lisp/mesh.lisp
        ${CMAKE_SOURCE_DIR}/lisp/boundaries.lisp
        ${CMAKE_SOURCE_DIR}/lisp/sdl-macros.lisp
        ${CMAKE_SOURCE_DIR}/lisp/simulation.lisp
        DESTINATION share/streamvorti/lisp)

    # Install demo SDL files
    INSTALL(FILES
        ${CMAKE_SOURCE_DIR}/demo/cavity.lisp
        ${CMAKE_SOURCE_DIR}/demo/channel.lisp
        DESTINATION share/streamvorti/demo)
ENDIF()

# -------------------------
# Testing
# -------------------------
IF(STREAMVORTI_BUILD_TESTS)
    enable_testing()
    find_package(GTest REQUIRED)
    add_subdirectory(${CMAKE_SOURCE_DIR}/tests)
ENDIF()

# Install the StreamVorti finding module in the cmake folder under the installation directory.
#INSTALL(FILES cmake/FindStreamVorti.cmake DESTINATION cmake)
