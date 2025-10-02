#--------------------------------------------------------------
CMAKE_MINIMUM_REQUIRED(VERSION 3.30)

IF(POLICY CMP0181)
  CMAKE_POLICY(SET CMP0181 NEW) # CMake 3.30
ENDIF()

PROJECT(StreamVorti VERSION 2.0.0)


# Set project's CMake modules path
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules")

# Library compilation options.
OPTION(STREAMVORTI_BUILD_SHARED_LIBS "Build as shared library" OFF )
OPTION(STREAMVORTI_BUILD_STATIC_LIBS "Build as static library" ON )
OPTION(ENABLE_PROFILING "Enable gprof profiling" OFF )

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
SET(HEADERS  ${CMAKE_SOURCE_DIR}/include/StreamVorti/stream_vorti.hpp )


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

# Use the same C++ compiler as MFEM. This is needed when MFEM was built using
# an MPI wrapper and we do not have explicitly the MPI compile and link flags.
if (NOT CMAKE_CXX_COMPILER AND MFEM_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER "${MFEM_CXX_COMPILER}")
endif()

# Include paths and libraries needed by MFEM
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MFEM_CXX_FLAGS}")
INCLUDE_DIRECTORIES(${MFEM_INCLUDE_DIRS})

# -------------------------

# Subdirectories of StreamVorti.
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/doc)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/src/approximants)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/src/support_domain)


# Compilation of the shared library.
IF( STREAMVORTI_BUILD_SHARED_LIBS )
    ADD_LIBRARY( ${PROJECT_NAME} SHARED  ${HEADERS}
        $<TARGET_OBJECTS:Approximants>
        $<TARGET_OBJECTS:SupportDomain>)

    # Link with dependencies
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${CGAL_LIBRARIES} ${MFEM_LIBRARIES})

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
    ADD_LIBRARY( ${PROJECT_NAME}_static STATIC  ${HEADERS}
        $<TARGET_OBJECTS:Approximants>
        $<TARGET_OBJECTS:SupportDomain>)

    # Link with dependencies
    TARGET_LINK_LIBRARIES(${PROJECT_NAME}_static ${CGAL_LIBRARIES} ${MFEM_LIBRARIES})

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

# Install the StreamVorti finding module in the cmake folder under the installation directory.
#INSTALL(FILES cmake/FindStreamVorti.cmake DESTINATION cmake)
