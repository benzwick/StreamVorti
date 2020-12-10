#--------------------------------------------------------------
CMAKE_MINIMUM_REQUIRED(VERSION 3.1.0 FATAL_ERROR)

IF(POLICY CMP0053)
  CMAKE_POLICY(SET CMP0053 NEW) # CMake 3.1
ENDIF()

PROJECT(StreamVorti VERSION 1.0.0)


# Set project's CMake modules path
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules")

# Library compilation options.
OPTION(STREAMVORTI_BUILD_SHARED_LIBS "Build as shared library" OFF )
OPTION(STREAMVORTI_BUILD_STATIC_LIBS "Build as static library" ON )

# Set the version of C++.
SET(CMAKE_CXX_STANDARD 17)


# Set additional C++ flags for supported compilers.
IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Ob2")
ELSEIF (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mfpmath=sse -O3 -march=native -msse2 -ffast-math -fPIC -Wall")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfpmath=sse -O3 -march=native -msse2 -ffast-math -fPIC -Wall")
ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ansi -frounding-math -fsignaling-nans -mfpmath=sse -march=native -msse2 -O3 -ffast-math -fPIC -Wall")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ansi -frounding-math -fsignaling-nans -mfpmath=sse -march=native -msse2 -O3 -ffast-math -fPIC -Wall")
ENDIF ()

# Include path to the headers of StreamVorti.
INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/include )


# Set the headers of the StreamVorti library.
SET(HEADERS  ${CMAKE_SOURCE_DIR}/include/StreamVorti/stream_vorti.hpp )


# Find dependencies of StreamVorti.
# -------------------------

# Find CGAL.
FIND_PACKAGE(CGAL  REQUIRED COMPONENTS Core)
INCLUDE(${CGAL_USE_FILE})

# Find Eigen
FIND_PACKAGE(Eigen3  REQUIRED)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIRS})

# From https://github.com/GLVis/glvis/blob/9f2df220342d13ef42c2985257ed96ba43f3bb70/CMakeLists.txt
# Import MFEM. The following variables can be used to help CMake find MFEM:
#  * MFEM_DIR - absolute path to the MFEM build or install prefix.
#  * mfem_DIR - absolute path to where MFEMConfig.cmake is.
message(STATUS "Looking for mfem ...")
set(MFEM_DIR "" CACHE PATH "Path to the MFEM build or install prefix.")
if (MFEM_DIR)
   find_package(mfem REQUIRED NAMES MFEM HINTS "${MFEM_DIR}"
                "${MFEM_DIR}/lib/cmake/mfem" NO_DEFAULT_PATH)
else()
   find_package(mfem REQUIRED NAMES MFEM)
endif()
message(STATUS "Found mfem config in: ${mfem_DIR} (version ${MFEM_VERSION})")
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
IF(STREAMVORTI_BUILD_SHARED_LIBS)
    TARGET_LINK_LIBRARIES(MfemRun  ${PROJECT_NAME})
ELSE()
    TARGET_LINK_LIBRARIES(MfemRun  ${PROJECT_NAME}_static)
ENDIF()
INSTALL(TARGETS MfemRun RUNTIME DESTINATION bin/StreamVorti)

# Install the header files to include/StreamVorti directory under CMAKE_INSTALL_PREFIX.
INSTALL(FILES ${HEADERS} DESTINATION include/StreamVorti)

# Install the StreamVorti finding module in the cmake folder under the installation directory.
#INSTALL(FILES cmake/FindStreamVorti.cmake DESTINATION cmake)
