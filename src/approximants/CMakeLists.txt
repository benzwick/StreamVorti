#----------------------------------------------
CMAKE_MINIMUM_REQUIRED(VERSION 3.30)

IF(POLICY CMP0181)
  CMAKE_POLICY(SET CMP0181 NEW) # CMake 3.30
ENDIF()

PROJECT(Approximants)

FIND_PACKAGE(Eigen3 CONFIG REQUIRED)

# Module source files.
SET(SOURCES dcpse.cpp dcpse_2d.cpp dcpse_3d.cpp)

# Add parallel DCPSE if MPI is available
IF(MFEM_USE_MPI)
    LIST(APPEND SOURCES par_dcpse.cpp par_dcpse_2d.cpp)
ENDIF()

# Include the headers directory.
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)

# Module header files.
SET(HEADERS ${CMAKE_SOURCE_DIR}/include/StreamVorti/approximants/dcpse.hpp
            ${CMAKE_SOURCE_DIR}/include/StreamVorti/approximants/dcpse_2d.hpp
            ${CMAKE_SOURCE_DIR}/include/StreamVorti/approximants/dcpse_3d.hpp)

# Add parallel DCPSE headers if MPI is available
IF(MFEM_USE_MPI)
    LIST(APPEND HEADERS ${CMAKE_SOURCE_DIR}/include/StreamVorti/approximants/par_dcpse.hpp
                        ${CMAKE_SOURCE_DIR}/include/StreamVorti/approximants/par_dcpse_2d.hpp)
ENDIF()

ADD_LIBRARY(${PROJECT_NAME}  OBJECT  ${SOURCES}  ${HEADERS})

# Link to CGAL and Eigen3 to inherit their include directories
TARGET_LINK_LIBRARIES(${PROJECT_NAME} PUBLIC CGAL::CGAL CGAL::CGAL_Core Eigen3::Eigen)

# Enable OpenMP for DCPSE node loop (use parent's OpenMP configuration)
IF(OpenMP_CXX_FOUND)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
ENDIF()

SET_TARGET_PROPERTIES(${PROJECT_NAME}  PROPERTIES  LINKER_LANGUAGE CXX)

# Install the header files to include/StreamVorti directory under CMAKE_INSTALL_PREFIX.
INSTALL(FILES  ${HEADERS}  DESTINATION  include/StreamVorti/approximants)
