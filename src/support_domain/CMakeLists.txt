#----------------------------------------------
CMAKE_MINIMUM_REQUIRED(VERSION 3.30)

IF(POLICY CMP0181)
  CMAKE_POLICY(SET CMP0181 NEW) # CMake 3.30
ENDIF()

PROJECT(SupportDomain)


# Module source files.
SET(SOURCES support_domain.cpp)

# Add parallel support domain if MPI is available
IF(MFEM_USE_MPI)
    LIST(APPEND SOURCES par_support_domain.cpp)
ENDIF()

# Include the headers directory.
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)

# Module header files.
SET(HEADERS ${CMAKE_SOURCE_DIR}/include/StreamVorti/support_domain/support_domain.hpp)

# Add parallel support domain header if MPI is available
IF(MFEM_USE_MPI)
    LIST(APPEND HEADERS ${CMAKE_SOURCE_DIR}/include/StreamVorti/support_domain/par_support_domain.hpp)
ENDIF()

ADD_LIBRARY(${PROJECT_NAME} OBJECT ${SOURCES} ${HEADERS} )

# Link to CGAL to inherit its include directories
TARGET_LINK_LIBRARIES(${PROJECT_NAME} PUBLIC CGAL::CGAL CGAL::CGAL_Core)

# Enable OpenMP for support domain loops (use parent's OpenMP configuration)
IF(OpenMP_CXX_FOUND)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
ENDIF()

SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)

# Install the header files to include/StreamVorti directory under CMAKE_INSTALL_PREFIX.
INSTALL( FILES ${HEADERS} DESTINATION include/StreamVorti/support_domain )
