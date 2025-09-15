#----------------------------------------------
CMAKE_MINIMUM_REQUIRED(VERSION 3.30)

IF(POLICY CMP0181)
  CMAKE_POLICY(SET CMP0181 NEW) # CMake 3.30
ENDIF()

PROJECT(Approximants)


# Module source files.
SET(SOURCES dcpse.cpp dcpse_2d.cpp dcpse_3d.cpp)

# Include the headers directory.
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)

# Module header files.
SET(HEADERS ${CMAKE_SOURCE_DIR}/include/StreamVorti/approximants/dcpse.hpp
            ${CMAKE_SOURCE_DIR}/include/StreamVorti/approximants/dcpse_2d.hpp
            ${CMAKE_SOURCE_DIR}/include/StreamVorti/approximants/dcpse_3d.hpp)

ADD_LIBRARY(${PROJECT_NAME}  OBJECT  ${SOURCES}  ${HEADERS})

SET_TARGET_PROPERTIES(${PROJECT_NAME}  PROPERTIES  LINKER_LANGUAGE CXX)

# Install the header files to include/StreamVorti directory under CMAKE_INSTALL_PREFIX.
INSTALL(FILES  ${HEADERS}  DESTINATION  include/StreamVorti/approximants)
