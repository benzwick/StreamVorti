#----------------------------------------------
CMAKE_MINIMUM_REQUIRED(VERSION 3.1.0 FATAL_ERROR)

IF(POLICY CMP0053)
  CMAKE_POLICY(SET CMP0053 NEW) # CMake 3.1
ENDIF()

PROJECT(Conditions)


# Module source files.
SET(SOURCES  conditions_handler.cpp  dirichlet.cpp  ebciem.cpp  load_curve.cpp  loading.cpp)

# Include the headers directory.
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)

# Module header files.
SET(HEADERS ${CMAKE_SOURCE_DIR}/include/StreamVorti/conditions/conditions.hpp
            ${CMAKE_SOURCE_DIR}/include/StreamVorti/conditions/conditions_handler.hpp
            ${CMAKE_SOURCE_DIR}/include/StreamVorti/conditions/dirichlet.hpp
            ${CMAKE_SOURCE_DIR}/include/StreamVorti/conditions/ebciem.hpp
            ${CMAKE_SOURCE_DIR}/include/StreamVorti/conditions/load_curve.hpp
            ${CMAKE_SOURCE_DIR}/include/StreamVorti/conditions/loading.hpp)

ADD_LIBRARY(${PROJECT_NAME} OBJECT ${SOURCES} ${HEADERS} )

SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)

# Install the header files to include/StreamVorti directory under CMAKE_INSTALL_PREFIX.
INSTALL( FILES ${HEADERS} DESTINATION include/StreamVorti/conditions )  
