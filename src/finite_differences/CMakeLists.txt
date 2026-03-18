#----------------------------------------------

# Module source files.
set(SOURCES fd.cpp fd_2d.cpp fd_3d.cpp)

# Module header files.
set(HEADERS ${CMAKE_SOURCE_DIR}/include/StreamVorti/finite_differences/fd.hpp
            ${CMAKE_SOURCE_DIR}/include/StreamVorti/finite_differences/fd_2d.hpp
            ${CMAKE_SOURCE_DIR}/include/StreamVorti/finite_differences/fd_3d.hpp)

add_library(FiniteDifferences OBJECT ${SOURCES} ${HEADERS})

target_include_directories(FiniteDifferences PRIVATE ${CMAKE_SOURCE_DIR}/include)

set_target_properties(FiniteDifferences PROPERTIES LINKER_LANGUAGE CXX)

# Install the header files to include/StreamVorti directory under CMAKE_INSTALL_PREFIX.
install(FILES ${HEADERS} DESTINATION include/StreamVorti/finite_differences)
