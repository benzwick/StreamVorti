#----------------------------------------------

# Module source files.
set(SOURCES dcpse.cpp dcpse_2d.cpp dcpse_3d.cpp)

# Add parallel DCPSE if MPI is available
if(MFEM_USE_MPI)
    list(APPEND SOURCES par_dcpse.cpp par_dcpse_2d.cpp)
endif()

# Module header files.
set(HEADERS ${CMAKE_SOURCE_DIR}/include/StreamVorti/approximants/dcpse.hpp
            ${CMAKE_SOURCE_DIR}/include/StreamVorti/approximants/dcpse_2d.hpp
            ${CMAKE_SOURCE_DIR}/include/StreamVorti/approximants/dcpse_3d.hpp)

# Add parallel DCPSE headers if MPI is available
if(MFEM_USE_MPI)
    list(APPEND HEADERS ${CMAKE_SOURCE_DIR}/include/StreamVorti/approximants/par_dcpse.hpp
                        ${CMAKE_SOURCE_DIR}/include/StreamVorti/approximants/par_dcpse_2d.hpp)
endif()

add_library(Approximants OBJECT ${SOURCES} ${HEADERS})

# Link to CGAL and Eigen3 to inherit their include directories
target_link_libraries(Approximants PUBLIC CGAL::CGAL CGAL::CGAL_Core Eigen3::Eigen)

# Enable OpenMP for DCPSE node loop (use parent's OpenMP configuration)
if(OpenMP_CXX_FOUND)
    target_link_libraries(Approximants PRIVATE OpenMP::OpenMP_CXX)
endif()

target_include_directories(Approximants PRIVATE ${CMAKE_SOURCE_DIR}/include)

set_target_properties(Approximants PROPERTIES LINKER_LANGUAGE CXX)

# Install the header files to include/StreamVorti directory under CMAKE_INSTALL_PREFIX.
install(FILES ${HEADERS} DESTINATION include/StreamVorti/approximants)
