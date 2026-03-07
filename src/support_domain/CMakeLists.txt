#----------------------------------------------

# Module source files.
set(SOURCES support_domain.cpp)

# Add parallel support domain if MPI is available
if(MFEM_USE_MPI)
    list(APPEND SOURCES par_support_domain.cpp)
endif()

# Module header files.
set(HEADERS ${CMAKE_SOURCE_DIR}/include/StreamVorti/support_domain/support_domain.hpp)

# Add parallel support domain header if MPI is available
if(MFEM_USE_MPI)
    list(APPEND HEADERS ${CMAKE_SOURCE_DIR}/include/StreamVorti/support_domain/par_support_domain.hpp)
endif()

add_library(SupportDomain OBJECT ${SOURCES} ${HEADERS})

# Link to CGAL to inherit its include directories
target_link_libraries(SupportDomain PUBLIC CGAL::CGAL CGAL::CGAL_Core)

# Enable OpenMP for support domain loops (use parent's OpenMP configuration)
if(OpenMP_CXX_FOUND)
    target_link_libraries(SupportDomain PRIVATE OpenMP::OpenMP_CXX)
endif()

target_include_directories(SupportDomain PRIVATE ${CMAKE_SOURCE_DIR}/include)

set_target_properties(SupportDomain PROPERTIES LINKER_LANGUAGE CXX)

# Install the header files to include/StreamVorti directory under CMAKE_INSTALL_PREFIX.
install(FILES ${HEADERS} DESTINATION include/StreamVorti/support_domain)
