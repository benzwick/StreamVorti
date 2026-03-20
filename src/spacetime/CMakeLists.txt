# Space-Time Navier-Stokes module
ADD_LIBRARY(SpaceTime OBJECT
    st_mesh.cpp
    st_integrators.cpp
    st_solver.cpp
)

target_include_directories(SpaceTime PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${MFEM_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIRS}
)

target_link_libraries(SpaceTime
    ${MFEM_LIBRARIES}
)

IF(OpenMP_CXX_FOUND)
    target_link_libraries(SpaceTime OpenMP::OpenMP_CXX)
ENDIF()
