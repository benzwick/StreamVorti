# CMakeLists.txt - Main selector
cmake_minimum_required(VERSION 3.18)

# Option to force which configuration to use
set(FORCE_CONFIG "AUTO" CACHE STRING "Force configuration: AUTO, HPC, or LOCAL")

# Auto-detect or use forced configuration
if(FORCE_CONFIG STREQUAL "HPC") # cmake .. -DFORCE_CONFIG=HPC
    message(STATUS "Forced HPC configuration")
    include(${CMAKE_SOURCE_DIR}/cmake/CMakeLists_HPC.cmake)
elseif(FORCE_CONFIG STREQUAL "LOCAL") # cmake .. -DFORCE_CONFIG=HPC
    message(STATUS "Forced LOCAL configuration")
    include(${CMAKE_SOURCE_DIR}/cmake/CMakeLists_Local.cmake)
else()
    # Auto-detect
    if(DEFINED ENV{SPACK_ENV} OR DEFINED ENV{SPACK_ROOT})
        message(STATUS "Auto-detected Spack/HPC environment")
        include(${CMAKE_SOURCE_DIR}/cmake/CMakeLists_HPC.cmake)
    else()
        message(STATUS "Auto-detected local environment")
        include(${CMAKE_SOURCE_DIR}/cmake/CMakeLists_Local.cmake)
    endif()
endif()