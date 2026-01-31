# tests/CMakeLists.txt
#
# StreamVorti - Unit Tests
# Copyright (C) 2020-2025 Benjamin F. Zwick

include(GoogleTest)

# Test sources
set(TEST_SOURCES
    test_main.cpp
)

# Add Lisp integration tests if ECL is enabled
if(STREAMVORTI_WITH_ECL)
    list(APPEND TEST_SOURCES
        test_ecl_runtime.cpp
        test_lisp_bridge.cpp
        test_lisp_mesh.cpp
        test_lisp_dcpse.cpp
        test_sdl_loader.cpp
    )
endif()

# Create test executable
add_executable(streamvorti_tests ${TEST_SOURCES})

# Include directories
target_include_directories(streamvorti_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${MFEM_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIRS}
)

if(STREAMVORTI_WITH_ECL)
    target_include_directories(streamvorti_tests PRIVATE ${ECL_INCLUDE_DIRS})
endif()

# Link libraries
target_link_libraries(streamvorti_tests
    GTest::gtest
    GTest::gtest_main
    ${MFEM_LIBRARIES}
    ${CGAL_LIBRARIES}
)

if(STREAMVORTI_BUILD_STATIC_LIBS)
    target_link_libraries(streamvorti_tests ${PROJECT_NAME}_static)
else()
    target_link_libraries(streamvorti_tests ${PROJECT_NAME})
endif()

if(STREAMVORTI_WITH_ECL)
    target_link_libraries(streamvorti_tests ${ECL_LIBRARIES})
    target_compile_definitions(streamvorti_tests PRIVATE STREAMVORTI_WITH_ECL)
endif()

# Discover and register tests
gtest_discover_tests(streamvorti_tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    PROPERTIES
        LABELS "unit"
)

# Copy test data files
configure_file(
    ${CMAKE_SOURCE_DIR}/demo/cavity.lisp
    ${CMAKE_CURRENT_BINARY_DIR}/cavity.lisp
    COPYONLY
)

# Code coverage (optional)
option(ENABLE_COVERAGE "Enable code coverage" OFF)
if(ENABLE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(streamvorti_tests PRIVATE --coverage -O0 -g)
    target_link_options(streamvorti_tests PRIVATE --coverage)
endif()
