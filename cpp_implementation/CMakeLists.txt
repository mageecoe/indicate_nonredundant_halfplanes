cmake_minimum_required(VERSION 3.16)
project(PolytopeRedundancy VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required libraries
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Optional dependencies
find_package(OpenMP QUIET)
find_package(Catch2 3 QUIET)

# Check for qhull availability (optional)
find_program(QHALF_EXECUTABLE qhalf DOC "Path to qhalf executable from qhull")
if(QHALF_EXECUTABLE)
    message(STATUS "Found qhalf: ${QHALF_EXECUTABLE}")
    set(HAVE_QHULL ON)
else()
    message(STATUS "qhalf not found - qhull integration will be disabled")
    set(HAVE_QHULL OFF)
endif()

# Include directories
include_directories(include)

# Source files
set(LIB_SOURCES
    src/matrix_utils.cpp
    src/find_bfs.cpp
    src/active_set_solver.cpp
    src/core.cpp
)

# Add qhull solver if available
if(HAVE_QHULL)
    list(APPEND LIB_SOURCES src/qhull_solver.cpp)
endif()

# Main library
add_library(polytope_redundancy ${LIB_SOURCES})

# Link libraries
target_link_libraries(polytope_redundancy 
    ${BLAS_LIBRARIES} 
    ${LAPACK_LIBRARIES}
)

# OpenMP support if available
if(OpenMP_CXX_FOUND)
    target_link_libraries(polytope_redundancy OpenMP::OpenMP_CXX)
    target_compile_definitions(polytope_redundancy PRIVATE USE_OPENMP)
endif()

# Qhull support if available
if(HAVE_QHULL)
    target_compile_definitions(polytope_redundancy PRIVATE HAVE_QHULL)
    target_compile_definitions(polytope_redundancy PRIVATE QHALF_PATH="${QHALF_EXECUTABLE}")
endif()

# Compiler flags for optimization
target_compile_options(polytope_redundancy PRIVATE
    $<$<CONFIG:Release>:-O3 -DNDEBUG>
    $<$<CONFIG:Debug>:-g -O0>
)

# Handle different BLAS implementations
if(BLAS_FOUND)
    target_compile_definitions(polytope_redundancy PRIVATE ${BLAS_DEFINITIONS})
endif()

# Example executable
add_executable(basic_example examples/basic_example.cpp)
target_link_libraries(basic_example polytope_redundancy)

# Tests if Catch2 is available
if(Catch2_FOUND)
    enable_testing()
    
    set(TEST_SOURCES
        tests/test_main.cpp
        tests/test_utils.cpp
        tests/test_core.cpp
    )
    
    # Add qhull tests if available
    if(HAVE_QHULL)
        list(APPEND TEST_SOURCES tests/test_qhull_comparison.cpp)
    endif()
    
    add_executable(tests ${TEST_SOURCES})
    
    target_link_libraries(tests 
        polytope_redundancy 
        Catch2::Catch2WithMain
    )
    
    include(CTest)
    include(Catch)
    catch_discover_tests(tests)
endif()

# Installation
install(TARGETS polytope_redundancy
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(DIRECTORY include/polytope_redundancy
    DESTINATION include
)