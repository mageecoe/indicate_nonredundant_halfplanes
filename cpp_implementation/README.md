# C++ Implementation of Polytope Redundancy Removal

This directory contains a high-performance C++ implementation of the polytope redundancy removal algorithm originally developed in MATLAB. The algorithm identifies and removes redundant halfplane constraints from polytope H-representations using a primal active-set approach.

## Overview

The implementation translates the MATLAB code from the main directory into modern C++ with BLAS/LAPACK acceleration for optimal performance. It maintains the same algorithmic approach while providing a clean, object-oriented interface.

## Features

- **High Performance**: Uses optimized BLAS/LAPACK routines for all linear algebra operations
- **Modern C++**: Written in C++17 with RAII, move semantics, and exception safety
- **Numerical Robustness**: Same tolerances and numerical safeguards as the original MATLAB implementation
- **Comprehensive Testing**: Unit tests covering all components with numerical accuracy verification
- **Memory Efficient**: Minimizes allocations through careful memory management and in-place operations

## Dependencies

### Required
- **C++17 compatible compiler** (GCC 7+, Clang 6+, MSVC 2017+)
- **CMake 3.16+** for building
- **BLAS library** (OpenBLAS recommended, Intel MKL, or system BLAS)
- **LAPACK library**

### Optional
- **Catch2** for running tests
- **OpenMP** for parallelization (automatic detection)

## Building

```bash
# Create build directory
mkdir build && cd build

# Configure with CMake
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build
make -j$(nproc)

# Run tests (if Catch2 is available)
ctest

# Run example
./basic_example
```

### Build Options

```bash
# Debug build with debugging symbols
cmake .. -DCMAKE_BUILD_TYPE=Debug

# Specify BLAS/LAPACK libraries explicitly
cmake .. -DBLAS_LIBRARIES=/path/to/blas -DLAPACK_LIBRARIES=/path/to/lapack

# Use Intel MKL
cmake .. -DBLA_VENDOR=Intel10_64lp
```

## Usage

### Basic Example

```cpp
#include "polytope_redundancy/core.hpp"

using namespace polytope_redundancy;

// Create constraint matrix A and vector b for polytope {x | Ax <= b}
Matrix A(4, 2);
A(0, 0) = -1.0; A(0, 1) =  1.0;  // -x + y <= 1
A(1, 0) =  1.0; A(1, 1) =  0.0;  //  x     <= 2  (redundant)
A(2, 0) =  0.0; A(2, 1) = -1.0;  //     -y <= 0.3
A(3, 0) =  1.0; A(3, 1) =  0.0;  //  x     <= 1

Vector b(4);
b[0] = 1.0; b[1] = 2.0; b[2] = 0.3; b[3] = 1.0;

// Remove redundant constraints
PolytopeRedundancyRemover solver;
auto result = solver.indicate_nonredundant_halfplanes(A, b);

if (result.success) {
    std::cout << "Original constraints: " << A.rows() << std::endl;
    std::cout << "Non-redundant: " << result.A_min.rows() << std::endl;
    std::cout << "Eliminated: " << A.rows() - result.A_min.rows() << std::endl;
}
```

### Advanced Usage

```cpp
// Specify interior point for numerical stability
Vector interior_point(2);
interior_point[0] = 0.5; interior_point[1] = 0.5;

// Check only specific constraints
std::vector<bool> indices_to_check(4, false);
indices_to_check[1] = true;  // Only check constraint 1

auto result = solver.indicate_nonredundant_halfplanes(A, b, indices_to_check, interior_point);

// Access detailed results
for (int i = 0; i < A.rows(); ++i) {
    if (result.redundant_indices[i]) {
        std::cout << "Constraint " << i << " is redundant" << std::endl;
    } else if (result.unverified_indices[i]) {
        std::cout << "Constraint " << i << " could not be verified" << std::endl;
    }
}
```

## Architecture

### Core Components

- **`Matrix` and `Vector` classes**: RAII wrappers around BLAS-compatible arrays with column-major storage
- **`PolytopeRedundancyRemover`**: Main algorithm implementation
- **`ActiveSetSolver`**: Linear programming solver using primal active-set method  
- **Matrix utilities**: Normalization, symmetry detection, duplicate removal
- **Geometry utilities**: Basic feasible solution finding via rayshots

### Algorithm Flow

1. **Input validation** and interior point verification
2. **Preprocessing**: Halfplane normalization, symmetry detection, duplicate removal
3. **Iterative redundancy checking**: Uses active-set LP solver with constraint selection heuristics
4. **Result assembly**: Maps results back to original constraint indices

### Performance Optimizations

- **BLAS acceleration**: All matrix operations use optimized BLAS routines
- **Memory efficiency**: Minimizes allocations through object reuse and in-place operations
- **Numerical stability**: Careful tolerance management and condition number monitoring
- **Symmetry exploitation**: Reduces computation for symmetric polytopes

## Testing

The test suite covers:

- **Unit tests**: Individual component functionality
- **Integration tests**: End-to-end algorithm verification
- **Numerical accuracy**: Comparison against known results
- **Edge cases**: Degenerate polytopes, empty sets, numerical limits

Run tests with:
```bash
cd build
ctest --verbose
```

## Performance

Typical performance on modern hardware:
- **Small problems** (< 100 constraints): Sub-millisecond
- **Medium problems** (100-1000 constraints): Milliseconds to seconds  
- **Large problems** (1000+ constraints): Seconds to minutes

Performance scales roughly as O(m²n) where m is the number of constraints and n is the dimension.

## Numerical Considerations

The implementation uses the same tolerances as the original MATLAB code:
- **Redundancy detection**: 1e-12
- **Interior point validation**: 1e-10  
- **Duplicate removal**: 1e-5
- **Symmetry detection**: 1e-6

These can be adjusted via the solver constructor for specific applications.

## License

This C++ implementation inherits the same MIT license as the original MATLAB code. See the COPYRIGHT.txt file in the parent directory for full license terms.