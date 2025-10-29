#pragma once

#include "matrix_utils.hpp"
#include <vector>

namespace polytope_redundancy {

struct BFSResult {
    Vector x;
    std::vector<bool> active_constraints;
    bool success;
};

// Find basic feasible solution using iterative rayshots
BFSResult find_basic_feasible_solution(const Matrix& A, const Vector& b, const Vector& interior_point = Vector());

} // namespace polytope_redundancy