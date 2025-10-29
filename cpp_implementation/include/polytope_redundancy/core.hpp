#pragma once

#include "matrix_utils.hpp"
#include "active_set_solver.hpp"
#include <vector>

namespace polytope_redundancy {

struct RedundancyResult {
    Matrix A_min;
    Vector b_min;
    std::vector<bool> redundant_indices;
    std::vector<bool> unverified_indices;
    bool success;
    int iterations;
    double solve_time;
};

class PolytopeRedundancyRemover {
private:
    double tolerance_;
    int max_iterations_;
    ActiveSetSolver solver_;
    
public:
    PolytopeRedundancyRemover(double tolerance = 1e-12, int max_iterations = 10000)
        : tolerance_(tolerance), max_iterations_(max_iterations), solver_(1e-6, max_iterations) {}  // MATLAB uses 1e-6 for active set solver
    
    // Main algorithm interface
    RedundancyResult indicate_nonredundant_halfplanes(
        const Matrix& A, 
        const Vector& b,
        const std::vector<bool>& indices_to_check = std::vector<bool>(),
        const Vector& interior_point = Vector()
    );
    
    // Convenience function matching MATLAB interface
    std::pair<Matrix, Vector> remove_redundant_constraints(
        const Matrix& A, 
        const Vector& b,
        const Vector& interior_point = Vector()
    ) {
        auto result = indicate_nonredundant_halfplanes(A, b, std::vector<bool>(), interior_point);
        return std::make_pair(result.A_min, result.b_min);
    }
};

} // namespace polytope_redundancy