#pragma once

#include "matrix_utils.hpp"
#include "geometry.hpp"
#include <vector>

namespace polytope_redundancy {

struct ActiveSetResult {
    Vector x;
    std::vector<bool> active_constraints;
    std::vector<bool> detected_nonredundant;
    bool optimal_found;
    int iterations;
};

class ActiveSetSolver {
private:
    double tolerance_;
    int max_iterations_;
    
    // QR-based solver for active set method (matching MATLAB approach)
    bool solve_qr_system(const Matrix& A_active, const Vector& f, 
                        Vector& x, Vector& mu, QRResult& qr_factors);
    
public:
    ActiveSetSolver(double tolerance = 1e-6, int max_iterations = 10000)  // Changed to match MATLAB
        : tolerance_(tolerance), max_iterations_(max_iterations) {}
    
    // Main solver interface
    ActiveSetResult solve(const Vector& f, const Matrix& A, const Vector& b,
                         const Vector& x_init = Vector(),
                         const std::vector<bool>& active_init = std::vector<bool>());
};

} // namespace polytope_redundancy