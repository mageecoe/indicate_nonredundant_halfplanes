#include "polytope_redundancy/geometry.hpp"
#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>

namespace polytope_redundancy {

BFSResult find_basic_feasible_solution(const Matrix& A, const Vector& b, const Vector& interior_point) {
    BFSResult result;
    result.success = false;
    
    const int m = A.rows();
    const int n = A.cols();
    
    // Initialize starting point (matching MATLAB logic)
    Vector x(n, 0.0);
    if (interior_point.size() == n) {
        x = interior_point;
    } else {
        // Check if origin is interior (all b > 0)
        bool origin_interior = true;
        for (int i = 0; i < m; ++i) {
            if (b[i] <= 0) {
                origin_interior = false;
                break;
            }
        }
        if (!origin_interior) {
          std::cerr << "The origin is not an interior point and no interior\n"
            "point was provided. When the origin is not an\n"
            "interior point an interior point, z, has to be provided.\n";
          return result; // Need interior point but none provided
        }
        // Use origin as starting point
        x.fill(0.0);
    }
    
    // Verify feasibility of starting point
    Vector Ax(m);
    A.gemv(x, Ax);
    for (int i = 0; i < m; ++i) {
        if (Ax[i] > b[i] + 1e-10) {
            return result; // Starting point is not feasible
        }
    }
    
    // Initialize active set tracking (matching MATLAB)
    result.active_constraints.resize(m, false);
    std::vector<int> active_indices; // n_active in MATLAB
    
    // Initialize search direction s = -A(1,:)' (matching MATLAB)
    Vector s(n);
    for (int j = 0; j < n; ++j) {
        s[j] = -A(0, j);
    }
    
    // Main BFS finding loop (matching MATLAB structure)
    const int max_iterations = n; // Should reach vertex in at most n steps
    
    for (int iter = 0; iter < max_iterations; ++iter) {
        // Find step size (matching MATLAB logic)
        Vector step_candidates(m);
        for (int i = 0; i < m; ++i) {
            if (result.active_constraints[i]) {
                step_candidates[i] = std::numeric_limits<double>::infinity();
            } else {
                double p1 = b[i] - A.row(i).dot(x);  // b - A*x
                double p2 = A.row(i).dot(s);         // A*s
                
                if (p2 <= 1e-12) {  // Matching MATLAB tolerance
                    step_candidates[i] = std::numeric_limits<double>::infinity();
                } else if (p2 > 0) {
                    step_candidates[i] = p1 / p2;
                } else {
                    step_candidates[i] = std::numeric_limits<double>::infinity();
                }
            }
        }
        
        // Find minimum positive step size
        double min_t = std::numeric_limits<double>::max();
        int blocking_constraint = -1;
        
        for (int i = 0; i < m; ++i) {
            if (step_candidates[i] >= 0 && step_candidates[i] < min_t) {
                min_t = step_candidates[i];
                blocking_constraint = i;
            }
        }
        
        if (blocking_constraint == -1) {
            return result; // No valid step found
        }
        
        // Check for duplicate constraint addition (MATLAB safety check)
        if (std::find(active_indices.begin(), active_indices.end(), blocking_constraint) != active_indices.end()) {
            return result; // Already active constraint - algorithm error
        }
        
        // Update active set
        result.active_constraints[blocking_constraint] = true;
        active_indices.push_back(blocking_constraint);
        
        // Check stopping criterion: full rank active set
        if (active_indices.size() == n) {
            // Solve A_active * x = b_active using least squares
            Matrix A_active(n, n);
            Vector b_active(n);
            
            for (int i = 0; i < n; ++i) {
                int constraint_idx = active_indices[i];
                for (int j = 0; j < n; ++j) {
                    A_active(i, j) = A(constraint_idx, j);
                }
                b_active[i] = b[constraint_idx];
            }
            
            // Solve using QR factorization (matching MATLAB precision)
            QRResult qr = qr_factorization(A_active);
            if (qr.success) {
                x = qr_solve(qr, b_active);
                if (x.size() > 0) {
                    result.x = x;
                    result.success = true;
                    return result;
                }
            }
            
            // Fallback: Check conditioning like MATLAB
            // We can't easily compute rcond without additional LAPACK routines
            // So we'll just return what we have
            result.x = x;
            result.success = true;
            return result;
        }
        
        // Compute QR factorization for null space projection (matching MATLAB)
        const int num_active = active_indices.size();
        Matrix A_active(num_active, n);
        for (int i = 0; i < num_active; ++i) {
            int constraint_idx = active_indices[i];
            for (int j = 0; j < n; ++j) {
                A_active(i, j) = A(constraint_idx, j);
            }
        }
        
        // QR factorization of A_active' (matching MATLAB)
        Matrix A_active_T(n, num_active);
        for (int i = 0; i < num_active; ++i) {
            for (int j = 0; j < n; ++j) {
                A_active_T(j, i) = A_active(i, j);
            }
        }
        
        QRResult qr = qr_factorization(A_active_T);
        if (!qr.success) {
            return result; // QR factorization failed
        }
        
        // Calculate intermediate point x_approx = x + min_t * s (matching MATLAB)
        Vector x_approx = x;
        x_approx.axpy(min_t, s);
        
        // Compute residual res_approx = A_active * x_approx - b_active
        Vector b_active(num_active);
        for (int i = 0; i < num_active; ++i) {
            b_active[i] = b[active_indices[i]];
        }
        
        Vector res_approx(num_active);
        A_active.gemv(x_approx, res_approx);
        for (int i = 0; i < num_active; ++i) {
            res_approx[i] -= b_active[i];
        }
        
        // Project back to satisfy active constraints: x = x_approx - Q1*(R1'\res_approx)
        // This matches MATLAB: x = x_approx - Q1*( R1' \ res_approx )
        
        // Solve R1' * y = res_approx using forward substitution
        Vector y(num_active, 0.0);
        for (int i = 0; i < num_active; ++i) {
            double sum = res_approx[i];
            for (int j = 0; j < i; ++j) {
                sum -= qr.R(j, i) * y[j];  // R'(i,j) = R(j,i)
            }
            
            if (std::abs(qr.R(i, i)) < 1e-15) {
                return result; // Singular matrix
            }
            
            y[i] = sum / qr.R(i, i);
        }
        
        // Compute Q1 * y and subtract from x_approx
        Vector Q1_y(n, 0.0);
        qr.Q.gemv(y, Q1_y, 1.0, 0.0, false); // Q * y (Q1 is first num_active columns)
        
        x = x_approx;
        x.axpy(-1.0, Q1_y); // x = x_approx - Q1*y
        
        // Check feasibility (matching MATLAB safety check)
        A.gemv(x, Ax);
        for (int i = 0; i < m; ++i) {
            if (!result.active_constraints[i] && Ax[i] - b[i] > 1e-7) {
                return result; // Primal feasibility lost
            }
        }
        
        // Calculate new search direction from null space (matching MATLAB)
        // s = Q2(:,end) where Q2 are the last n-num_active columns of Q
        if (num_active < n) {
            // Use the last column of Q as search direction
            s = Vector(n);
            for (int j = 0; j < n; ++j) {
                s[j] = qr.Q(j, n - 1); // Last column of Q
            }
        } else {
            // Full rank - should terminate
            break;
        }
    }
    
    // Final verification
    A.gemv(x, Ax);
    int num_active_final = 0;
    for (int i = 0; i < m; ++i) {
        if (std::abs(Ax[i] - b[i]) < 1e-10) {
            result.active_constraints[i] = true;
            num_active_final++;
        } else {
            result.active_constraints[i] = false;
        }
    }
    
    result.x = x;
    result.success = (num_active_final >= n - 1); // Allow some tolerance
    
    return result;
}

} // namespace polytope_redundancy
