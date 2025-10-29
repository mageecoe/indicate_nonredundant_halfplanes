#include "polytope_redundancy/active_set_solver.hpp"
#include <cmath>
#include <algorithm>
#include <limits>

namespace polytope_redundancy {

bool ActiveSetSolver::solve_qr_system(const Matrix& A_active, const Vector& f, 
                                      Vector& x, Vector& mu, QRResult& qr_factors) {
    const int n = A_active.cols();
    const int m_active = A_active.rows();
    
    if (m_active == 0) {
        // Unconstrained problem: x = 0 (since we're minimizing f'x with f = -constraint_normal)
        x = Vector(n, 0.0);
        mu = Vector(0);
        qr_factors = QRResult{};
        qr_factors.success = true;
        return true;
    }
    
    // Compute QR factorization of A_active' (matching MATLAB's approach)
    Matrix A_active_T(n, m_active);
    for (int i = 0; i < m_active; ++i) {
        for (int j = 0; j < n; ++j) {
            A_active_T(j, i) = A_active(i, j);
        }
    }
    
    qr_factors = qr_factorization(A_active_T);
    if (!qr_factors.success) {
        return false;
    }
    
    // Calculate dual variables: mu = R\(-Q'*f)  [MATLAB approach]
    // First compute Q'*f
    Vector Qf(m_active);
    qr_factors.Q.gemv(f, Qf, 1.0, 0.0, true); // Q' * f
    
    // Negate: -Q'*f 
    Qf.scale(-1.0);
    
    // Solve R*mu = -Q'*f using back substitution
    mu = Vector(m_active, 0.0);
    for (int i = m_active - 1; i >= 0; --i) {
        double sum = Qf[i];
        for (int j = i + 1; j < m_active; ++j) {
            sum -= qr_factors.R(i, j) * mu[j];
        }
        
        if (std::abs(qr_factors.R(i, i)) < 1e-15) {
            return false; // Singular matrix
        }
        
        mu[i] = sum / qr_factors.R(i, i);
    }
    
    return true;
}

ActiveSetResult ActiveSetSolver::solve(const Vector& f, const Matrix& A, const Vector& b,
                                      const Vector& x_init, const std::vector<bool>& active_init) {
    ActiveSetResult result;
    result.optimal_found = false;
    result.iterations = 0;
    
    const int m = A.rows();
    const int n = A.cols();
    
    // Initialize detected nonredundant constraints (matching MATLAB)
    result.detected_nonredundant.resize(m, false);
    
    // Initialize active set - use BFS if not provided (matching MATLAB logic)
    Vector x = x_init;
    std::vector<bool> active = active_init;
    
    if (x_init.size() != n || std::all_of(active.begin(), active.end(), [](bool a) { return !a; })) {
        BFSResult bfs = find_basic_feasible_solution(A, b);
        if (!bfs.success) {
            return result; // Cannot find initial feasible solution
        }
        x = bfs.x;
        active = bfs.active_constraints;
    }
    
    // Track active constraint indices (matching MATLAB's n_active)
    std::vector<int> n_active;
    for (int i = 0; i < m; ++i) {
        if (active[i]) {
            n_active.push_back(i);
        }
    }
    
    // QR factorization variables (matching MATLAB)
    QRResult Q_R;  // Main QR factorization
    QRResult Q_reduced_R_reduced; // Reduced QR factorization
    std::vector<int> n_reduced; // Indices for reduced system
    Matrix I = Matrix::identity(n);       // Identity matrix for QR updates
    Matrix I_reduced; // Reduced identity matrix - initialize later with bounds check
    
    // Main iteration loop (matching MATLAB structure)
    for (result.iterations = 0; result.iterations < max_iterations_; ++result.iterations) {
        
        // Mark all active constraints as detected nonredundant (MATLAB line 97)
        for (int i = 0; i < m; ++i) {
            if (active[i]) {
                result.detected_nonredundant[i] = true;
            }
        }
        
        // Build active constraint matrix (MATLAB lines 101-103)
        const int m_active = n_active.size();
        Matrix A_active(m_active, n);
        Vector b_active(m_active);
        
        for (int i = 0; i < m_active; ++i) {
            int orig_idx = n_active[i];
            for (int j = 0; j < n; ++j) {
                A_active(i, j) = A(orig_idx, j);
            }
            b_active[i] = b[orig_idx];
        }
        
        // Find QR factorization on first iteration (MATLAB lines 100-104)
        if (result.iterations == 0) {
            Matrix A_active_T(n, m_active);
            for (int i = 0; i < m_active; ++i) {
                for (int j = 0; j < n; ++j) {
                    A_active_T(j, i) = A_active(i, j);
                }
            }
            Q_R = qr_factorization(A_active_T);
            if (!Q_R.success) {
                return result; // QR factorization failed
            }
        }
        
        // Calculate dual variables using QR factors (MATLAB lines 107-109)
        Vector Qf(m_active);
        Q_R.Q.gemv(f, Qf, 1.0, 0.0, true); // Q' * f
        Qf.scale(-1.0); // -Q'*f
        
        Vector mu(m_active, 0.0);
        for (int i = m_active - 1; i >= 0; --i) {
            double sum = Qf[i];
            for (int j = i + 1; j < m_active; ++j) {
                sum -= Q_R.R(i, j) * mu[j];
            }
            
            if (std::abs(Q_R.R(i, i)) < 1e-15) {
                result.optimal_found = false;
                std::fill(active.begin(), active.end(), false);
                result.active_constraints = active;
                return result;
            }
            
            mu[i] = sum / Q_R.R(i, i);
        }
        
        // Check stopping criterion (MATLAB line 112)
        bool all_mu_positive = true;
        for (int i = 0; i < m_active; ++i) {
            if (mu[i] < -tolerance_) {  // Using MATLAB's -1e-6 tolerance
                all_mu_positive = false;
                break;
            }
        }
        
        if (all_mu_positive) {
            result.optimal_found = true;
            break;
        }
        
        // Find and remove constraint with most negative dual variable (MATLAB lines 120-122)
        int n_remove_idx = 0;
        double min_mu = mu[0];
        for (int i = 1; i < m_active; ++i) {
            if (mu[i] < min_mu) {
                min_mu = mu[i];
                n_remove_idx = i;
            }
        }
        
        int n_remove = n_active[n_remove_idx];
        active[n_remove] = false;
        
        // Get constraint to remove (MATLAB line 123)
        Vector a_remove(n);
        for (int j = 0; j < n; ++j) {
            a_remove[j] = A(n_remove, j);
        }
        
        // Update reduced system for null space calculation (MATLAB lines 128-154)
        if (result.iterations == 0) {
            // Find A_reduced (MATLAB lines 129-136)
            Matrix A_reduced_temp(m_active, n);
            n_reduced.clear();
            for (int i = 0; i < m_active; ++i) {
                if (n_active[i] != n_remove) {
                    for (int j = 0; j < n; ++j) {
                        A_reduced_temp(n_reduced.size(), j) = A(n_active[i], j);
                    }
                    n_reduced.push_back(n_active[i]);
                }
            }
            
            if (n_reduced.size() > 0) {
                Matrix A_reduced(n_reduced.size(), n);
                for (int i = 0; i < n_reduced.size(); ++i) {
                    for (int j = 0; j < n; ++j) {
                        A_reduced(i, j) = A_reduced_temp(i, j);
                    }
                }
                
                Matrix A_reduced_T(n, n_reduced.size());
                for (int i = 0; i < n_reduced.size(); ++i) {
                    for (int j = 0; j < n; ++j) {
                        A_reduced_T(j, i) = A_reduced(i, j);
                    }
                }
                Q_reduced_R_reduced = qr_factorization(A_reduced_T);
            }
        } else {
            // Update QR factorization for reduced system if needed (MATLAB lines 138-153)
            // This would require implementing qrupdate functionality
            // For now, recompute (less efficient but correct)
            n_reduced.clear();
            for (int i = 0; i < m_active; ++i) {
                if (n_active[i] != n_remove) {
                    n_reduced.push_back(n_active[i]);
                }
            }
            
            if (n_reduced.size() > 0) {
                Matrix A_reduced(n_reduced.size(), n);
                for (int i = 0; i < n_reduced.size(); ++i) {
                    for (int j = 0; j < n; ++j) {
                        A_reduced(i, j) = A(n_reduced[i], j);
                    }
                }
                
                Matrix A_reduced_T(n, n_reduced.size());
                for (int i = 0; i < n_reduced.size(); ++i) {
                    for (int j = 0; j < n; ++j) {
                        A_reduced_T(j, i) = A_reduced(i, j);
                    }
                }
                Q_reduced_R_reduced = qr_factorization(A_reduced_T);
            }
        }
        
        // Calculate search direction in null space (MATLAB line 157)
        Vector s(n);
        if (n_reduced.size() < n && Q_reduced_R_reduced.success) {
            // Use last column of Q as search direction (null space)
            for (int j = 0; j < n; ++j) {
                s[j] = Q_reduced_R_reduced.Q(j, n - 1);
            }
        } else {
            // Fallback: use negative gradient
            s = f;
            s.scale(-1.0);
        }
        
        // Check direction (MATLAB lines 160-162)
        double a_remove_dot_s = 0.0;
        for (int j = 0; j < n; ++j) {
            a_remove_dot_s += a_remove[j] * s[j];
        }
        if (a_remove_dot_s > 0) {
            s.scale(-1.0);
        }
        
        // Find step sizes to activate inactive constraints (MATLAB lines 169-173)
        Vector p1(m), p2(m);
        A.gemv(x, p1);  // A*x
        A.gemv(s, p2);  // A*s
        
        for (int i = 0; i < m; ++i) {
            p1[i] = b[i] - p1[i];  // b - A*x
        }
        
        Vector step_sizes(m, std::numeric_limits<double>::infinity());
        for (int i = 0; i < m; ++i) {
            if (active[i]) {
                step_sizes[i] = std::numeric_limits<double>::infinity();
            } else if (p2[i] <= 1e-12) {  // MATLAB's tolerance
                step_sizes[i] = std::numeric_limits<double>::infinity();
            } else {
                step_sizes[i] = p1[i] / p2[i];
            }
        }
        
        // Find minimum positive step size (MATLAB line 178)
        double min_t = std::numeric_limits<double>::max();
        int n_add = -1;
        for (int i = 0; i < m; ++i) {
            if (step_sizes[i] >= 0 && step_sizes[i] < min_t) {
                min_t = step_sizes[i];
                n_add = i;
            }
        }
        
        // Check for infinite step (MATLAB lines 181-183)
        if (std::isinf(min_t)) {
            result.optimal_found = false;
            std::fill(active.begin(), active.end(), false);
            result.active_constraints = active;
            return result;
        }
        
        // Check for duplicate constraint (MATLAB lines 185-187)
        if (std::find(n_active.begin(), n_active.end(), n_add) != n_active.end()) {
            result.optimal_found = false;
            std::fill(active.begin(), active.end(), false);
            result.active_constraints = active;
            return result;
        }
        
        // Update active constraints (MATLAB lines 194-199)
        Vector a_add(n);
        for (int j = 0; j < n; ++j) {
            a_add[j] = A(n_add, j);
        }
        double b_add = b[n_add];
        
        // Replace removed constraint with added constraint in n_active
        n_active[n_remove_idx] = n_add;
        active[n_add] = true;
        
        // Update A_active and b_active
        for (int j = 0; j < n; ++j) {
            A_active(n_remove_idx, j) = a_add[j];
        }
        b_active[n_remove_idx] = b_add;
        
        // Update QR factorization (MATLAB line 202)
        // For now, recompute QR factorization (would need qrupdate for efficiency)
        Matrix A_active_T_updated(n, m_active);
        for (int i = 0; i < m_active; ++i) {
            for (int j = 0; j < n; ++j) {
                A_active_T_updated(j, i) = A_active(i, j);
            }
        }
        Q_R = qr_factorization(A_active_T_updated);
        
        if (!Q_R.success) {
            result.optimal_found = false;
            std::fill(active.begin(), active.end(), false);
            result.active_constraints = active;
            return result;
        }
        
        // Check condition number (MATLAB lines 206-210)
        // Approximate condition check by looking at diagonal elements
        double min_diag = std::numeric_limits<double>::max();
        for (int i = 0; i < std::min(Q_R.R.rows(), Q_R.R.cols()); ++i) {
            min_diag = std::min(min_diag, std::abs(Q_R.R(i, i)));
        }
        if (min_diag < 1e-15) {
            result.optimal_found = false;
            std::fill(active.begin(), active.end(), false);
            result.active_constraints = active;
            return result;
        }
        
        // Update solution x = Q * (R' \ b_active) (MATLAB line 213)
        Vector y(m_active, 0.0);
        for (int i = 0; i < m_active; ++i) {
            double sum = b_active[i];
            for (int j = 0; j < i; ++j) {
                sum -= Q_R.R(j, i) * y[j];  // R'(i,j) = R(j,i)
            }
            
            if (std::abs(Q_R.R(i, i)) < 1e-15) {
                result.optimal_found = false;
                std::fill(active.begin(), active.end(), false);
                result.active_constraints = active;
                return result;
            }
            
            y[i] = sum / Q_R.R(i, i);
        }
        
        x = Vector(n, 0.0);
        Q_R.Q.gemv(y, x, 1.0, 0.0, false); // Q * y
        
        // Correct variables if inaccurate (MATLAB lines 216-226)
        Vector Ax_check(m_active);
        A_active.gemv(x, Ax_check);
        double max_res = 0.0;
        for (int i = 0; i < m_active; ++i) {
            max_res = std::max(max_res, std::abs(Ax_check[i] - b_active[i]));
        }
        
        if (max_res > 1e-9) {
            // Try direct solve (fallback)
            // This would require implementing A_active \ b_active
            // For now, keep current solution
        }
        
        // Verify feasibility (MATLAB lines 229-234)
        Vector Ax_full(m);
        A.gemv(x, Ax_full);
        bool feasible = true;
        for (int i = 0; i < m; ++i) {
            if (Ax_full[i] - b[i] > max_res + 1e-7) {
                feasible = false;
                break;
            }
        }
        
        if (!feasible) {
            result.optimal_found = false;
            std::fill(active.begin(), active.end(), false);
            result.active_constraints = active;
            return result;
        }
    }
    
    // Check if maximum iterations reached (MATLAB lines 239-243)
    if (result.iterations >= max_iterations_) {
        result.optimal_found = false;
        std::fill(active.begin(), active.end(), false);
    }
    
    result.x = x;
    result.active_constraints = active;
    
    return result;
}

} // namespace polytope_redundancy