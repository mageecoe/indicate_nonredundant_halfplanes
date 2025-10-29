#include "polytope_redundancy/core.hpp"
#include <chrono>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace polytope_redundancy {

RedundancyResult PolytopeRedundancyRemover::indicate_nonredundant_halfplanes(
    const Matrix& A, const Vector& b,
    const std::vector<bool>& indices_to_check,
    const Vector& interior_point) {
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    RedundancyResult result;
    result.success = false;
    result.iterations = 0;
    
    const int m = A.rows();
    const int n = A.cols();
    
    if (m != b.size()) {
        return result; // Dimension mismatch
    }
    
    // Initialize indices to check
    std::vector<bool> ind_to_check = indices_to_check.size() == m ? 
        indices_to_check : std::vector<bool>(m, true);
    
    // Check if origin is interior point or use provided point
    Vector z = interior_point;
    bool should_shift_back = false;
    const double b_tol = 1e-10;
    
    // Check feasibility
    Vector Az(m);
    if (z.size() == n) {
        A.gemv(z, Az);
        for (int i = 0; i < m; ++i) {
            if (b[i] - Az[i] < b_tol) {
                return result; // Interior point not valid
            }
        }
        should_shift_back = true;
    } else {
        // Check if origin is interior
        for (int i = 0; i < m; ++i) {
            if (b[i] < b_tol) {
                return result; // Origin not interior and no point provided
            }
        }
        z = Vector(n, 0.0);
    }
    
    // Shift constraints if needed
    Matrix A_work = A;
    Vector b_work = b;
    
    if (should_shift_back) {
        for (int i = 0; i < m; ++i) {
            b_work[i] -= Az[i];
        }
    }
    
    // Normalize halfplane description
    auto [A_norm, b_norm] = normalize_halfplane_description(A_work, b_work, false);
    
    // Make set symmetric if possible
    auto sym_result = make_set_symmetric(A_norm, b_norm, 1e-6);
    A_work = sym_result.H;
    b_work = sym_result.h;
    bool is_symmetric = sym_result.is_symmetric;
    
    // Update indices_to_check for symmetry permutation
    std::vector<bool> ind_to_check_perm(m);
    for (int i = 0; i < m; ++i) {
        ind_to_check_perm[i] = ind_to_check[sym_result.permutation[i]];
    }
    ind_to_check = ind_to_check_perm;
    
    // Remove duplicate halfplanes
    Matrix A_combined(m, n + 1);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            A_combined(i, j) = A_work(i, j);
        }
        A_combined(i, n) = b_work[i];
    }
    
    auto [A_unique, ind_notdup] = unique_with_tolerance(A_combined, 1e-5);
    
    // Build reduced system
    int m_unique = A_unique.rows();
    Matrix A_reduced(m_unique, n);
    Vector b_reduced(m_unique);
    std::vector<bool> ind_to_check_reduced;
    
    int unique_idx = 0;
    for (int i = 0; i < m; ++i) {
        if (ind_notdup[i]) {
            for (int j = 0; j < n; ++j) {
                A_reduced(unique_idx, j) = A_work(i, j);
            }
            b_reduced[unique_idx] = b_work[i];
            ind_to_check_reduced.push_back(ind_to_check[i]);
            unique_idx++;
        }
    }
    
    // Calculate row norms for heuristic
    Vector norms = A_reduced.row_norms();
    
    // Initialize result arrays
    std::vector<bool> ind_remain(m_unique, false);
    std::vector<bool> ind_nred(m_unique, false);
    std::vector<bool> ind_red(m_unique, false);
    std::vector<bool> ind_failed(m_unique, false);
    std::vector<bool> ind_detected_nred(m_unique, false);
    
    for (int i = 0; i < m_unique; ++i) {
        ind_remain[i] = ind_to_check_reduced[i];
    }
    
    // Active set tracking
    std::vector<bool> ind_active(m_unique, false);
    Vector x;
    
    // Main iteration loop
    const int MAX_IT = std::count(ind_remain.begin(), ind_remain.end(), true);
    
    for (result.iterations = 0; result.iterations < MAX_IT && result.iterations < max_iterations_; ++result.iterations) {
        // Choose halfplane to identify
        int j = -1;
        
        if (std::any_of(ind_remain.begin(), ind_remain.end(), [](bool b) { return b; })) {
            if (result.iterations == 0 || x.size() == 0) {
                // First iteration or no solution yet
                for (int i = 0; i < m_unique; ++i) {
                    if (ind_remain[i]) {
                        j = i;
                        break;
                    }
                }
            } else {
                // Use heuristic: choose constraint with largest inner product with current solution
                Vector inner_products(m_unique);
                A_reduced.gemv(x, inner_products);
                
                double max_cos = -std::numeric_limits<double>::max();
                double x_norm = x.norm();
                
                for (int i = 0; i < m_unique; ++i) {
                    if (ind_remain[i]) {
                        double cos_val = inner_products[i] / (x_norm * norms[i]);
                        if (cos_val > max_cos) {
                            max_cos = cos_val;
                            j = i;
                        }
                    }
                }
            }
        } else {
            break;
        }
        
        if (j == -1) break;
        
        // Create reduced system excluding already identified redundant constraints
        std::vector<int> active_indices;
        for (int i = 0; i < m_unique; ++i) {
            if (!ind_red[i]) {
                active_indices.push_back(i);
            }
        }
        
        const int m_active = active_indices.size();
        Matrix A_active(m_active, n);
        Vector b_active(m_active);
        std::vector<bool> ind_active_active(m_active, false);
        
        for (int i = 0; i < m_active; ++i) {
            int orig_idx = active_indices[i];
            for (int k = 0; k < n; ++k) {
                A_active(i, k) = A_reduced(orig_idx, k);
            }
            b_active[i] = b_reduced[orig_idx];
            ind_active_active[i] = ind_active[orig_idx];
        }
        
        // Find j in reduced system
        int j_active = -1;
        for (int i = 0; i < m_active; ++i) {
            if (active_indices[i] == j) {
                j_active = i;
                break;
            }
        }
        
        if (j_active == -1) {
            ind_remain[j] = false;
            continue;
        }
        
        // Solve optimization problem: minimize -A(j,:)'x subject to A_active x <= b_active
        Vector f(n);
        for (int i = 0; i < n; ++i) {
            f[i] = -A_reduced(j, i);
        }
        
        ActiveSetResult solver_result = solver_.solve(f, A_active, b_active, x, ind_active_active);
        
        // Update active set for full system
        std::fill(ind_active.begin(), ind_active.end(), false);
        for (int i = 0; i < m_active; ++i) {
            if (solver_result.active_constraints[i]) {
                ind_active[active_indices[i]] = true;
            }
        }
        
        // Update detected nonredundant constraints
        std::fill(ind_detected_nred.begin(), ind_detected_nred.end(), false);
        for (int i = 0; i < m_active; ++i) {
            if (solver_result.detected_nonredundant[i]) {
                ind_detected_nred[active_indices[i]] = true;
            }
        }
        
        // Check if constraint j is redundant
        if (!ind_active[j] && solver_result.optimal_found) {
            // Compute A(j,:) * x
            double Ax_j = 0.0;
            for (int k = 0; k < n; ++k) {
                Ax_j += A_reduced(j, k) * solver_result.x[k];
            }
            
            if (std::abs(Ax_j - b_reduced[j]) >= tolerance_) {
                ind_red[j] = true;
            } else {
                ind_nred[j] = true;
            }
        } else {
            ind_nred[j] = true;
        }
        
        if (!solver_result.optimal_found) {
            ind_failed[j] = true;
        }
        
        // Remove j from remaining constraints
        ind_remain[j] = false;
        
        // Update solution
        x = solver_result.x;
        
        // Check if we need to restart (insufficient active constraints)
        int num_active = std::count(ind_active.begin(), ind_active.end(), true);
        if (num_active < n) {
            std::fill(ind_active.begin(), ind_active.end(), false);
            x = Vector(n, 0.0);
        }
        
        // Update nonredundant and remaining based on detected constraints
        for (int i = 0; i < m_unique; ++i) {
            if (ind_detected_nred[i]) {
                ind_nred[i] = true;
                ind_remain[i] = false;
            }
        }
        
        // Handle symmetry
        if (is_symmetric) {
            const int half = m / 2;
            int j_mirror = (j + half) % m;
            if (j_mirror < m_unique) {
                ind_red[j_mirror] = ind_red[j];
                ind_nred[j_mirror] = ind_nred[j];
                ind_remain[j_mirror] = ind_remain[j];
                ind_failed[j_mirror] = ind_failed[j];
            }
        }
    }
    
    // Combine results
    std::vector<bool> ind_minrep(m_unique);
    std::vector<bool> ind_not_verified(m_unique);
    
    for (int i = 0; i < m_unique; ++i) {
        ind_minrep[i] = ind_failed[i] || ind_nred[i];
        ind_not_verified[i] = ind_failed[i] && !ind_nred[i];
    }
    
    // Map back to original indices
    result.redundant_indices.resize(m, false);
    result.unverified_indices.resize(m, false);
    
    unique_idx = 0;
    for (int i = 0; i < m; ++i) {
        if (ind_notdup[i]) {
            result.redundant_indices[i] = !ind_minrep[unique_idx];
            result.unverified_indices[i] = ind_not_verified[unique_idx];
            unique_idx++;
        }
    }
    
    // Reverse symmetry permutation
    std::vector<bool> redundant_orig(m);
    std::vector<bool> unverified_orig(m);
    for (int i = 0; i < m; ++i) {
        redundant_orig[sym_result.permutation[i]] = result.redundant_indices[i];
        unverified_orig[sym_result.permutation[i]] = result.unverified_indices[i];
    }
    result.redundant_indices = redundant_orig;
    result.unverified_indices = unverified_orig;
    
    // Build minimal representation
    int num_nonredundant = 0;
    for (int i = 0; i < m; ++i) {
        if (!result.redundant_indices[i]) {
            num_nonredundant++;
        }
    }
    
    result.A_min = Matrix(num_nonredundant, n);
    result.b_min = Vector(num_nonredundant);
    
    int min_idx = 0;
    for (int i = 0; i < m; ++i) {
        if (!result.redundant_indices[i]) {
            for (int j = 0; j < n; ++j) {
                result.A_min(min_idx, j) = A(i, j);
            }
            result.b_min[min_idx] = b[i];
            min_idx++;
        }
    }
    
    // Shift back if needed
    if (should_shift_back) {
        Vector Az_min(num_nonredundant);
        result.A_min.gemv(z, Az_min);
        result.b_min.axpy(1.0, Az_min);
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    result.solve_time = std::chrono::duration<double>(end_time - start_time).count();
    result.success = true;
    
    return result;
}

} // namespace polytope_redundancy