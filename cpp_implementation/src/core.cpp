#include "polytope_redundancy/core.hpp"
#include <chrono>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>

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
    bool should_shift = false;
    const double b_tol = 1e-10;

    Vector Az(m);

    // Check if origin is already interior (matching Matlab: any(b < b_tol))
    bool origin_interior = true;
    for (int i = 0; i < m; ++i) {
        if (b[i] < b_tol) { origin_interior = false; break; }
    }

    if (!origin_interior) {
        // Origin is not interior; a valid interior point z must be provided
        if (z.size() != n) return result;
        A.gemv(z, Az);
        for (int i = 0; i < m; ++i) {
            if (b[i] - Az[i] < b_tol) return result; // z not interior
        }
        should_shift = true;
    } else if (z.size() == n) {
        // Origin is interior; validate z is feasible but do not shift
        A.gemv(z, Az);
        for (int i = 0; i < m; ++i) {
            if (b[i] - Az[i] < b_tol) return result;
        }
    }

    // Shift constraints if needed (only when origin is not interior)
    Matrix A_work = A;
    Vector b_work = b;

    if (should_shift) {
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
    
    // Remove duplicate halfplanes. Like Matlab, rows are kept in their original
    // positions; ind_notdup masks out all but the last copy of each duplicate,
    // and the masked-out copies are reported as redundant.
    Matrix A_combined(m, n + 1);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            A_combined(i, j) = A_work(i, j);
        }
        A_combined(i, n) = b_work[i];
    }
    
    const std::vector<bool> ind_notdup = unique_with_tolerance(A_combined, 1e-5).second;
    
    // Calculate row norms for heuristic
    Vector norms = A_work.row_norms();
    
    // Initialize result arrays
    std::vector<bool> ind_remain(m, false);
    std::vector<bool> ind_nred(m, false);
    std::vector<bool> ind_red(m, false);
    std::vector<bool> ind_failed(m, false);
    std::vector<bool> ind_detected_nred(m, false);
    
    for (int i = 0; i < m; ++i) {
        ind_remain[i] = ind_notdup[i] && ind_to_check[i];
    }
    
    // Active set tracking
    std::vector<bool> ind_active(m, false);
    Vector x;
    
    // Main iteration loop. Every iteration removes at least j from ind_remain,
    // so MAX_IT iterations always suffice to check every halfplane.
    const int MAX_IT = std::count(ind_remain.begin(), ind_remain.end(), true);
    
    for (result.iterations = 0; result.iterations < MAX_IT; ++result.iterations) {
        // Choose halfplane to identify
        int j = -1;
        
        if (x.size() == n) {
            // Use heuristic: choose constraint with largest inner product with current solution
            Vector inner_products(m);
            A_work.gemv(x, inner_products);
            
            double max_cos = -std::numeric_limits<double>::infinity();
            double x_norm = x.norm();
            
            for (int i = 0; i < m; ++i) {
                if (ind_remain[i]) {
                    double cos_val = inner_products[i] / (x_norm * norms[i]);
                    if (cos_val > max_cos) {
                        max_cos = cos_val;
                        j = i;
                    }
                }
            }
        }
        
        if (j == -1) {
            // First iteration, no solution yet, or only NaN cosines (Matlab's
            // max ignores NaN): take the first remaining halfplane
            for (int i = 0; i < m; ++i) {
                if (ind_remain[i]) {
                    j = i;
                    break;
                }
            }
        }
        
        if (j == -1) break;
        
        // Create reduced system excluding duplicates and already identified redundant constraints
        std::vector<int> active_indices;
        for (int i = 0; i < m; ++i) {
            if (ind_notdup[i] && !ind_red[i]) {
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
                A_active(i, k) = A_work(orig_idx, k);
            }
            b_active[i] = b_work[orig_idx];
            ind_active_active[i] = ind_active[orig_idx];
        }
        
        // Solve optimization problem: minimize -A(j,:)'x subject to A_active x <= b_active
        Vector f(n);
        for (int i = 0; i < n; ++i) {
            f[i] = -A_work(j, i);
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
                Ax_j += A_work(j, k) * solver_result.x[k];
            }
            
            if (std::abs(Ax_j - b_work[j]) >= tolerance_) {
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
            x = Vector();
        }
        
        // Update nonredundant and remaining based on detected constraints
        for (int i = 0; i < m; ++i) {
            if (ind_detected_nred[i]) {
                ind_nred[i] = true;
                ind_remain[i] = false;
            }
        }
        
        // Handle symmetry (mirror pairs are (i, i + m/2) in the full indexing)
        if (is_symmetric) {
            const int half = m / 2;
            const int j_mirror = (j + half) % m;

            // One-directional copy for ind_red and ind_failed (matching Matlab)
            ind_red[j_mirror] = ind_red[j];
            ind_failed[j_mirror] = ind_failed[j];

            // Bidirectional OR for ind_nred over all mirror pairs (matching Matlab)
            for (int i = 0; i < half; ++i) {
                bool combined = ind_nred[i] || ind_nred[i + half];
                ind_nred[i] = combined;
                ind_nred[i + half] = combined;
            }

            // One-directional copy for ind_remain, then bidirectional AND (matching Matlab)
            ind_remain[j_mirror] = ind_remain[j];
            for (int i = 0; i < half; ++i) {
                bool combined = ind_remain[i] && ind_remain[i + half];
                ind_remain[i] = combined;
                ind_remain[i + half] = combined;
            }
        }
    }
    
    // Combine results; map back to original order (reverse symmetry permutation)
    result.redundant_indices.assign(m, false);
    result.unverified_indices.assign(m, false);
    for (int i = 0; i < m; ++i) {
        result.redundant_indices[sym_result.permutation[i]] = !(ind_failed[i] || ind_nred[i]);
        result.unverified_indices[sym_result.permutation[i]] = ind_failed[i] && !ind_nred[i];
    }
    
    // Build minimal representation from the original (unshifted, unnormalized) rows
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
    
    auto end_time = std::chrono::high_resolution_clock::now();
    result.solve_time = std::chrono::duration<double>(end_time - start_time).count();
    result.success = true;
    
    return result;
}

} // namespace polytope_redundancy