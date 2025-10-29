#include "polytope_redundancy/matrix_utils.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace polytope_redundancy {

std::pair<Matrix, Vector> normalize_halfplane_description(const Matrix& H, const Vector& h, bool normalize_H) {
    Matrix H_norm = H;
    Vector h_norm = h;
    
    if (normalize_H) {
        // Normalize by row norms
        Vector row_norms = H.row_norms();
        
        for (int i = 0; i < H.rows(); ++i) {
            if (row_norms[i] > 1e-10) {
                const double norm_inv = 1.0 / row_norms[i];
                for (int j = 0; j < H.cols(); ++j) {
                    H_norm(i, j) *= norm_inv;
                }
                h_norm[i] *= norm_inv;
            }
        }
    } else {
        // Normalize by absolute value of h
        for (int i = 0; i < h.size(); ++i) {
            if (std::abs(h[i]) > 1e-10 && !std::isinf(h[i])) {
                const double abs_h = std::abs(h[i]);
                const double scale = 1.0 / abs_h;
                
                for (int j = 0; j < H.cols(); ++j) {
                    H_norm(i, j) *= scale;
                }
                h_norm[i] = h[i] / abs_h;
            }
        }
    }
    
    return std::make_pair(std::move(H_norm), std::move(h_norm));
}

SymmetricResult make_set_symmetric(const Matrix& H, const Vector& h, double tol) {
    SymmetricResult result;
    result.is_symmetric = false;
    result.H = H;
    result.h = h;
    
    const int nrows = h.size();
    
    // Initialize permutation as identity (MATLAB line 33)
    result.permutation.resize(nrows);
    for (int i = 0; i < nrows; ++i) {
        result.permutation[i] = i;  // C++ uses 0-based indexing unlike MATLAB's 1-based
    }
    
    // Check if even number of rows (MATLAB line 34)
    if (nrows % 2 != 0) {
        return result; // Odd number of rows, cannot be symmetric
    }
    
    // Normalize first (MATLAB line 35)
    auto [H_norm, h_norm] = normalize_halfplane_description(H, h, false);
    
    // Create combined matrix [H h] (MATLAB line 36)
    Matrix Hh(nrows, H.cols() + 1);
    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < H.cols(); ++j) {
            Hh(i, j) = H_norm(i, j);
        }
        Hh(i, H.cols()) = h_norm[i];
    }
    
    // Sort rows and track permutation (MATLAB line 38: [Hhs, J] = sortrows(Hh))
    std::vector<int> J(nrows);
    for (int i = 0; i < nrows; ++i) {
        J[i] = i;
    }
    
    std::sort(J.begin(), J.end(), [&](int a, int b) {
        for (int j = 0; j <= H.cols(); ++j) {
            if (std::abs(Hh(a, j) - Hh(b, j)) > 1e-12) {
                return Hh(a, j) < Hh(b, j);
            }
        }
        return false;  // Equal rows
    });
    
    // Create sorted matrix Hhs (MATLAB style)
    Matrix Hhs(nrows, H.cols() + 1);
    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j <= H.cols(); ++j) {
            Hhs(i, j) = Hh(J[i], j);
        }
    }
    
    // Flip lower part (MATLAB lines 40-43)
    const int half = nrows / 2;
    std::vector<int> upper_indices, lower_indices;
    for (int i = 0; i < half; ++i) {
        upper_indices.push_back(i);
        lower_indices.push_back(half + i);
    }
    
    // Flip lower part: Hhs(lower,:) = flipud(Hhs(lower, :))
    for (int i = 0; i < half; ++i) {
        int src_idx = nrows - 1 - i;   // flipud source index
        int dst_idx = half + i;        // destination index
        
        // Swap rows in Hhs
        for (int j = 0; j <= H.cols(); ++j) {
            std::swap(Hhs(dst_idx, j), Hhs(src_idx, j));
        }
        // Swap corresponding indices in J
        std::swap(J[dst_idx], J[src_idx]);
    }
    
    // Check if lower part is symmetric with upper part (MATLAB lines 45-49)
    Matrix Hs_upper(half, H.cols());
    Vector hs_upper(half);
    Matrix Hs_lower(half, H.cols());
    Vector hs_lower(half);
    
    for (int i = 0; i < half; ++i) {
        // Extract upper part
        for (int j = 0; j < H.cols(); ++j) {
            Hs_upper(i, j) = Hhs(i, j);
        }
        hs_upper[i] = Hhs(i, H.cols());
        
        // Extract lower part  
        for (int j = 0; j < H.cols(); ++j) {
            Hs_lower(i, j) = Hhs(half + i, j);
        }
        hs_lower[i] = Hhs(half + i, H.cols());
    }
    
    // Check symmetry: all(all(abs([Hs_upper, hs_upper] - [-Hs_lower, hs_lower]) < tol))
    bool symmetric = true;
    for (int i = 0; i < half && symmetric; ++i) {
        // Check H matrix symmetry: Hs_upper == -Hs_lower
        for (int j = 0; j < H.cols(); ++j) {
            if (std::abs(Hs_upper(i, j) - (-Hs_lower(i, j))) > tol) {
                symmetric = false;
                break;
            }
        }
        if (symmetric) {
            // Check h vector symmetry: hs_upper == hs_lower
            if (std::abs(hs_upper[i] - hs_lower[i]) > tol) {
                symmetric = false;
            }
        }
    }
    
    if (symmetric) {
        // Force set to be numerically symmetric (MATLAB lines 52-54)
        result.is_symmetric = true;
        result.permutation = J;  // Send permutation information 
        
        result.H = Matrix(nrows, H.cols());
        result.h = Vector(nrows);
        
        // H = [Hs_upper; -Hs_upper]
        for (int i = 0; i < half; ++i) {
            for (int j = 0; j < H.cols(); ++j) {
                result.H(i, j) = Hs_upper(i, j);
                result.H(half + i, j) = -Hs_upper(i, j);
            }
        }
        
        // h = [hs_upper; hs_upper]
        for (int i = 0; i < half; ++i) {
            result.h[i] = hs_upper[i];
            result.h[half + i] = hs_upper[i];
        }
    }
    
    return result;
}

std::pair<Matrix, std::vector<bool>> unique_with_tolerance(const Matrix& A, double tol) {
    const int m = A.rows();
    const int n = A.cols();
    
    // Initialize all rows as unique (MATLAB approach)
    std::vector<bool> ind_nred(m, true);
    std::vector<bool> rows_to_check(m, true);
    
    // Identify duplicate rows (matching MATLAB's unique_tol algorithm exactly)
    for (int it = 0; it < m; ++it) {
        
        // Find next row to check (MATLAB lines 41-46)
        int j = -1;
        for (int i = 0; i < m; ++i) {
            if (rows_to_check[i]) {
                j = i;
                rows_to_check[i] = false;
                break;
            }
        }
        
        if (j == -1) break;
        
        // Indices of possible duplicate rows (MATLAB line 49)
        std::vector<bool> ind_duplicate_candidates = rows_to_check;
        
        // Loop over columns (MATLAB line 52)
        for (int col = 0; col < n; ++col) {
            
            // Check stopping criterion (MATLAB lines 54-57)
            bool any_candidates = false;
            for (int i = 0; i < m; ++i) {
                if (ind_duplicate_candidates[i]) {
                    any_candidates = true;
                    break;
                }
            }
            if (!any_candidates) break;
            
            // Check column differences (MATLAB lines 59-67)
            for (int i = 0; i < m; ++i) {
                if (ind_duplicate_candidates[i]) {
                    double diff = std::abs(A(i, col) - A(j, col));
                    if (diff > tol) {
                        ind_duplicate_candidates[i] = false;
                    }
                }
            }
        }
        
        // Mark as non-unique if duplicates found (MATLAB lines 71-73)
        bool found_duplicate = false;
        for (int i = 0; i < m; ++i) {
            if (ind_duplicate_candidates[i]) {
                found_duplicate = true;
                break;
            }
        }
        
        if (found_duplicate) {
            ind_nred[j] = false;  // Mark current row as non-unique (keep later occurrences)
        }
    }
    
    // Count unique rows
    int unique_count = 0;
    for (bool unique : ind_nred) {
        if (unique) unique_count++;
    }
    
    // Create result matrix with unique rows (MATLAB line 77)
    Matrix result(unique_count, n);
    int result_row = 0;
    for (int i = 0; i < m; ++i) {
        if (ind_nred[i]) {
            for (int j = 0; j < n; ++j) {
                result(result_row, j) = A(i, j);
            }
            result_row++;
        }
    }
    
    return std::make_pair(std::move(result), ind_nred);
}

Vector row_norms(const Matrix& A) {
    return A.row_norms();
}

QRResult qr_factorization(const Matrix& A) {
    QRResult result;
    result.success = false;
    
    const int m = A.rows();
    const int n = A.cols();
    const int min_mn = std::min(m, n);
    
    // Copy A since LAPACK modifies the input matrix
    result.Q = A; // Will be overwritten with Q
    result.R = Matrix(min_mn, n); // R is min(m,n) x n
    result.tau = Vector(min_mn);
    
    // Workspace query
    int lwork = -1;
    double work_query;
    int info = 0;
    
    dgeqrf_(&m, &n, result.Q.data(), &m, result.tau.data(), &work_query, &lwork, &info);
    
    if (info != 0) {
        return result;
    }
    
    lwork = static_cast<int>(work_query);
    Vector work(lwork);
    
    // Actual QR factorization
    dgeqrf_(&m, &n, result.Q.data(), &m, result.tau.data(), work.data(), &lwork, &info);
    
    if (info != 0) {
        return result;
    }
    
    // Extract R matrix (upper triangular part)
    for (int i = 0; i < min_mn; ++i) {
        for (int j = i; j < n; ++j) {
            result.R(i, j) = result.Q(i, j);
        }
    }
    
    // Generate explicit Q matrix
    dorgqr_(&m, &min_mn, &min_mn, result.Q.data(), &m, result.tau.data(), work.data(), &lwork, &info);
    
    result.success = (info == 0);
    return result;
}

Vector qr_solve(const QRResult& qr, const Vector& b) {
    if (!qr.success) {
        return Vector(); // Empty vector indicates failure
    }
    
    const int m = qr.Q.rows();
    const int n = qr.R.cols();
    const int min_mn = std::min(m, n);
    
    if (b.size() != m) {
        return Vector(); // Dimension mismatch
    }
    
    // Compute Q'*b
    Vector Qtb(min_mn);
    Vector b_copy = b;
    qr.Q.gemv(b_copy, Qtb, 1.0, 0.0, true); // Q'*b
    
    // Solve R*x = Q'*b using back substitution
    Vector x(n, 0.0);
    
    for (int i = min_mn - 1; i >= 0; --i) {
        double sum = Qtb[i];
        for (int j = i + 1; j < n; ++j) {
            sum -= qr.R(i, j) * x[j];
        }
        
        if (std::abs(qr.R(i, i)) < 1e-15) {
            return Vector(); // Singular matrix
        }
        
        x[i] = sum / qr.R(i, i);
    }
    
    return x;
}

QRResult qr_update(const QRResult& qr_old, const Vector& new_row, bool add_row) {
    // Simplified implementation: just recompute QR factorization
    // For optimal performance, would implement proper QR updates (Givens rotations)
    
    if (!add_row || !qr_old.success) {
        return QRResult{}; // Not implemented for row removal or invalid input
    }
    
    const int m_old = qr_old.Q.rows();
    const int n = qr_old.Q.cols();
    
    if (new_row.size() != qr_old.R.cols()) {
        return QRResult{}; // Dimension mismatch
    }
    
    // Build new matrix with added row
    Matrix A_new(m_old + 1, qr_old.R.cols());
    
    // Copy original matrix (reconstruct from Q*R, but this is inefficient)
    // For now, we'll assume the caller provides the full matrix
    // In practice, this would be optimized with proper QR updates
    
    return QRResult{}; // Placeholder - would need full implementation
}

} // namespace polytope_redundancy