#pragma once

#include "matrix_utils.hpp"
#include "core.hpp"  // For RedundancyResult
#include <vector>
#include <string>

namespace polytope_redundancy {

/**
 * @brief Qhull-based polytope redundancy removal using qhalf
 * 
 * This class provides an interface to qhull's qhalf program for finding
 * minimal representations of polytopes. It serves as a reference implementation
 * to validate the custom active-set algorithm.
 */
class QhullSolver {
private:
    double tolerance_;
    std::string qhalf_path_;
    
    /**
     * @brief Convert from Ax <= b format to qhull's ax + by + c <= 0 format
     */
    std::pair<Matrix, Vector> convert_to_qhull_format(const Matrix& A, const Vector& b) const;
    
    /**
     * @brief Execute qhalf and parse output
     */
    std::vector<int> run_qhalf(const Matrix& A_qhull, const Vector& b_qhull, const Vector& interior_point) const;
    
    /**
     * @brief Find a valid interior point if none provided
     */
    Vector find_interior_point(const Matrix& A, const Vector& b) const;

public:
    /**
     * @brief Constructor
     * @param tolerance Numerical tolerance (not used by qhull, but kept for consistency)
     * @param qhalf_path Path to qhalf executable (default: "qhalf")
     */
    QhullSolver(double tolerance = 1e-12, const std::string& qhalf_path = "qhalf")
        : tolerance_(tolerance), qhalf_path_(qhalf_path) {}
    
    /**
     * @brief Find minimal representation using qhull
     * @param A Constraint matrix (m x n)
     * @param b Constraint vector (m x 1)
     * @param interior_point Optional interior point (n x 1). If empty, will attempt to find one.
     * @return RedundancyResult with minimal representation
     */
    RedundancyResult indicate_nonredundant_halfplanes(
        const Matrix& A, 
        const Vector& b,
        const Vector& interior_point = Vector()
    );
    
    /**
     * @brief Check if qhalf is available on the system
     */
    bool is_qhalf_available() const;
    
    /**
     * @brief Get the path to qhalf executable
     */
    const std::string& get_qhalf_path() const { return qhalf_path_; }
    
    /**
     * @brief Set the path to qhalf executable
     */
    void set_qhalf_path(const std::string& path) { qhalf_path_ = path; }
};

} // namespace polytope_redundancy