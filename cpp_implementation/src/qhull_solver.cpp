#include "polytope_redundancy/qhull_solver.hpp"
#include <sstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <unistd.h>  // for getpid()
#include <set>       // for removing duplicates

namespace polytope_redundancy {

std::pair<Matrix, Vector> QhullSolver::convert_to_qhull_format(const Matrix& A, const Vector& b) const {
    // Qhull uses ax + by + c <= 0 format, we have Ax <= b
    // Convert: Ax <= b  =>  Ax - b <= 0  =>   A, -b  format for qhull
    
    int m = A.rows();
    int n = A.cols();
    
    Matrix A_qhull(m, n);
    Vector b_qhull(m);
    
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            A_qhull(i, j) = A(i, j);  // Keep A the same
        }
        b_qhull[i] = -b[i];  // Negate b for qhull format
    }
    
    return std::make_pair(std::move(A_qhull), std::move(b_qhull));
}

Vector QhullSolver::find_interior_point(const Matrix& A, const Vector& b) const {
    // Simple approach: try origin first
    int n = A.cols();
    Vector candidate(n, 0.0);
    
    bool is_interior = true;
    for (int i = 0; i < A.rows(); ++i) {
        double dot_product = 0.0;
        for (int j = 0; j < n; ++j) {
            dot_product += A(i, j) * candidate[j];
        }
        if (dot_product >= b[i] - tolerance_) {
            is_interior = false;
            break;
        }
    }
    
    if (is_interior) {
        return candidate;
    }
    
    // If origin doesn't work, try to find a point using a simple heuristic
    // This is a simplified approach - in practice, you might want to solve an LP
    Vector point(n);
    for (int j = 0; j < n; ++j) {
        point[j] = 0.1;  // Small positive values
    }
    
    // Check if this works
    is_interior = true;
    for (int i = 0; i < A.rows(); ++i) {
        double dot_product = 0.0;
        for (int j = 0; j < n; ++j) {
            dot_product += A(i, j) * point[j];
        }
        if (dot_product >= b[i] - tolerance_) {
            is_interior = false;
            break;
        }
    }
    
    if (is_interior) {
        return point;
    }
    
    // If that doesn't work either, return empty vector to signal failure
    return Vector();
}

std::vector<int> QhullSolver::run_qhalf(const Matrix& A_qhull, const Vector& b_qhull, const Vector& interior_point) const {
    // Create input for qhalf
    std::ostringstream input;
    
    int n = A_qhull.cols();  // dimension
    int m = A_qhull.rows();  // number of halfspaces
    
    // First line: dimension 1 interior_point
    input << n << " 1";
    for (int j = 0; j < n; ++j) {
        input << " " << interior_point[j];
    }
    input << "\n";
    
    // Second line: (dimension+1) number_of_halfspaces
    input << (n + 1) << " " << m << "\n";
    
    // Halfspace definitions: coefficients offset
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            input << A_qhull(i, j) << " ";
        }
        input << b_qhull[i] << "\n";
    }
    
    std::string input_str = input.str();
    
    // Use temporary file approach
    std::string temp_filename = "/tmp/qhull_input_" + std::to_string(getpid()) + ".txt";
    FILE* temp_file = fopen(temp_filename.c_str(), "w");
    if (!temp_file) {
        throw std::runtime_error("Failed to create temporary file");
    }
    
    fwrite(input_str.c_str(), 1, input_str.length(), temp_file);
    fclose(temp_file);
    
    // Run qhalf with input redirection
    std::string full_command = qhalf_path_ + " Fx < " + temp_filename;
    FILE* qhalf_pipe = popen(full_command.c_str(), "r");
    if (!qhalf_pipe) {
        remove(temp_filename.c_str());
        throw std::runtime_error("Failed to run qhalf command");
    }
    
    // Read output
    std::vector<int> nonredundant_indices;
    char buffer[256];
    std::string all_output;
    while (fgets(buffer, sizeof(buffer), qhalf_pipe) != nullptr) {
        all_output += buffer;
    }
    
    int exit_code = pclose(qhalf_pipe);
    remove(temp_filename.c_str());
    
    if (exit_code != 0) {
        throw std::runtime_error("qhalf failed with exit code: " + std::to_string(exit_code));
    }
    
    // Parse the output - qhalf Fx returns indices of non-redundant halfspaces
    std::istringstream iss(all_output);
    std::string line;
    while (std::getline(iss, line)) {
        line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1); // trim
        if (line.empty()) continue;
        
        try {
            int index = std::stoi(line);
            if (index >= 0 && index < m) {
                nonredundant_indices.push_back(index);
            }
        } catch (const std::exception&) {
            // Skip non-integer lines
        }
    }
    
    return nonredundant_indices;
}

bool QhullSolver::is_qhalf_available() const {
    std::string command = "which " + qhalf_path_ + " > /dev/null 2>&1";
    int result = system(command.c_str());
    return result == 0;
}

RedundancyResult QhullSolver::indicate_nonredundant_halfplanes(
    const Matrix& A, 
    const Vector& b,
    const Vector& interior_point) {
    
    RedundancyResult result;
    result.success = false;
    result.iterations = 0;
    result.solve_time = 0.0;
    
    int m = A.rows();
    int n = A.cols();
    
    if (m == 0 || n == 0) {
        result.success = true;
        result.A_min = Matrix(0, n);
        result.b_min = Vector(0);
        result.redundant_indices.resize(m, false);
        result.unverified_indices.resize(m, false);
        return result;
    }
    
    // Check if qhalf is available
    if (!is_qhalf_available()) {
        throw std::runtime_error("qhalf not found in PATH. Please install qhull.");
    }
    
    try {
        // Find interior point if not provided
        Vector interior = interior_point;
        if (interior.size() == 0) {
            interior = find_interior_point(A, b);
            if (interior.size() == 0) {
                throw std::runtime_error("Could not find a suitable interior point");
            }
        }
        
        // Validate interior point
        bool is_valid_interior = true;
        for (int i = 0; i < m; ++i) {
            double dot_product = 0.0;
            for (int j = 0; j < n; ++j) {
                dot_product += A(i, j) * interior[j];
            }
            if (dot_product >= b[i] - tolerance_) {
                is_valid_interior = false;
                break;
            }
        }
        
        if (!is_valid_interior) {
            throw std::runtime_error("Provided interior point is not in the interior of the polytope");
        }
        
        // Convert to qhull format
        auto [A_qhull, b_qhull] = convert_to_qhull_format(A, b);
        
        // Run qhalf
        std::vector<int> nonredundant_indices = run_qhalf(A_qhull, b_qhull, interior);
        
        // Initialize result vectors
        result.redundant_indices.resize(m, true);  // Assume all redundant initially
        result.unverified_indices.resize(m, false);
        
        // Mark non-redundant constraints
        std::set<int> unique_indices(nonredundant_indices.begin(), nonredundant_indices.end());
        for (int idx : unique_indices) {
            if (idx >= 0 && idx < m) {
                result.redundant_indices[idx] = false;
            }
        }
        
        // Build minimal representation
        int num_nonredundant = unique_indices.size();
        result.A_min = Matrix(num_nonredundant, n);
        result.b_min = Vector(num_nonredundant);
        
        int out_row = 0;
        for (int i = 0; i < m; ++i) {
            if (!result.redundant_indices[i]) {
                for (int j = 0; j < n; ++j) {
                    result.A_min(out_row, j) = A(i, j);
                }
                result.b_min[out_row] = b[i];
                ++out_row;
            }
        }
        
        result.success = true;
        
    } catch (const std::exception& e) {
        std::cerr << "QhullSolver error: " << e.what() << std::endl;
        // Return original polytope on failure
        result.A_min = A;
        result.b_min = b;
        result.redundant_indices.resize(m, false);
        result.unverified_indices.resize(m, false);
        result.success = false;
    }
    
    return result;
}

} // namespace polytope_redundancy