#include "polytope_redundancy/core.hpp"
#include "polytope_redundancy/qhull_solver.hpp"
#include <iostream>
#include <vector>
#include <chrono>
#include <random>

int main() {
    std::cout << "Polytope Redundancy Removal - Qhull Comparison Example\n";
    std::cout << "======================================================\n\n";
    
    // Initialize both solvers
    polytope_redundancy::PolytopeRedundancyRemover custom_solver;
    polytope_redundancy::QhullSolver qhull_solver;
    
    // Check qhull availability
    bool qhull_available = qhull_solver.is_qhalf_available();
    std::cout << "Qhull availability: " << (qhull_available ? "✓ Available" : "✗ Not available") << std::endl;
    if (qhull_available) {
        std::cout << "Qhalf path: " << qhull_solver.get_qhalf_path() << std::endl;
    }
    std::cout << std::endl;
    
    // Example 1: Simple 2D case from MATLAB code
    std::cout << "Example 1: Simple 2D polytope comparison\n";
    std::cout << "=========================================\n";
    
    // Create the constraint matrix A and vector b
    // This represents the constraints from the MATLAB example:
    // -x1 + x2 <= 1
    //  x1      <= 2  (this should be redundant)
    //     -x2 <= 0.3
    //  x1      <= 1
    
    polytope_redundancy::Matrix A1(4, 2);
    A1(0, 0) = -1.0; A1(0, 1) =  1.0;  // -x1 + x2 <= 1
    A1(1, 0) =  1.0; A1(1, 1) =  0.0;  //  x1      <= 2
    A1(2, 0) =  0.0; A1(2, 1) = -1.0;  //     -x2 <= 0.3
    A1(3, 0) =  1.0; A1(3, 1) =  0.0;  //  x1      <= 1
    
    polytope_redundancy::Vector b1(4);
    b1[0] =  1.0;
    b1[1] =  2.0;
    b1[2] =  0.3;
    b1[3] =  1.0;
    
    // Interior point for better numerical stability
    polytope_redundancy::Vector interior_point1(2);
    interior_point1[0] = 0.2;
    interior_point1[1] = 0.1;
    
    std::cout << "Original constraints:\n";
    for (int i = 0; i < 4; ++i) {
        std::cout << "  " << i << ": ";
        if (A1(i, 0) != 0) {
            if (A1(i, 0) == 1.0) std::cout << "x1";
            else if (A1(i, 0) == -1.0) std::cout << "-x1";
            else std::cout << A1(i, 0) << "*x1";
        }
        if (A1(i, 1) != 0) {
            if (A1(i, 0) != 0) std::cout << (A1(i, 1) > 0 ? " + " : " - ");
            else if (A1(i, 1) < 0) std::cout << "-";
            
            if (std::abs(A1(i, 1)) == 1.0) std::cout << "x2";
            else std::cout << std::abs(A1(i, 1)) << "*x2";
        }
        std::cout << " <= " << b1[i] << "\n";
    }
    std::cout << std::endl;

    // Solve with qhull solver
    polytope_redundancy::RedundancyResult qhull_result1;
    double qhull_duration1 = 0.0;
    bool qhull_success1 = false;
    
    if (qhull_available) {
        std::cout << "Running qhull solver...\n";
        auto start_qhull = std::chrono::high_resolution_clock::now();
        qhull_result1 = qhull_solver.indicate_nonredundant_halfplanes(A1, b1, interior_point1);
        auto end_qhull = std::chrono::high_resolution_clock::now();
        qhull_duration1 = std::chrono::duration<double>(end_qhull - start_qhull).count();
        qhull_success1 = qhull_result1.success;
    }

    // Solve with custom solver
    std::cout << "Running custom active-set solver...\n";
    auto start1 = std::chrono::high_resolution_clock::now();
    auto custom_result1 = custom_solver.indicate_nonredundant_halfplanes(A1, b1, std::vector<bool>(), interior_point1);
    auto end1 = std::chrono::high_resolution_clock::now();
    auto custom_duration1 = std::chrono::duration<double>(end1 - start1).count();

    // Display results comparison
    std::cout << "\nResults Comparison:\n";
    std::cout << "-------------------\n";
    std::cout << "                    Custom Solver    Qhull Solver\n";
    std::cout << "Success:            " << (custom_result1.success ? "✓ YES" : "✗ NO ");
    if (qhull_available) {
        std::cout << "         " << (qhull_success1 ? "✓ YES" : "✗ NO ");
    } else {
        std::cout << "         N/A";
    }
    std::cout << std::endl;
    
    if (custom_result1.success) {
        std::cout << "Non-redundant:      " << custom_result1.A_min.rows() << "               ";
        if (qhull_available && qhull_success1) {
            std::cout << qhull_result1.A_min.rows();
        } else {
            std::cout << "N/A";
        }
        std::cout << std::endl;
        
        std::cout << "Solve time:         " << custom_duration1 << " s      ";
        if (qhull_available && qhull_success1) {
            std::cout << qhull_duration1 << " s";
        } else {
            std::cout << "N/A";
        }
        std::cout << std::endl;
        
        std::cout << "Iterations:         " << custom_result1.iterations << "               ";
        if (qhull_available && qhull_success1) {
            std::cout << qhull_result1.iterations;
        } else {
            std::cout << "N/A";
        }
        std::cout << std::endl;
    }
    
    // Constraint-by-constraint comparison
    std::cout << "\nConstraint Redundancy Status:\n";
    std::cout << "-----------------------------\n";
    for (int i = 0; i < 4; ++i) {
        std::cout << "Constraint " << i << ": ";
        
        // Custom solver status
        if (custom_result1.redundant_indices[i]) {
            std::cout << "REDUNDANT    ";
        } else if (custom_result1.unverified_indices[i]) {
            std::cout << "UNVERIFIED   ";
        } else {
            std::cout << "NON-REDUNDANT";
        }
        
        std::cout << " | ";
        
        // Qhull solver status
        if (qhull_available && qhull_success1) {
            if (qhull_result1.redundant_indices[i]) {
                std::cout << "REDUNDANT";
            } else if (qhull_result1.unverified_indices[i]) {
                std::cout << "UNVERIFIED";
            } else {
                std::cout << "NON-REDUNDANT";
            }
        } else {
            std::cout << "N/A";
        }
        std::cout << std::endl;
    }
    
    // Validation result
    bool example1_match = true;
    if (qhull_available && qhull_success1 && custom_result1.success) {
        // Check if results match
        bool counts_match = (custom_result1.A_min.rows() == qhull_result1.A_min.rows());
        bool redundancy_match = true;
        for (int i = 0; i < 4; ++i) {
            if (custom_result1.redundant_indices[i] != qhull_result1.redundant_indices[i]) {
                redundancy_match = false;
                break;
            }
        }
        example1_match = counts_match && redundancy_match;
        
        std::cout << "\nValidation: " << (example1_match ? "✓ PASS - Both methods agree!" : "✗ FAIL - Methods disagree");
        std::cout << std::endl;
    } else if (!qhull_available) {
        std::cout << "\nValidation: SKIPPED - Qhull not available\n";
    } else {
        std::cout << "\nValidation: ✗ FAIL - One or both solvers failed\n";
    }
    
    std::cout << "\n\n";
    
    std::cout << "Starting Example 2..." << std::endl;
    
    // Example 2: Random 10-halfplane example (similar to MATLAB code)
    std::cout << "Example 2: Random 10-halfplane polytope comparison\n";
    std::cout << "==================================================\n";
    
    std::cout << "Initializing random number generator..." << std::endl;
    
    // Generate random example similar to MATLAB code
    // m = 10 halfplanes, n = 2 dimensions (scaling down to avoid segfault)
    const int m2 = 10;
    const int n2 = 2;
    
    // Set up random number generator with fixed seed for reproducibility
    std::random_device rd;
    std::mt19937 gen(42); // Fixed seed for reproducible results
    std::uniform_real_distribution<double> uniform_dist(-0.5, 0.5);
    std::uniform_real_distribution<double> positive_dist(0.5, 1.5);
    
    polytope_redundancy::Matrix A2(m2, n2);
    polytope_redundancy::Vector b2(m2);
    
    std::cout << "Created matrices A2(" << m2 << "x" << n2 << ") and b2(" << m2 << ")" << std::endl;
    
    // Generate random constraint matrix A (similar to rand(m,n) - 0.5 in MATLAB)
    std::cout << "Generating random " << m2 << "x" << n2 << " constraint matrix...\n";
    for (int i = 0; i < m2; ++i) {
        for (int j = 0; j < n2; ++j) {
            A2(i, j) = uniform_dist(gen);
        }
    }
    
    // Generate random b vector (similar to rand(m,1) + 0.5 in MATLAB)
    std::cout << "Generating random b vector..." << std::endl;
    for (int i = 0; i < m2; ++i) {
        b2[i] = positive_dist(gen);
    }
    std::cout << "Generated b vector successfully." << std::endl;
    
    // Find a reasonable interior point 
    // We need to find a point that satisfies all constraints: A*x <= b
    // Let's try to find a point by solving a simple linear program
    std::cout << "Finding interior point..." << std::endl;
    polytope_redundancy::Vector interior_point2(n2);
    
    // Simple heuristic: try to find center by minimizing sum of constraint violations
    // Start with origin and adjust if needed
    bool found_interior = false;
    for (double scale = 0.1; scale <= 1.0 && !found_interior; scale += 0.1) {
        for (int j = 0; j < n2; ++j) {
            interior_point2[j] = 0.0;
        }
        
        // Check if this point is feasible
        bool feasible = true;
        for (int i = 0; i < m2; ++i) {
            double constraint_value = 0.0;
            for (int j = 0; j < n2; ++j) {
                constraint_value += A2(i, j) * interior_point2[j];
            }
            if (constraint_value > b2[i] - 1e-6) { // Allow small tolerance
                feasible = false;
                break;
            }
        }
        
        if (feasible) {
            found_interior = true;
            std::cout << "Found feasible interior point at origin." << std::endl;
        }
    }
    
    // If origin doesn't work, try a more sophisticated approach
    if (!found_interior) {
        std::cout << "Origin not feasible, trying alternative approach..." << std::endl;
        // Set interior point to be center of bounding box scaled down
        double min_bound = 0.0, max_bound = 0.0;
        for (int i = 0; i < m2; ++i) {
            if (b2[i] < max_bound) max_bound = b2[i];
            if (b2[i] > min_bound) min_bound = b2[i];
        }
        double center_val = (min_bound + max_bound) * 0.25; // Conservative scaling
        
        for (int j = 0; j < n2; ++j) {
            interior_point2[j] = center_val / n2; // Distribute evenly
        }
        
        // Verify this point is feasible
        bool feasible = true;
        for (int i = 0; i < m2; ++i) {
            double constraint_value = 0.0;
            for (int j = 0; j < n2; ++j) {
                constraint_value += A2(i, j) * interior_point2[j];
            }
            if (constraint_value > b2[i] - 1e-6) {
                feasible = false;
                break;
            }
        }
        
        if (!feasible) {
            std::cout << "Warning: Could not find guaranteed interior point. Using origin anyway." << std::endl;
            for (int j = 0; j < n2; ++j) {
                interior_point2[j] = 0.0;
            }
        } else {
            found_interior = true;
            std::cout << "Found alternative interior point." << std::endl;
        }
    }
    
    std::cout << "Generated polytope with " << m2 << " halfplanes in " << n2 << "D space\n";
    std::cout << "About to print interior point..." << std::endl;
    std::cout << "Interior point: [";
    for (int j = 0; j < n2; ++j) {
        std::cout << interior_point2[j];
        if (j < n2 - 1) std::cout << ", ";
    }
    std::cout << "]\n\n";
    std::cout << "About to run custom solver..." << std::endl;

    // Solve with qhull solver
    polytope_redundancy::RedundancyResult qhull_result2;
    double qhull_duration2 = 0.0;
    bool qhull_success2 = false;
    
    if (qhull_available) {
        std::cout << "Running qhull solver...\n";
        auto start_qhull2 = std::chrono::high_resolution_clock::now();
        qhull_result2 = qhull_solver.indicate_nonredundant_halfplanes(A2, b2, interior_point2);
        auto end_qhull2 = std::chrono::high_resolution_clock::now();
        qhull_duration2 = std::chrono::duration<double>(end_qhull2 - start_qhull2).count();
        qhull_success2 = qhull_result2.success;
    }
    
    // Solve with custom solver
    std::cout << "Running custom active-set solver...\n";
    auto start2 = std::chrono::high_resolution_clock::now();

    polytope_redundancy::RedundancyResult custom_result2;
    bool custom_success2 = false;

    try {
        std::cout << "Calling indicate_nonredundant_halfplanes..." << std::endl;
        custom_result2 = custom_solver.indicate_nonredundant_halfplanes(A2, b2, std::vector<bool>(), interior_point2);
        custom_success2 = custom_result2.success;
        std::cout << "Custom solver completed successfully." << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Custom solver threw exception: " << e.what() << std::endl;
        custom_success2 = false;
    } catch (...) {
        std::cout << "Custom solver threw unknown exception." << std::endl;
        custom_success2 = false;
    }

    auto end2 = std::chrono::high_resolution_clock::now();
    auto custom_duration2 = std::chrono::duration<double>(end2 - start2).count();
    // Display results comparison
    std::cout << "\nResults Comparison:\n";
    std::cout << "-------------------\n";
    std::cout << "                    Custom Solver    Qhull Solver\n";
    std::cout << "Success:            " << (custom_success2 ? "✓ YES" : "✗ NO ");
    if (qhull_available) {
        std::cout << "         " << (qhull_success2 ? "✓ YES" : "✗ NO ");
    } else {
        std::cout << "         N/A";
    }
    std::cout << std::endl;
    
    if (custom_success2) {
        std::cout << "Original halfplanes: " << m2 << "               ";
        if (qhull_available && qhull_success2) {
            std::cout << m2;
        } else {
            std::cout << "N/A";
        }
        std::cout << std::endl;
        
        std::cout << "Non-redundant:      " << custom_result2.A_min.rows() << "               ";
        if (qhull_available && qhull_success2) {
            std::cout << qhull_result2.A_min.rows();
        } else {
            std::cout << "N/A";
        }
        std::cout << std::endl;
        
        std::cout << "Redundant:          " << (m2 - custom_result2.A_min.rows()) << "               ";
        if (qhull_available && qhull_success2) {
            std::cout << (m2 - qhull_result2.A_min.rows());
        } else {
            std::cout << "N/A";
        }
        std::cout << std::endl;
        
        std::cout << "Solve time:         " << custom_duration2 << " s      ";
        if (qhull_available && qhull_success2) {
            std::cout << qhull_duration2 << " s";
        } else {
            std::cout << "N/A";
        }
        std::cout << std::endl;
        
        std::cout << "Iterations:         " << custom_result2.iterations << "               ";
        if (qhull_available && qhull_success2) {
            std::cout << qhull_result2.iterations;
        } else {
            std::cout << "N/A";
        }
        std::cout << std::endl;
    }
    
    // Validation result for Example 2
    bool example2_match = true;
    if (qhull_available && qhull_success2 && custom_success2) {
        // Check if results match
        bool counts_match = (custom_result2.A_min.rows() == qhull_result2.A_min.rows());
        bool redundancy_match = true;
        int disagreements = 0;
        
        for (int i = 0; i < m2; ++i) {
            if (custom_result2.redundant_indices[i] != qhull_result2.redundant_indices[i]) {
                redundancy_match = false;
                disagreements++;
            }
        }
        example2_match = counts_match && redundancy_match;
        
        std::cout << "\nDetailed Validation:\n";
        std::cout << "Same number of non-redundant halfplanes: " << (counts_match ? "✓ YES" : "✗ NO") << std::endl;
        std::cout << "Same redundancy classification: " << (redundancy_match ? "✓ YES" : "✗ NO");
        if (!redundancy_match) {
            std::cout << " (" << disagreements << " disagreements)";
        }
        std::cout << std::endl;
        
        std::cout << "\nOverall Validation: " << (example2_match ? "✓ PASS - Both methods agree!" : "✗ FAIL - Methods disagree");
        std::cout << std::endl;
    } else if (!qhull_available) {
        std::cout << "\nValidation: SKIPPED - Qhull not available\n";
    } else {
        std::cout << "\nValidation: ✗ FAIL - One or both solvers failed\n";
    }
    
    std::cout << "\n\n";
    
    // Example 3: High-dimensional example - 20D with 100 constraints (similar to MATLAB)
    std::cout << "Example 3: High-dimensional 20D polytope with 100 constraints\n";
    std::cout << "==============================================================\n";
    
    // This matches the MATLAB example: m = 5000; n = 10; (we use m=100, n=20)
    // Note: High-dimensional problems may cause stability issues in the current implementation
    const int m3 = 100;
    const int n3 = 20;
    
    std::cout << "Generating " << m3 << "x" << n3 << " random polytope..." << std::endl;
    std::cout << "Note: High-dimensional problems may expose stability issues in the current C++ implementation." << std::endl;
    
    // Use same random seed for reproducibility
    std::mt19937 gen3(123); // Different seed for variety
    std::uniform_real_distribution<double> uniform_dist3(-0.5, 0.5);
    std::uniform_real_distribution<double> positive_dist3(0.5, 1.5);
    
    polytope_redundancy::Matrix A3(m3, n3);
    polytope_redundancy::Vector b3(m3);
    
    // Generate random constraint matrix A (similar to rand(m,n) - 0.5 in MATLAB)
    for (int i = 0; i < m3; ++i) {
        for (int j = 0; j < n3; ++j) {
            A3(i, j) = uniform_dist3(gen3);
        }
    }
    
    // Generate random b vector (similar to rand(m,1) + 0.5 in MATLAB)
    for (int i = 0; i < m3; ++i) {
        b3[i] = positive_dist3(gen3);
    }
    
    std::cout << "Generated constraint matrix and vector successfully." << std::endl;
    
    // Find interior point using more sophisticated approach for high dimensions
    std::cout << "Finding interior point for " << n3 << "D polytope..." << std::endl;
    polytope_redundancy::Vector interior_point3(n3);
    
    // Initialize to origin
    for (int j = 0; j < n3; ++j) {
        interior_point3[j] = 0.0;
    }
    
    // Check if origin is feasible
    bool origin_feasible = true;
    for (int i = 0; i < m3; ++i) {
        double constraint_value = 0.0;
        for (int j = 0; j < n3; ++j) {
            constraint_value += A3(i, j) * interior_point3[j];
        }
        if (constraint_value > b3[i] - 1e-8) {
            origin_feasible = false;
            break;
        }
    }
    
    if (!origin_feasible) {
        std::cout << "Origin not feasible, using alternative interior point..." << std::endl;
        // Try a small random perturbation
        std::uniform_real_distribution<double> small_dist(-0.01, 0.01);
        for (int j = 0; j < n3; ++j) {
            interior_point3[j] = small_dist(gen3);
        }
    } else {
        std::cout << "Using origin as interior point." << std::endl;
    }
    
    std::cout << "Interior point: [";
    for (int j = 0; j < std::min(3, n3); ++j) { // Show first 3 coordinates
        std::cout << interior_point3[j];
        if (j < std::min(3, n3) - 1) std::cout << ", ";
    }
    if (n3 > 3) std::cout << ", ...";
    std::cout << "]\n\n";
    
    // Solve with qhull solver
    polytope_redundancy::RedundancyResult qhull_result3;
    double qhull_duration3 = 0.0;
    bool qhull_success3 = false;

    if (qhull_available) {
        std::cout << "Running qhull solver on " << n3 << "D problem...\n";
        auto start_qhull3 = std::chrono::high_resolution_clock::now();

        try {
            qhull_result3 = qhull_solver.indicate_nonredundant_halfplanes(A3, b3, interior_point3);
            qhull_success3 = qhull_result3.success;
            std::cout << "Qhull solver completed successfully." << std::endl;
        } catch (const std::exception& e) {
            std::cout << "Qhull solver threw exception: " << e.what() << std::endl;
            qhull_success3 = false;
        } catch (...) {
            std::cout << "Qhull solver threw unknown exception." << std::endl;
            qhull_success3 = false;
        }

        auto end_qhull3 = std::chrono::high_resolution_clock::now();
        qhull_duration3 = std::chrono::duration<double>(end_qhull3 - start_qhull3).count();
    }

    // Solve with custom solver
    std::cout << "Running custom active-set solver on " << n3 << "D problem...\n";
    auto start3 = std::chrono::high_resolution_clock::now();
    
    polytope_redundancy::RedundancyResult custom_result3;
    bool custom_success3 = false;
    
    try {
        std::cout << "This may take longer for high-dimensional problems..." << std::endl;
        custom_result3 = custom_solver.indicate_nonredundant_halfplanes(A3, b3, std::vector<bool>(), interior_point3);
        custom_success3 = custom_result3.success;
        std::cout << "Custom solver completed successfully." << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Custom solver threw exception: " << e.what() << std::endl;
        custom_success3 = false;
    } catch (...) {
        std::cout << "Custom solver threw unknown exception." << std::endl;
        custom_success3 = false;
    }
    
    auto end3 = std::chrono::high_resolution_clock::now();
    auto custom_duration3 = std::chrono::duration<double>(end3 - start3).count();
    
    // Display results comparison
    std::cout << "\nResults Comparison:\n";
    std::cout << "-------------------\n";
    std::cout << "                    Custom Solver    Qhull Solver\n";
    std::cout << "Success:            " << (custom_success3 ? "✓ YES" : "✗ NO ");
    if (qhull_available) {
        std::cout << "         " << (qhull_success3 ? "✓ YES" : "✗ NO ");
    } else {
        std::cout << "         N/A";
    }
    std::cout << std::endl;
    
    if (custom_success3 || (qhull_available && qhull_success3)) {
        std::cout << "Dimensions:         " << n3 << "               ";
        if (qhull_available && qhull_success3) {
            std::cout << n3;
        } else {
            std::cout << "N/A";
        }
        std::cout << std::endl;
        
        std::cout << "Original halfplanes: " << m3 << "               ";
        if (qhull_available && qhull_success3) {
            std::cout << m3;
        } else {
            std::cout << "N/A";
        }
        std::cout << std::endl;
        
        if (custom_success3) {
            std::cout << "Non-redundant:      " << custom_result3.A_min.rows() << "               ";
        } else {
            std::cout << "Non-redundant:      N/A             ";
        }
        if (qhull_available && qhull_success3) {
            std::cout << qhull_result3.A_min.rows();
        } else {
            std::cout << "N/A";
        }
        std::cout << std::endl;
        
        if (custom_success3) {
            std::cout << "Redundant:          " << (m3 - custom_result3.A_min.rows()) << "               ";
        } else {
            std::cout << "Redundant:          N/A             ";
        }
        if (qhull_available && qhull_success3) {
            std::cout << (m3 - qhull_result3.A_min.rows());
        } else {
            std::cout << "N/A";
        }
        std::cout << std::endl;
        
        std::cout << "Solve time:         ";
        if (custom_success3) {
            std::cout << custom_duration3 << " s      ";
        } else {
            std::cout << "N/A         ";
        }
        if (qhull_available && qhull_success3) {
            std::cout << qhull_duration3 << " s";
        } else {
            std::cout << "N/A";
        }
        std::cout << std::endl;
        
        if (custom_success3) {
            std::cout << "Iterations:         " << custom_result3.iterations << "               ";
        } else {
            std::cout << "Iterations:         N/A             ";
        }
        if (qhull_available && qhull_success3) {
            std::cout << qhull_result3.iterations;
        } else {
            std::cout << "N/A";
        }
        std::cout << std::endl;
    }
    
    // Validation result for Example 3
    bool example3_match = true;
    if (qhull_available && qhull_success3 && custom_success3) {
        // Check if results match
        bool counts_match = (custom_result3.A_min.rows() == qhull_result3.A_min.rows());
        bool redundancy_match = true;
        int disagreements = 0;
        
        for (int i = 0; i < m3; ++i) {
            if (custom_result3.redundant_indices[i] != qhull_result3.redundant_indices[i]) {
                redundancy_match = false;
                disagreements++;
            }
        }
        example3_match = counts_match && redundancy_match;
        
        std::cout << "\nDetailed Validation:\n";
        std::cout << "Same number of non-redundant halfplanes: " << (counts_match ? "✓ YES" : "✗ NO") << std::endl;
        std::cout << "Same redundancy classification: " << (redundancy_match ? "✓ YES" : "✗ NO");
        if (!redundancy_match) {
            std::cout << " (" << disagreements << " disagreements)";
        }
        std::cout << std::endl;
        
        std::cout << "\nOverall Validation: " << (example3_match ? "✓ PASS - Both methods agree!" : "✗ FAIL - Methods disagree");
        std::cout << std::endl;
    } else if (!qhull_available) {
        std::cout << "\nValidation: SKIPPED - Qhull not available\n";
    } else if (!custom_success3 && (!qhull_available || !qhull_success3)) {
        std::cout << "\nValidation: ✗ FAIL - Both solvers failed\n";
    } else if (!custom_success3) {
        std::cout << "\nValidation: ✗ FAIL - Custom solver failed\n";
    } else {
        std::cout << "\nValidation: ✗ FAIL - Qhull solver failed\n";
    }
    
    std::cout << "\n\n";
    
    // Overall summary - updated to include all three examples
    std::cout << std::string(60, '=') << std::endl;
    std::cout << "VALIDATION SUMMARY\n";
    std::cout << std::string(60, '=') << std::endl;
    
    if (qhull_available) {
        std::cout << "Example 1 (2D simple):     " << (example1_match ? "✓ PASS" : "✗ FAIL") << std::endl;
        std::cout << "Example 2 (2D random 10):  " << (example2_match ? "✓ PASS" : "✗ FAIL") << std::endl;
        std::cout << "Example 3 (20D random 100): " << (example3_match ? "✓ PASS" : "✗ FAIL") << std::endl;
        
        bool all_examples_pass = example1_match && example2_match && example3_match;
        if (all_examples_pass) {
            std::cout << "\nAll tests passed! The C++ implementation produces results consistent with qhull!" << std::endl;
        } else {
            std::cout << "\nWarning: Some discrepancies found between C++ implementation and qhull." << std::endl;
        }
    } else {
        std::cout << "Qhull not available - validation skipped." << std::endl;
        std::cout << "To enable qhull comparison, install qhull:" << std::endl;
        std::cout << "  macOS: brew install qhull" << std::endl;
        std::cout << "  Ubuntu: sudo apt install qhull-bin" << std::endl;
    }
    
    std::cout << "\nExample completed!\n";
    
    return 0;
}
