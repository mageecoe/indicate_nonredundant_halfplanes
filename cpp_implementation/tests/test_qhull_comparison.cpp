#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "polytope_redundancy/core.hpp"
#include "polytope_redundancy/qhull_solver.hpp"
#include <chrono>
#include <random>

using namespace polytope_redundancy;
using Catch::Approx;

TEST_CASE("QhullSolver basic functionality", "[qhull]") {
    QhullSolver qhull_solver;
    
    // Skip test if qhalf is not available
    if (!qhull_solver.is_qhalf_available()) {
        SKIP("qhalf not available, skipping qhull tests");
    }
    
    SECTION("Simple 2D example") {
        // Create simple 2D polytope with one redundant constraint
        Matrix A(4, 2);
        A(0, 0) = -1.0; A(0, 1) = 1.0;   // -x + y <= 1
        A(1, 0) = 1.0;  A(1, 1) = 0.0;   //  x     <= 2  (redundant)
        A(2, 0) = 0.0;  A(2, 1) = -1.0;  //     -y <= 0.3
        A(3, 0) = 1.0;  A(3, 1) = 0.0;   //  x     <= 1
        
        Vector b(4);
        b[0] = 1.0; b[1] = 2.0; b[2] = 0.3; b[3] = 1.0;
        
        Vector interior_point(2);
        interior_point[0] = 0.0; interior_point[1] = 0.0;
        
        auto result = qhull_solver.indicate_nonredundant_halfplanes(A, b, interior_point);
        
        REQUIRE(result.success);
        REQUIRE(result.A_min.rows() == 3);  // Should remove one redundant constraint
        REQUIRE(result.redundant_indices[1] == true);  // Constraint 1 should be redundant
    }
}

TEST_CASE("Compare QhullSolver vs PolytopeRedundancyRemover", "[qhull][comparison]") {
    QhullSolver qhull_solver;
    PolytopeRedundancyRemover custom_solver;
    
    // Skip test if qhalf is not available
    if (!qhull_solver.is_qhalf_available()) {
        SKIP("qhalf not available, skipping comparison tests");
    }
    
    SECTION("2D polytope comparison") {
        Matrix A(4, 2);
        A(0, 0) = -1.0; A(0, 1) = 1.0;   // -x + y <= 1
        A(1, 0) = 1.0;  A(1, 1) = 0.0;   //  x     <= 2  (redundant)
        A(2, 0) = 0.0;  A(2, 1) = -1.0;  //     -y <= 0.3
        A(3, 0) = 1.0;  A(3, 1) = 0.0;   //  x     <= 1
        
        Vector b(4);
        b[0] = 1.0; b[1] = 2.0; b[2] = 0.3; b[3] = 1.0;
        
        Vector interior_point(2);
        interior_point[0] = 0.0; interior_point[1] = 0.0;
        
        auto qhull_result = qhull_solver.indicate_nonredundant_halfplanes(A, b, interior_point);
        auto custom_result = custom_solver.indicate_nonredundant_halfplanes(A, b, std::vector<bool>(), interior_point);
        
        REQUIRE(qhull_result.success);
        REQUIRE(custom_result.success);
        
        // Both should find the same number of non-redundant constraints
        REQUIRE(qhull_result.A_min.rows() == custom_result.A_min.rows());
        
        // Compare redundant indices
        for (int i = 0; i < A.rows(); ++i) {
            INFO("Checking constraint " << i);
            REQUIRE(qhull_result.redundant_indices[i] == custom_result.redundant_indices[i]);
        }
    }
    
    SECTION("3D polytope comparison") {
        // Create a 3D cube with extra redundant constraints
        Matrix A(8, 3);
        Vector b(8);
        
        // Standard cube constraints: -1 <= x,y,z <= 1
        A(0, 0) = 1.0;  A(0, 1) = 0.0;  A(0, 2) = 0.0;  b[0] = 1.0;   // x <= 1
        A(1, 0) = -1.0; A(1, 1) = 0.0;  A(1, 2) = 0.0;  b[1] = 1.0;   // -x <= 1
        A(2, 0) = 0.0;  A(2, 1) = 1.0;  A(2, 2) = 0.0;  b[2] = 1.0;   // y <= 1
        A(3, 0) = 0.0;  A(3, 1) = -1.0; A(3, 2) = 0.0;  b[3] = 1.0;   // -y <= 1
        A(4, 0) = 0.0;  A(4, 1) = 0.0;  A(4, 2) = 1.0;  b[4] = 1.0;   // z <= 1
        A(5, 0) = 0.0;  A(5, 1) = 0.0;  A(5, 2) = -1.0; b[5] = 1.0;   // -z <= 1
        
        // Add redundant constraints
        A(6, 0) = 1.0;  A(6, 1) = 0.0;  A(6, 2) = 0.0;  b[6] = 2.0;   // x <= 2 (redundant)
        A(7, 0) = 0.0;  A(7, 1) = 1.0;  A(7, 2) = 0.0;  b[7] = 1.5;   // y <= 1.5 (redundant)
        
        Vector interior_point(3, 0.0);  // Origin
        
        auto qhull_result = qhull_solver.indicate_nonredundant_halfplanes(A, b, interior_point);
        auto custom_result = custom_solver.indicate_nonredundant_halfplanes(A, b, std::vector<bool>(), interior_point);
        
        REQUIRE(qhull_result.success);
        REQUIRE(custom_result.success);
        
        // Both should identify the same redundant constraints
        REQUIRE(qhull_result.A_min.rows() == custom_result.A_min.rows());
        REQUIRE(qhull_result.A_min.rows() == 6);  // Should keep only the 6 cube faces
        
        // Constraints 6 and 7 should be redundant
        REQUIRE(qhull_result.redundant_indices[6] == true);
        REQUIRE(qhull_result.redundant_indices[7] == true);
        REQUIRE(custom_result.redundant_indices[6] == true);
        REQUIRE(custom_result.redundant_indices[7] == true);
    }
    
    SECTION("Random polytope comparison") {
        std::mt19937 gen(42);  // Fixed seed for reproducibility
        std::uniform_real_distribution<double> dis(-1.0, 1.0);
        
        int m = 20;  // number of constraints
        int n = 5;   // dimension
        
        Matrix A(m, n);
        Vector b(m);
        
        // Generate random constraints
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                A(i, j) = dis(gen);
            }
            b[i] = 1.0 + std::abs(dis(gen));  // Ensure positive b values
        }
        
        Vector interior_point(n, 0.0);  // Try origin
        
        // Check if origin is interior
        bool origin_interior = true;
        for (int i = 0; i < m; ++i) {
            double dot_product = 0.0;
            for (int j = 0; j < n; ++j) {
                dot_product += A(i, j) * interior_point[j];
            }
            if (dot_product >= b[i] - 1e-6) {
                origin_interior = false;
                break;
            }
        }
        
        if (!origin_interior) {
            // Try a different interior point
            for (int j = 0; j < n; ++j) {
                interior_point[j] = 0.1;
            }
        }
        
        auto qhull_result = qhull_solver.indicate_nonredundant_halfplanes(A, b, interior_point);
        auto custom_result = custom_solver.indicate_nonredundant_halfplanes(A, b, std::vector<bool>(), interior_point);
        
        if (qhull_result.success && custom_result.success) {
            // If both succeed, they should agree on the number of non-redundant constraints
            INFO("Qhull found " << qhull_result.A_min.rows() << " non-redundant constraints");
            INFO("Custom found " << custom_result.A_min.rows() << " non-redundant constraints");
            
            // Allow for small differences in numerical tolerances
            REQUIRE(std::abs(static_cast<int>(qhull_result.A_min.rows()) - 
                           static_cast<int>(custom_result.A_min.rows())) <= 2);
        }
    }
}

TEST_CASE("Performance comparison", "[qhull][performance][.]") {  // '.' tag excludes from default runs
    QhullSolver qhull_solver;
    PolytopeRedundancyRemover custom_solver;
    
    if (!qhull_solver.is_qhalf_available()) {
        SKIP("qhalf not available, skipping performance tests");
    }
    
    SECTION("Medium-sized problem performance") {
        std::mt19937 gen(123);
        std::uniform_real_distribution<double> dis(-1.0, 1.0);
        
        int m = 100;  // constraints
        int n = 10;   // dimension
        
        Matrix A(m, n);
        Vector b(m);
        
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                A(i, j) = dis(gen);
            }
            b[i] = 2.0;  // Large enough to likely include origin
        }
        
        Vector interior_point(n, 0.0);
        
        // Time qhull solver
        auto start_qhull = std::chrono::high_resolution_clock::now();
        auto qhull_result = qhull_solver.indicate_nonredundant_halfplanes(A, b, interior_point);
        auto end_qhull = std::chrono::high_resolution_clock::now();
        auto qhull_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_qhull - start_qhull);
        
        // Time custom solver
        auto start_custom = std::chrono::high_resolution_clock::now();
        auto custom_result = custom_solver.indicate_nonredundant_halfplanes(A, b, std::vector<bool>(), interior_point);
        auto end_custom = std::chrono::high_resolution_clock::now();
        auto custom_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_custom - start_custom);
        
        std::cout << "Performance comparison (m=" << m << ", n=" << n << "):\n";
        std::cout << "  Qhull solver: " << qhull_time.count() << " ms\n";
        std::cout << "  Custom solver: " << custom_time.count() << " ms\n";
        
        if (qhull_result.success && custom_result.success) {
            std::cout << "  Qhull found: " << qhull_result.A_min.rows() << " non-redundant constraints\n";
            std::cout << "  Custom found: " << custom_result.A_min.rows() << " non-redundant constraints\n";
        }
        
        // Both should succeed for this test case
        REQUIRE(qhull_result.success);
        REQUIRE(custom_result.success);
    }
}

TEST_CASE("QhullSolver edge cases", "[qhull][edge]") {
    QhullSolver qhull_solver;
    
    if (!qhull_solver.is_qhalf_available()) {
        SKIP("qhalf not available, skipping edge case tests");
    }
    
    SECTION("Empty constraint matrix") {
        Matrix A(0, 2);
        Vector b(0);
        Vector interior_point(2, 0.0);
        
        auto result = qhull_solver.indicate_nonredundant_halfplanes(A, b, interior_point);
        
        REQUIRE(result.success);
        REQUIRE(result.A_min.rows() == 0);
        REQUIRE(result.redundant_indices.empty());
    }
    
    SECTION("Single constraint") {
        Matrix A(1, 2);
        A(0, 0) = 1.0; A(0, 1) = 0.0;  // x <= 1
        
        Vector b(1);
        b[0] = 1.0;
        
        Vector interior_point(2, 0.0);  // Origin should be interior
        
        auto result = qhull_solver.indicate_nonredundant_halfplanes(A, b, interior_point);
        
        REQUIRE(result.success);
        REQUIRE(result.A_min.rows() == 1);
        REQUIRE(result.redundant_indices[0] == false);
    }
}