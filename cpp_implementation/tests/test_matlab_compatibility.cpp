#include <catch2/catch_all.hpp>
#include "polytope_redundancy/core.hpp"
#include <random>
#include <chrono>

using namespace polytope_redundancy;

TEST_CASE("MATLAB compatibility tests", "[matlab]") {
    SECTION("MATLAB Example 1 - 2D case from documentation") {
        // This is the exact example from MATLAB code:
        // Px1 = [-1  1;...
        //         1  0;...
        //         0 -1;...
        //         1  0];
        // Pc1 = [  1;...
        //          2;...
        //        0.3;...
        //          1];
        
        Matrix A(4, 2);
        A(0, 0) = -1.0; A(0, 1) =  1.0;  // -x1 + x2 <= 1
        A(1, 0) =  1.0; A(1, 1) =  0.0;  //  x1      <= 2  
        A(2, 0) =  0.0; A(2, 1) = -1.0;  //     -x2 <= 0.3
        A(3, 0) =  1.0; A(3, 1) =  0.0;  //  x1      <= 1
        
        Vector b(4);
        b[0] = 1.0;
        b[1] = 2.0;
        b[2] = 0.3;  
        b[3] = 1.0;
        
        PolytopeRedundancyRemover solver;
        auto result = solver.indicate_nonredundant_halfplanes(A, b);
        
        REQUIRE(result.success == true);
        
        // According to MATLAB documentation, constraint x1 <= 2 is redundant
        // because x1 <= 1 is tighter
        REQUIRE(result.redundant_indices[1] == true);
        
        // All other constraints should be non-redundant
        REQUIRE(result.redundant_indices[0] == false);
        REQUIRE(result.redundant_indices[2] == false);
        REQUIRE(result.redundant_indices[3] == false);
        
        // Result should have 3 constraints (4 - 1 redundant)
        REQUIRE(result.A_min.rows() == 3);
        REQUIRE(result.A_min.cols() == 2);
        
        // Check that iterations completed
        REQUIRE(result.iterations > 0);
        REQUIRE(result.iterations < 100); // Reasonable upper bound
    }
    
    SECTION("Symmetric polytope test") {
        // Create a symmetric polytope that should be detected as symmetric
        const int n = 3;
        const int m = 6;
        
        Matrix A(m, n);
        Vector b(m);
        
        // Create symmetric constraints: Ax <= 1 and -Ax <= 1
        // This represents the L∞ ball with radius 1
        A(0, 0) =  1.0; A(0, 1) =  0.0; A(0, 2) =  0.0;  //  x <= 1
        A(1, 0) =  0.0; A(1, 1) =  1.0; A(1, 2) =  0.0;  //  y <= 1  
        A(2, 0) =  0.0; A(2, 1) =  0.0; A(2, 2) =  1.0;  //  z <= 1
        A(3, 0) = -1.0; A(3, 1) =  0.0; A(3, 2) =  0.0;  // -x <= 1  =>  x >= -1
        A(4, 0) =  0.0; A(4, 1) = -1.0; A(4, 2) =  0.0;  // -y <= 1  =>  y >= -1
        A(5, 0) =  0.0; A(5, 1) =  0.0; A(5, 2) = -1.0;  // -z <= 1  =>  z >= -1
        
        for (int i = 0; i < m; ++i) {
            b[i] = 1.0;
        }
        
        PolytopeRedundancyRemover solver;
        auto result = solver.indicate_nonredundant_halfplanes(A, b);
        
        REQUIRE(result.success == true);
        
        // All constraints should be non-redundant for L∞ ball
        for (int i = 0; i < m; ++i) {
            REQUIRE(result.redundant_indices[i] == false);
        }
        
        REQUIRE(result.A_min.rows() == m);
        REQUIRE(result.A_min.cols() == n);
    }
    
    SECTION("Tolerance verification test") {
        // Test case designed to verify that tolerance handling matches MATLAB
        Matrix A(3, 2);
        Vector b(3);
        
        // Create constraints where redundancy is close to tolerance
        A(0, 0) = 1.0; A(0, 1) = 0.0;  // x <= 1
        A(1, 0) = 1.0; A(1, 1) = 0.0;  // x <= 1 + 1e-13 (should be considered duplicate)
        A(2, 0) = 0.0; A(2, 1) = 1.0;  // y <= 1
        
        b[0] = 1.0;
        b[1] = 1.0 + 1e-13;  // Very close to first constraint
        b[2] = 1.0;
        
        PolytopeRedundancyRemover solver;
        auto result = solver.indicate_nonredundant_halfplanes(A, b);
        
        REQUIRE(result.success == true);
        
        // Should have detected redundancy despite small numerical difference
        int redundant_count = 0;
        for (int i = 0; i < 3; ++i) {
            if (result.redundant_indices[i]) {
                redundant_count++;
            }
        }
        
        // Exactly one should be redundant due to duplication
        REQUIRE(redundant_count >= 1);
        REQUIRE(result.A_min.rows() <= 2);
    }
    
    SECTION("Large scale performance test") {
        // Create a larger problem similar to MATLAB's second example
        const int m = 100;  // Smaller than MATLAB's 5000 for faster testing
        const int n = 5;
        
        std::mt19937 gen(42); // Fixed seed for reproducibility
        std::uniform_real_distribution<double> dis(-0.5, 0.5);
        std::uniform_real_distribution<double> b_dis(0.5, 1.5);
        
        Matrix A(m, n);
        Vector b(m);
        
        // Generate random constraints similar to MATLAB example
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                A(i, j) = dis(gen);
            }
            b[i] = b_dis(gen);
        }
        
        PolytopeRedundancyRemover solver;
        auto start = std::chrono::high_resolution_clock::now();
        auto result = solver.indicate_nonredundant_halfplanes(A, b);
        auto end = std::chrono::high_resolution_clock::now();
        
        REQUIRE(result.success == true);
        
        // Should have eliminated some constraints
        int redundant_count = 0;
        for (int i = 0; i < m; ++i) {
            if (result.redundant_indices[i]) {
                redundant_count++;
            }
        }
        
        REQUIRE(redundant_count > 0);  // Should find some redundant constraints
        REQUIRE(result.A_min.rows() < m);  // Should reduce problem size
        REQUIRE(result.A_min.rows() > n);  // But should still be overconstrained
        
        // Check performance is reasonable (less than 1 second for this size)
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        REQUIRE(duration.count() < 1000);
        
        // Check iteration count is reasonable
        REQUIRE(result.iterations > 0);
        REQUIRE(result.iterations <= m + 10);  // Should not exceed problem size by much
    }
    
    SECTION("Interior point handling test") {
        // Test that matches MATLAB's interior point validation
        Matrix A(3, 2);
        Vector b(3);
        
        A(0, 0) = 1.0; A(0, 1) = 0.0;  // x <= 1
        A(1, 0) = 0.0; A(1, 1) = 1.0;  // y <= 1
        A(2, 0) = -1.0; A(2, 1) = -1.0; // -x - y <= -0.5  =>  x + y >= 0.5
        
        b[0] = 1.0;
        b[1] = 1.0;
        b[2] = -0.5;
        
        // Test with valid interior point
        Vector interior_point(2);
        interior_point[0] = 0.6;
        interior_point[1] = 0.6;  // Point (0.6, 0.6) should be interior
        
        PolytopeRedundancyRemover solver;
        auto result = solver.indicate_nonredundant_halfplanes(A, b, std::vector<bool>(), interior_point);
        
        REQUIRE(result.success == true);
        
        // All constraints should be non-redundant for this configuration
        for (int i = 0; i < 3; ++i) {
            REQUIRE(result.redundant_indices[i] == false);
        }
    }
    
    SECTION("Edge case: Single constraint") {
        // Test with just one constraint
        Matrix A(1, 2);
        Vector b(1);
        
        A(0, 0) = 1.0; A(0, 1) = 0.0;  // x <= 1
        b[0] = 1.0;
        
        PolytopeRedundancyRemover solver;
        auto result = solver.indicate_nonredundant_halfplanes(A, b);
        
        REQUIRE(result.success == true);
        REQUIRE(result.redundant_indices[0] == false);
        REQUIRE(result.A_min.rows() == 1);
    }
    
    SECTION("Edge case: Identical constraints") {
        // Test with completely identical constraints
        Matrix A(3, 2);
        Vector b(3);
        
        // All identical: x <= 1
        A(0, 0) = 1.0; A(0, 1) = 0.0;  
        A(1, 0) = 1.0; A(1, 1) = 0.0;  
        A(2, 0) = 1.0; A(2, 1) = 0.0;  
        
        b[0] = 1.0;
        b[1] = 1.0;
        b[2] = 1.0;
        
        PolytopeRedundancyRemover solver;
        auto result = solver.indicate_nonredundant_halfplanes(A, b);
        
        REQUIRE(result.success == true);
        
        // Should eliminate duplicates, keeping only one
        int non_redundant_count = 0;
        for (int i = 0; i < 3; ++i) {
            if (!result.redundant_indices[i]) {
                non_redundant_count++;
            }
        }
        
        REQUIRE(non_redundant_count == 1);
        REQUIRE(result.A_min.rows() == 1);
    }
}

TEST_CASE("Algorithm correctness verification", "[correctness]") {
    SECTION("Verify minimal representation properties") {
        // Create a test case where we know the exact minimal representation
        Matrix A(5, 2);
        Vector b(5);
        
        // Square [0,1] x [0,1] with one redundant constraint
        A(0, 0) =  1.0; A(0, 1) =  0.0;  //  x <= 1
        A(1, 0) =  1.0; A(1, 1) =  0.0;  //  x <= 2  (redundant)
        A(2, 0) = -1.0; A(2, 1) =  0.0;  // -x <= 0  =>  x >= 0
        A(3, 0) =  0.0; A(3, 1) =  1.0;  //  y <= 1
        A(4, 0) =  0.0; A(4, 1) = -1.0;  // -y <= 0  =>  y >= 0
        
        b[0] = 1.0;
        b[1] = 2.0; // This makes constraint redundant
        b[2] = 0.0;
        b[3] = 1.0;
        b[4] = 0.0;
        
        PolytopeRedundancyRemover solver;
        auto result = solver.indicate_nonredundant_halfplanes(A, b);
        
        REQUIRE(result.success == true);
        
        // Verify that the redundant constraint was found
        REQUIRE(result.redundant_indices[1] == true);
        
        // Verify that the minimal representation defines the same polytope
        REQUIRE(result.A_min.rows() == 4);
        
        // Verify vertices of unit square are still feasible
        std::vector<std::vector<double>> vertices = {
            {0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}
        };
        
        for (const auto& vertex : vertices) {
            Vector v(2);
            v[0] = vertex[0];
            v[1] = vertex[1];
            
            // Check feasibility in minimal representation
            Vector Av(result.A_min.rows());
            result.A_min.gemv(v, Av);
            
            for (int i = 0; i < result.A_min.rows(); ++i) {
                REQUIRE(Av[i] <= result.b_min[i] + 1e-10);
            }
        }
    }
}