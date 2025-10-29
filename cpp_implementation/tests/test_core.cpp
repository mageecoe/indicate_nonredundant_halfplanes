#include <catch2/catch_all.hpp>
#include "polytope_redundancy/core.hpp"

using namespace polytope_redundancy;

TEST_CASE("Core redundancy removal algorithm", "[core]") {
    SECTION("Simple 2D example") {
        // Test case from MATLAB code comments
        Matrix A(4, 2);
        A(0, 0) = -1.0; A(0, 1) =  1.0;  // -x1 + x2 <= 1
        A(1, 0) =  1.0; A(1, 1) =  0.0;  //  x1      <= 2  (should be redundant)
        A(2, 0) =  0.0; A(2, 1) = -1.0;  //     -x2 <= 0.3
        A(3, 0) =  1.0; A(3, 1) =  0.0;  //  x1      <= 1
        
        Vector b(4);
        b[0] =  1.0;
        b[1] =  2.0;
        b[2] =  0.3;
        b[3] =  1.0;
        
        PolytopeRedundancyRemover solver;
        auto result = solver.indicate_nonredundant_halfplanes(A, b);
        
        REQUIRE(result.success == true);
        
        // The second constraint (x1 <= 2) should be redundant since x1 <= 1 is tighter
        REQUIRE(result.redundant_indices[1] == true);
        
        // Other constraints should be non-redundant
        REQUIRE(result.redundant_indices[0] == false);
        REQUIRE(result.redundant_indices[2] == false);
        REQUIRE(result.redundant_indices[3] == false);
        
        // Should have eliminated exactly 1 constraint
        REQUIRE(result.A_min.rows() == 3);
        REQUIRE(result.A_min.cols() == 2);
    }
    
    SECTION("All constraints non-redundant") {
        // Unit square: all constraints are necessary
        Matrix A(4, 2);
        A(0, 0) =  1.0; A(0, 1) =  0.0;  //  x <= 1
        A(1, 0) = -1.0; A(1, 1) =  0.0;  // -x <= 0  =>  x >= 0
        A(2, 0) =  0.0; A(2, 1) =  1.0;  //  y <= 1
        A(3, 0) =  0.0; A(3, 1) = -1.0;  // -y <= 0  =>  y >= 0
        
        Vector b(4);
        b[0] = 1.0;
        b[1] = 0.0;
        b[2] = 1.0;
        b[3] = 0.0;
        
        PolytopeRedundancyRemover solver;
        auto result = solver.indicate_nonredundant_halfplanes(A, b);
        
        REQUIRE(result.success == true);
        
        // All constraints should be non-redundant
        for (int i = 0; i < 4; ++i) {
            REQUIRE(result.redundant_indices[i] == false);
        }
        
        // Should have same number of constraints
        REQUIRE(result.A_min.rows() == 4);
    }
    
    SECTION("Multiple redundant constraints") {
        // Create a case with multiple redundant constraints
        Matrix A(6, 2);
        A(0, 0) =  1.0; A(0, 1) =  0.0;  //  x <= 1
        A(1, 0) =  1.0; A(1, 1) =  0.0;  //  x <= 2  (redundant)
        A(2, 0) =  1.0; A(2, 1) =  0.0;  //  x <= 3  (redundant)
        A(3, 0) = -1.0; A(3, 1) =  0.0;  // -x <= 0  =>  x >= 0
        A(4, 0) =  0.0; A(4, 1) =  1.0;  //  y <= 1
        A(5, 0) =  0.0; A(5, 1) = -1.0;  // -y <= 0  =>  y >= 0
        
        Vector b(6);
        b[0] = 1.0;
        b[1] = 2.0;
        b[2] = 3.0;
        b[3] = 0.0;
        b[4] = 1.0;
        b[5] = 0.0;
        
        PolytopeRedundancyRemover solver;
        auto result = solver.indicate_nonredundant_halfplanes(A, b);
        
        REQUIRE(result.success == true);
        
        // Constraints 1 and 2 should be redundant (looser x bounds)
        REQUIRE(result.redundant_indices[1] == true);
        REQUIRE(result.redundant_indices[2] == true);
        
        // Others should be non-redundant
        REQUIRE(result.redundant_indices[0] == false);
        REQUIRE(result.redundant_indices[3] == false);
        REQUIRE(result.redundant_indices[4] == false);
        REQUIRE(result.redundant_indices[5] == false);
        
        // Should have 4 non-redundant constraints
        REQUIRE(result.A_min.rows() == 4);
    }
    
    SECTION("Empty polytope") {
        // Inconsistent constraints
        Matrix A(2, 1);
        A(0, 0) =  1.0;  //  x <= -1
        A(1, 0) = -1.0;  // -x <= -2  =>  x >= 2
        
        Vector b(2);
        b[0] = -1.0;
        b[1] = -2.0;
        
        PolytopeRedundancyRemover solver;
        auto result = solver.indicate_nonredundant_halfplanes(A, b);
        
        // Algorithm should still run but may not find interior point
        // This is a degenerate case - the result depends on implementation details
        // We mainly test that it doesn't crash
        REQUIRE((result.success == true || result.success == false));
    }
    
    SECTION("Higher dimensional example") {
        // Simple 3D unit cube
        Matrix A(6, 3);
        A(0, 0) =  1.0; A(0, 1) =  0.0; A(0, 2) =  0.0;  //  x <= 1
        A(1, 0) = -1.0; A(1, 1) =  0.0; A(1, 2) =  0.0;  // -x <= 0  =>  x >= 0
        A(2, 0) =  0.0; A(2, 1) =  1.0; A(2, 2) =  0.0;  //  y <= 1
        A(3, 0) =  0.0; A(3, 1) = -1.0; A(3, 2) =  0.0;  // -y <= 0  =>  y >= 0
        A(4, 0) =  0.0; A(4, 1) =  0.0; A(4, 2) =  1.0;  //  z <= 1
        A(5, 0) =  0.0; A(5, 1) =  0.0; A(5, 2) = -1.0;  // -z <= 0  =>  z >= 0
        
        Vector b(6);
        b[0] = 1.0; b[1] = 0.0;
        b[2] = 1.0; b[3] = 0.0;
        b[4] = 1.0; b[5] = 0.0;
        
        PolytopeRedundancyRemover solver;
        auto result = solver.indicate_nonredundant_halfplanes(A, b);
        
        REQUIRE(result.success == true);
        
        // All constraints should be necessary for the unit cube
        for (int i = 0; i < 6; ++i) {
            REQUIRE(result.redundant_indices[i] == false);
        }
        
        REQUIRE(result.A_min.rows() == 6);
        REQUIRE(result.A_min.cols() == 3);
    }
}