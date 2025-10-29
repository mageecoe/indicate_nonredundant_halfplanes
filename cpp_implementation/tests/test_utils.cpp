#include <catch2/catch_all.hpp>
#include "polytope_redundancy/matrix_utils.hpp"
#include "polytope_redundancy/geometry.hpp"

using namespace polytope_redundancy;

TEST_CASE("Vector operations", "[vector]") {
    SECTION("Construction and basic operations") {
        Vector v1(3, 2.0);
        REQUIRE(v1.size() == 3);
        REQUIRE(v1[0] == 2.0);
        REQUIRE(v1[1] == 2.0);
        REQUIRE(v1[2] == 2.0);
        
        Vector v2(3);
        v2[0] = 1.0;
        v2[1] = 2.0;
        v2[2] = 3.0;
        
        REQUIRE(v2.norm() == Approx(std::sqrt(14.0)).epsilon(1e-10));
    }
    
    SECTION("Dot product") {
        Vector v1(3);
        v1[0] = 1.0; v1[1] = 2.0; v1[2] = 3.0;
        
        Vector v2(3);
        v2[0] = 4.0; v2[1] = 5.0; v2[2] = 6.0;
        
        REQUIRE(v1.dot(v2) == Approx(32.0).epsilon(1e-10));
    }
    
    SECTION("AXPY operation") {
        Vector v1(3);
        v1[0] = 1.0; v1[1] = 2.0; v1[2] = 3.0;
        
        Vector v2(3);
        v2[0] = 2.0; v2[1] = 3.0; v2[2] = 4.0;
        
        v1.axpy(2.0, v2); // v1 = v1 + 2*v2
        
        REQUIRE(v1[0] == Approx(5.0).epsilon(1e-10));
        REQUIRE(v1[1] == Approx(8.0).epsilon(1e-10));
        REQUIRE(v1[2] == Approx(11.0).epsilon(1e-10));
    }
}

TEST_CASE("Matrix operations", "[matrix]") {
    SECTION("Construction and indexing") {
        Matrix A(2, 3);
        A(0, 0) = 1.0; A(0, 1) = 2.0; A(0, 2) = 3.0;
        A(1, 0) = 4.0; A(1, 1) = 5.0; A(1, 2) = 6.0;
        
        REQUIRE(A.rows() == 2);
        REQUIRE(A.cols() == 3);
        REQUIRE(A(0, 1) == 2.0);
        REQUIRE(A(1, 2) == 6.0);
    }
    
    SECTION("Matrix-vector multiplication") {
        Matrix A(2, 3);
        A(0, 0) = 1.0; A(0, 1) = 2.0; A(0, 2) = 3.0;
        A(1, 0) = 4.0; A(1, 1) = 5.0; A(1, 2) = 6.0;
        
        Vector x(3);
        x[0] = 1.0; x[1] = 2.0; x[2] = 3.0;
        
        Vector y(2);
        A.gemv(x, y);
        
        // Expected: [1*1 + 2*2 + 3*3, 4*1 + 5*2 + 6*3] = [14, 32]
        REQUIRE(y[0] == Approx(14.0).epsilon(1e-10));
        REQUIRE(y[1] == Approx(32.0).epsilon(1e-10));
    }
    
    SECTION("Row norms") {
        Matrix A(2, 3);
        A(0, 0) = 3.0; A(0, 1) = 4.0; A(0, 2) = 0.0; // Norm = 5
        A(1, 0) = 1.0; A(1, 1) = 1.0; A(1, 2) = 1.0; // Norm = sqrt(3)
        
        Vector norms = A.row_norms();
        
        REQUIRE(norms[0] == Approx(5.0).epsilon(1e-10));
        REQUIRE(norms[1] == Approx(std::sqrt(3.0)).epsilon(1e-10));
    }
}

TEST_CASE("Matrix utilities", "[matrix_utils]") {
    SECTION("Normalize halfplane description") {
        Matrix H(2, 2);
        H(0, 0) = 3.0; H(0, 1) = 4.0;
        H(1, 0) = 6.0; H(1, 1) = 8.0;
        
        Vector h(2);
        h[0] = 10.0;
        h[1] = 20.0;
        
        auto [H_norm, h_norm] = normalize_halfplane_description(H, h, false);
        
        // After normalization by |h|: first row should be [3/10, 4/10] and h[0] = 1
        REQUIRE(H_norm(0, 0) == Approx(0.3).epsilon(1e-10));
        REQUIRE(H_norm(0, 1) == Approx(0.4).epsilon(1e-10));
        REQUIRE(h_norm[0] == Approx(1.0).epsilon(1e-10));
        
        REQUIRE(H_norm(1, 0) == Approx(0.3).epsilon(1e-10));
        REQUIRE(H_norm(1, 1) == Approx(0.4).epsilon(1e-10));
        REQUIRE(h_norm[1] == Approx(1.0).epsilon(1e-10));
    }
    
    SECTION("Unique with tolerance") {
        Matrix A(4, 2);
        A(0, 0) = 1.0; A(0, 1) = 2.0;
        A(1, 0) = 3.0; A(1, 1) = 4.0;
        A(2, 0) = 1.0 + 1e-6; A(2, 1) = 2.0 + 1e-6; // Should be considered duplicate
        A(3, 0) = 5.0; A(3, 1) = 6.0;
        
        auto [A_unique, is_unique] = unique_with_tolerance(A, 1e-5);
        
        // Should have 3 unique rows (removing the near-duplicate)
        REQUIRE(A_unique.rows() == 3);
        
        // Check which rows are considered unique
        int unique_count = 0;
        for (bool unique : is_unique) {
            if (unique) unique_count++;
        }
        REQUIRE(unique_count == 3);
    }
}

TEST_CASE("Geometry utilities", "[geometry]") {
    SECTION("Basic feasible solution finder") {
        // Simple 2D example: unit square [0,1] x [0,1]
        Matrix A(4, 2);
        A(0, 0) = -1.0; A(0, 1) =  0.0; // -x >= -1  =>  x <= 1
        A(1, 0) =  1.0; A(1, 1) =  0.0; //  x >= 0   =>  -x <= 0
        A(2, 0) =  0.0; A(2, 1) = -1.0; // -y >= -1  =>  y <= 1
        A(3, 0) =  0.0; A(3, 1) =  1.0; //  y >= 0   =>  -y <= 0
        
        Vector b(4);
        b[0] = 1.0;  // x <= 1
        b[1] = 0.0;  // -x <= 0  => x >= 0
        b[2] = 1.0;  // y <= 1
        b[3] = 0.0;  // -y <= 0  => y >= 0
        
        Vector interior_point(2);
        interior_point[0] = 0.5;
        interior_point[1] = 0.5;
        
        BFSResult result = find_basic_feasible_solution(A, b, interior_point);
        
        REQUIRE(result.success == true);
        REQUIRE(result.x.size() == 2);
        
        // Verify feasibility
        Vector Ax(4);
        A.gemv(result.x, Ax);
        for (int i = 0; i < 4; ++i) {
            REQUIRE(Ax[i] <= b[i] + 1e-10);
        }
    }
}