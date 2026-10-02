#include <catch2/catch_all.hpp>
#include "polytope_redundancy/core.hpp"
#include <random>
#include <functional>

using namespace polytope_redundancy;

namespace {

Matrix make_matrix(const std::vector<std::vector<double>>& rows) {
    Matrix M(rows.size(), rows[0].size());
    for (size_t i = 0; i < rows.size(); ++i)
        for (size_t j = 0; j < rows[i].size(); ++j)
            M(i, j) = rows[i][j];
    return M;
}

// Ground truth by vertex enumeration: for polytopes in general position, a
// halfplane is nonredundant iff it is active at some vertex.
std::vector<bool> nonredundant_by_vertex_enumeration(const Matrix& A, const Vector& b) {
    const int m = A.rows(), n = A.cols();
    std::vector<bool> nonred(m, false);
    std::vector<int> idx(n);
    std::function<void(int, int)> recurse = [&](int start, int depth) {
        if (depth == n) {
            Matrix As(n, n);
            Vector x(n);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) As(i, j) = A(idx[i], j);
                x[i] = b[idx[i]];
            }
            std::vector<int> ipiv(n);
            int nrhs = 1, info = 0;
            dgesv_(&n, &nrhs, As.data(), &n, ipiv.data(), x.data(), &n, &info);
            if (info != 0) return;
            Vector Ax(m);
            A.gemv(x, Ax);
            for (int i = 0; i < m; ++i)
                if (Ax[i] > b[i] + 1e-9) return;
            for (int i = 0; i < n; ++i) nonred[idx[i]] = true;
            return;
        }
        for (int i = start; i < m; ++i) {
            idx[depth] = i;
            recurse(i + 1, depth + 1);
        }
    };
    recurse(0, 0);
    return nonred;
}

}  // namespace

TEST_CASE("Regression: QR column replacement", "[regression][qr]") {
    std::mt19937 gen(1);
    std::normal_distribution<> N;
    const int n = 6;

    SECTION("Square Q, every column position") {
        for (int p = 0; p < n; ++p) {
            Matrix A(n, n);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j) A(i, j) = N(gen);
            QRResult qr = qr_factorization(A);
            Vector a_new(n), a_old(n);
            for (int i = 0; i < n; ++i) { a_new[i] = N(gen); a_old[i] = A(i, p); A(i, p) = a_new[i]; }
            QRResult up = qr_column_replace(qr, a_new, a_old, p);
            REQUIRE(up.success);

            Matrix QR(n, n);
            up.Q.gemm(up.R, QR);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    REQUIRE(std::abs(QR(i, j) - A(i, j)) < 1e-12);
                    if (i > j) REQUIRE(std::abs(up.R(i, j)) < 1e-12);
                }
            }
        }
    }

    SECTION("Full Q of n x (n-1) matrix keeps last column in null space") {
        for (int p = 0; p < n - 1; ++p) {
            Matrix A(n, n - 1);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n - 1; ++j) A(i, j) = N(gen);
            QRResult qr = qr_factorization(A, /*full_q=*/true);
            Vector a_new(n), a_old(n);
            for (int i = 0; i < n; ++i) { a_new[i] = N(gen); a_old[i] = A(i, p); A(i, p) = a_new[i]; }
            QRResult up = qr_column_replace(qr, a_new, a_old, p);
            REQUIRE(up.success);
            for (int j = 0; j < n - 1; ++j) {
                double dot = 0.0;
                for (int i = 0; i < n; ++i) dot += up.Q(i, n - 1) * A(i, j);
                REQUIRE(std::abs(dot) < 1e-12);
            }
        }
    }
}

TEST_CASE("Regression: interior point shift", "[regression]") {
    // Square [1,3]^2 plus redundant x1 <= 5; origin is not interior
    Matrix A = make_matrix({{1, 0}, {-1, 0}, {0, 1}, {0, -1}, {1, 0}});
    Vector b(std::vector<double>{3, -1, 3, -1, 5});
    Vector z(std::vector<double>{2, 2});

    PolytopeRedundancyRemover solver;
    auto result = solver.indicate_nonredundant_halfplanes(A, b, {}, z);

    REQUIRE(result.success);
    REQUIRE(result.redundant_indices == std::vector<bool>{false, false, false, false, true});
    REQUIRE(result.b_min.size() == 4);
    for (int i = 0; i < 4; ++i) REQUIRE(result.b_min[i] == b[i]);
}

TEST_CASE("Regression: duplicate halfplanes", "[regression]") {
    // Like Matlab, only the last copy of a duplicate is kept
    Matrix A = make_matrix({{1, 0}, {-1, 0}, {0, 1}, {0, -1}, {1, 0}});
    Vector b(5, 1.0);

    PolytopeRedundancyRemover solver;
    auto result = solver.indicate_nonredundant_halfplanes(A, b);

    REQUIRE(result.success);
    REQUIRE(result.redundant_indices == std::vector<bool>{true, false, false, false, false});
    REQUIRE(result.A_min.rows() == 4);
}

TEST_CASE("Regression: symmetry detection", "[regression]") {
    SECTION("1D") {
        Matrix H = make_matrix({{-2}, {-1}, {1}, {2}});
        REQUIRE(make_set_symmetric(H, Vector(4, 1.0)).is_symmetric);
    }
    SECTION("2D hexagon") {
        Matrix H = make_matrix({{1, 0}, {0, 1}, {1, 1}, {-1, 0}, {0, -1}, {-1, -1}});
        REQUIRE(make_set_symmetric(H, Vector(6, 1.0)).is_symmetric);
    }
    SECTION("Symmetric set with duplicates") {
        // Square with x1 <= 1 and -x1 <= 1 duplicated, plus redundant +-x1 <= 2
        Matrix A = make_matrix({{1, 0}, {-1, 0}, {0, 1}, {0, -1}, {1, 0}, {-1, 0}, {2, 0}, {-2, 0}});
        Vector b(std::vector<double>{1, 1, 1, 1, 1, 1, 4, 4});

        PolytopeRedundancyRemover solver;
        auto result = solver.indicate_nonredundant_halfplanes(A, b);

        REQUIRE(result.success);
        REQUIRE(result.A_min.rows() == 4);
        // Kept rows must describe the square: each of +-e1, +-e2 exactly once
        for (int i = 0; i < 4; ++i) {
            double norm1 = std::abs(result.A_min(i, 0)) + std::abs(result.A_min(i, 1));
            REQUIRE(result.b_min[i] == Catch::Approx(norm1));
        }
    }
}

TEST_CASE("Regression: random polytopes vs vertex enumeration", "[regression][random]") {
    struct Config { int m, n; bool symmetric; };
    for (Config cfg : {Config{30, 2, false}, Config{60, 3, false}, Config{80, 4, false},
                       Config{40, 3, true}, Config{60, 4, true}}) {
        for (int seed = 0; seed < 10; ++seed) {
            std::mt19937 gen(1000 * cfg.m + seed);
            std::normal_distribution<> N;
            std::uniform_real_distribution<> U(0.5, 1.5);

            Matrix A(cfg.m, cfg.n);
            Vector b(cfg.m);
            const int rows = cfg.symmetric ? cfg.m / 2 : cfg.m;
            for (int i = 0; i < rows; ++i) {
                b[i] = cfg.symmetric ? 1.0 : U(gen);
                for (int j = 0; j < cfg.n; ++j) A(i, j) = N(gen);
                if (cfg.symmetric) {
                    b[i + rows] = b[i];
                    for (int j = 0; j < cfg.n; ++j) A(i + rows, j) = -A(i, j);
                }
            }

            PolytopeRedundancyRemover solver;
            auto result = solver.indicate_nonredundant_halfplanes(A, b);
            auto truth = nonredundant_by_vertex_enumeration(A, b);

            INFO("m=" << cfg.m << " n=" << cfg.n << " symmetric=" << cfg.symmetric << " seed=" << seed);
            REQUIRE(result.success);
            for (int i = 0; i < cfg.m; ++i) {
                INFO("row " << i);
                REQUIRE(result.redundant_indices[i] == !truth[i]);
            }
        }
    }
}
