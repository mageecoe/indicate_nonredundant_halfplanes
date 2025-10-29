#pragma once

#include <memory>
#include <vector>
#include <stdexcept>

// Forward declarations for BLAS/LAPACK
extern "C" {
    // BLAS Level 1
    double dnrm2_(const int* n, const double* x, const int* incx);
    double ddot_(const int* n, const double* x, const int* incx, const double* y, const int* incy);
    void daxpy_(const int* n, const double* alpha, const double* x, const int* incx, double* y, const int* incy);
    void dscal_(const int* n, const double* alpha, double* x, const int* incx);
    void dcopy_(const int* n, const double* x, const int* incx, double* y, const int* incy);
    
    // BLAS Level 2
    void dgemv_(const char* trans, const int* m, const int* n, const double* alpha,
                const double* a, const int* lda, const double* x, const int* incx,
                const double* beta, double* y, const int* incy);
    
    // BLAS Level 3
    void dgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
                const double* alpha, const double* a, const int* lda, const double* b, const int* ldb,
                const double* beta, double* c, const int* ldc);
    
    // LAPACK
    void dgesv_(const int* n, const int* nrhs, double* a, const int* lda, int* ipiv,
                double* b, const int* ldb, int* info);
    void dgeqrf_(const int* m, const int* n, double* a, const int* lda, double* tau,
                 double* work, const int* lwork, int* info);
    void dorgqr_(const int* m, const int* n, const int* k, double* a, const int* lda,
                 const double* tau, double* work, const int* lwork, int* info);
}

namespace polytope_redundancy {

class Vector {
private:
    std::unique_ptr<double[]> data_;
    int size_;
    bool owns_data_;
    double* raw_ptr_;

public:
    // Constructors
    Vector() : data_(nullptr), size_(0), owns_data_(false), raw_ptr_(nullptr) {}
    
    explicit Vector(int size) 
        : data_(std::make_unique<double[]>(size)), size_(size), owns_data_(true), raw_ptr_(data_.get()) {
        std::fill_n(raw_ptr_, size_, 0.0);
    }
    
    Vector(int size, double value) 
        : data_(std::make_unique<double[]>(size)), size_(size), owns_data_(true), raw_ptr_(data_.get()) {
        std::fill_n(raw_ptr_, size_, value);
    }
    
    Vector(const std::vector<double>& vec) 
        : data_(std::make_unique<double[]>(vec.size())), size_(static_cast<int>(vec.size())), 
          owns_data_(true), raw_ptr_(data_.get()) {
        std::copy(vec.begin(), vec.end(), raw_ptr_);
    }
    
    // Non-owning constructor (wraps existing data)
    Vector(double* data, int size) 
        : data_(nullptr), size_(size), owns_data_(false), raw_ptr_(data) {}
    
    // Copy constructor
    Vector(const Vector& other) 
        : data_(std::make_unique<double[]>(other.size_)), size_(other.size_), 
          owns_data_(true), raw_ptr_(data_.get()) {
        const int inc = 1;
        dcopy_(&size_, other.raw_ptr_, &inc, raw_ptr_, &inc);
    }
    
    // Move constructor
    Vector(Vector&& other) noexcept 
        : data_(std::move(other.data_)), size_(other.size_), 
          owns_data_(other.owns_data_), raw_ptr_(other.raw_ptr_) {
        other.raw_ptr_ = nullptr;
        other.size_ = 0;
    }
    
    // Assignment operators
    Vector& operator=(const Vector& other) {
        if (this != &other) {
            if (size_ != other.size_) {
                data_ = std::make_unique<double[]>(other.size_);
                size_ = other.size_;
                owns_data_ = true;
                raw_ptr_ = data_.get();
            }
            const int inc = 1;
            dcopy_(&size_, other.raw_ptr_, &inc, raw_ptr_, &inc);
        }
        return *this;
    }
    
    Vector& operator=(Vector&& other) noexcept {
        if (this != &other) {
            data_ = std::move(other.data_);
            size_ = other.size_;
            owns_data_ = other.owns_data_;
            raw_ptr_ = other.raw_ptr_;
            other.raw_ptr_ = nullptr;
            other.size_ = 0;
        }
        return *this;
    }
    
    // Accessors
    int size() const { return size_; }
    double* data() { return raw_ptr_; }
    const double* data() const { return raw_ptr_; }
    
    double& operator[](int i) { return raw_ptr_[i]; }
    const double& operator[](int i) const { return raw_ptr_[i]; }
    
    // BLAS operations
    double norm() const {
        const int inc = 1;
        return dnrm2_(&size_, raw_ptr_, &inc);
    }
    
    double dot(const Vector& other) const {
        if (size_ != other.size_) {
            throw std::invalid_argument("Vector sizes must match for dot product");
        }
        const int inc = 1;
        return ddot_(&size_, raw_ptr_, &inc, other.raw_ptr_, &inc);
    }
    
    void axpy(double alpha, const Vector& x) {
        if (size_ != x.size_) {
            throw std::invalid_argument("Vector sizes must match for axpy operation");
        }
        const int inc = 1;
        daxpy_(&size_, &alpha, x.raw_ptr_, &inc, raw_ptr_, &inc);
    }
    
    void scale(double alpha) {
        const int inc = 1;
        dscal_(&size_, &alpha, raw_ptr_, &inc);
    }
    
    // Utility functions
    void fill(double value) {
        std::fill_n(raw_ptr_, size_, value);
    }
    
    std::vector<double> to_vector() const {
        return std::vector<double>(raw_ptr_, raw_ptr_ + size_);
    }
};

class Matrix {
private:
    std::unique_ptr<double[]> data_;
    int rows_, cols_;
    bool owns_data_;
    double* raw_ptr_;

public:
    // Constructors
    Matrix() : data_(nullptr), rows_(0), cols_(0), owns_data_(false), raw_ptr_(nullptr) {}
    
    Matrix(int rows, int cols) 
        : data_(std::make_unique<double[]>(rows * cols)), rows_(rows), cols_(cols),
          owns_data_(true), raw_ptr_(data_.get()) {
        std::fill_n(raw_ptr_, rows_ * cols_, 0.0);
    }
    
    Matrix(int rows, int cols, double value) 
        : data_(std::make_unique<double[]>(rows * cols)), rows_(rows), cols_(cols),
          owns_data_(true), raw_ptr_(data_.get()) {
        std::fill_n(raw_ptr_, rows_ * cols_, value);
    }
    
    // Non-owning constructor
    Matrix(double* data, int rows, int cols) 
        : data_(nullptr), rows_(rows), cols_(cols), owns_data_(false), raw_ptr_(data) {}
    
    // Copy constructor
    Matrix(const Matrix& other) 
        : data_(std::make_unique<double[]>(other.rows_ * other.cols_)), 
          rows_(other.rows_), cols_(other.cols_), owns_data_(true), raw_ptr_(data_.get()) {
        const int size = rows_ * cols_;
        const int inc = 1;
        dcopy_(&size, other.raw_ptr_, &inc, raw_ptr_, &inc);
    }
    
    // Move constructor
    Matrix(Matrix&& other) noexcept 
        : data_(std::move(other.data_)), rows_(other.rows_), cols_(other.cols_),
          owns_data_(other.owns_data_), raw_ptr_(other.raw_ptr_) {
        other.raw_ptr_ = nullptr;
        other.rows_ = other.cols_ = 0;
    }
    
    // Assignment operators
    Matrix& operator=(const Matrix& other) {
        if (this != &other) {
            if (rows_ * cols_ != other.rows_ * other.cols_) {
                data_ = std::make_unique<double[]>(other.rows_ * other.cols_);
                owns_data_ = true;
                raw_ptr_ = data_.get();
            }
            rows_ = other.rows_;
            cols_ = other.cols_;
            const int size = rows_ * cols_;
            const int inc = 1;
            dcopy_(&size, other.raw_ptr_, &inc, raw_ptr_, &inc);
        }
        return *this;
    }
    
    Matrix& operator=(Matrix&& other) noexcept {
        if (this != &other) {
            data_ = std::move(other.data_);
            rows_ = other.rows_;
            cols_ = other.cols_;
            owns_data_ = other.owns_data_;
            raw_ptr_ = other.raw_ptr_;
            other.raw_ptr_ = nullptr;
            other.rows_ = other.cols_ = 0;
        }
        return *this;
    }
    
    // Accessors
    int rows() const { return rows_; }
    int cols() const { return cols_; }
    double* data() { return raw_ptr_; }
    const double* data() const { return raw_ptr_; }
    
    // Column-major indexing (FORTRAN style)
    double& operator()(int i, int j) { return raw_ptr_[j * rows_ + i]; }
    const double& operator()(int i, int j) const { return raw_ptr_[j * rows_ + i]; }
    
    // BLAS operations
    void gemv(const Vector& x, Vector& y, double alpha = 1.0, double beta = 0.0, bool transpose = false) const {
        if (transpose) {
            if (rows_ != x.size() || cols_ != y.size()) {
                throw std::invalid_argument("Matrix and vector dimensions incompatible for transposed gemv");
            }
        } else {
            if (cols_ != x.size() || rows_ != y.size()) {
                throw std::invalid_argument("Matrix and vector dimensions incompatible for gemv");
            }
        }
        
        const char trans = transpose ? 'T' : 'N';
        const int inc = 1;
        dgemv_(&trans, &rows_, &cols_, &alpha, raw_ptr_, &rows_, 
               x.data(), &inc, &beta, y.data(), &inc);
    }
    
    void gemm(const Matrix& B, Matrix& C, double alpha = 1.0, double beta = 0.0, 
              bool transpose_A = false, bool transpose_B = false) const {
        const char transa = transpose_A ? 'T' : 'N';
        const char transb = transpose_B ? 'T' : 'N';
        
        const int m = transpose_A ? cols_ : rows_;
        const int n = transpose_B ? B.rows_ : B.cols_;
        const int k = transpose_A ? rows_ : cols_;
        
        if (C.rows_ != m || C.cols_ != n) {
            throw std::invalid_argument("Result matrix dimensions incorrect for gemm");
        }
        
        dgemm_(&transa, &transb, &m, &n, &k, &alpha, raw_ptr_, &rows_, 
               B.raw_ptr_, &B.rows_, &beta, C.raw_ptr_, &C.rows_);
    }
    
    // Get row as non-owning vector
    Vector row(int i) {
        if (i < 0 || i >= rows_) {
            throw std::out_of_range("Row index out of range");
        }
        // For column-major storage, we need to stride through the matrix
        // This creates a view with stride = rows_
        // Note: This is not directly BLAS-compatible due to stride
        // For now, we'll copy the row
        Vector result(cols_);
        for (int j = 0; j < cols_; ++j) {
            result[j] = (*this)(i, j);
        }
        return result;
    }
    
    // Const version of row
    Vector row(int i) const {
        if (i < 0 || i >= rows_) {
            throw std::out_of_range("Row index out of range");
        }
        Vector result(cols_);
        for (int j = 0; j < cols_; ++j) {
            result[j] = (*this)(i, j);
        }
        return result;
    }
    
    // Get column as non-owning vector
    Vector col(int j) {
        if (j < 0 || j >= cols_) {
            throw std::out_of_range("Column index out of range");
        }
        return Vector(raw_ptr_ + j * rows_, rows_);
    }
    
    // Utility functions
    void fill(double value) {
        std::fill_n(raw_ptr_, rows_ * cols_, value);
    }
    
    Vector row_norms() const {
        Vector norms(rows_);
        for (int i = 0; i < rows_; ++i) {
            double sum = 0.0;
            for (int j = 0; j < cols_; ++j) {
                const double val = (*this)(i, j);
                sum += val * val;
            }
            norms[i] = std::sqrt(sum);
        }
        return norms;
    }
    
    // Static method to create identity matrix
    static Matrix identity(int n) {
        Matrix I(n, n);
        for (int i = 0; i < n; ++i) {
            I(i, i) = 1.0;
        }
        return I;
    }
};

// Function declarations for matrix utilities
std::pair<Matrix, Vector> normalize_halfplane_description(const Matrix& H, const Vector& h, bool normalize_H = false);

struct SymmetricResult {
    Matrix H;
    Vector h;
    bool is_symmetric;
    std::vector<int> permutation;
};

SymmetricResult make_set_symmetric(const Matrix& H, const Vector& h, double tol = 1e-6);

std::pair<Matrix, std::vector<bool>> unique_with_tolerance(const Matrix& A, double tol = 1e-5);

Vector row_norms(const Matrix& A);

// QR decomposition functionality
struct QRResult {
    Matrix Q;
    Matrix R;
    Vector tau;  // Householder reflector scalars
    bool success;
};

// Compute QR factorization of matrix A
QRResult qr_factorization(const Matrix& A);

// Solve linear system using QR factorization: solve R*x = Q'*b  
Vector qr_solve(const QRResult& qr, const Vector& b);

// Update QR factorization when adding/removing a row (simplified version)
QRResult qr_update(const QRResult& qr_old, const Vector& new_row, bool add_row = true);

} // namespace polytope_redundancy