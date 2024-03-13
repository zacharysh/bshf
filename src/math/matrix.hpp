#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include <cstddef>
#include <iostream>
#include <cassert> // assert
#include <vector>
//#include <utility> // pair ?
#include <memory> // std::unique_ptr

#include "math_other.hpp"

#include "../io.hpp"

// dsygv_ is a symbol in the LAPACK library files.
// Documentation: http://www.netlib.org/lapack/explore-html/index.html
extern "C"
{
    void dsygv_(int *ITYPE, char *JOBZ, char *UPLO,
            int *N, double *A, int *LDA, double *B, int *LDB,
            double *W, double *WORK , int *LWORK,  int *INFO);
}

template <typename T = double>
class Matrix
{
    public:
    
    std::size_t size_x;
    std::size_t size_y;
    std::vector<T> m_data; // row major order

    Matrix();
    Matrix(const std::size_t size_x_, const std::size_t size_y_);
    Matrix(const std::size_t size_x_, const std::size_t size_y_, std::vector<T> data_);

    Matrix(const Matrix<T> &) = default;
    Matrix(Matrix<T> &&) = default;
    auto operator=(const Matrix<T> &) -> Matrix<T>& = default;
    auto operator=(Matrix<T> &&) -> Matrix<T>& = default;

    auto operator()(std::size_t x, std::size_t y) -> T&;
    auto operator()(std::size_t x, std::size_t y) const -> T;

    //auto operator= (const Matrix<T> &rhs) -> Matrix<T>&;
    auto operator+=(const Matrix<T> &rhs) -> Matrix<T>&;
    auto operator-=(const Matrix<T> &rhs) -> Matrix<T>&;
    auto operator*=(const T rhs)          -> Matrix<T>&;

    auto data() -> T* { return m_data.data(); }

    friend auto operator<<(std::ostream &os, const Matrix<T> &matrix) -> std::ostream&
    {
        for(std::size_t i = 0; i < matrix.m_size_x; ++i)
        {
            os << "|";
            for(std::size_t j = 0; j < matrix.m_size_y - 1; ++j)
            {
                //os << matrix(i,j) << ", ";
                printf("%+1.4e, ", matrix(i,j));
            }
            os << matrix(i, matrix.m_size_y - 1) << "|\n";
        }
        return os;
    }

    auto get_row(std::size_t row) -> std::vector<T>;
};

template <typename T>
class SquareMatrix : public Matrix<T>
{
    public:
    SquareMatrix(std::size_t size_)
    : Matrix<T>(size_, size_) {};

    SquareMatrix() : Matrix<T>() {};

    auto get_size() -> std::size_t { return this->size_x; }
};

namespace MatrixTools
{
    auto solve_eigen_system(SquareMatrix<double> A, SquareMatrix<double> B) -> std::pair<SquareMatrix<double>, std::vector<double> >;

}; // namespace MatrixTools

#endif