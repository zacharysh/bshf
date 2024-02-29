#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include <cstddef>
#include <iostream>

#include <cassert> // assert

#include <vector>

#include "math_other.hpp"

// dsygv_ is a symbol in the LAPACK library files.
// Documentation: http://www.netlib.org/lapack/explore-html/index.html
extern "C"
{
    int dsygv_(int *ITYPE, char *JOBZ, char *UPLO,
            int *N, double *A, int *LDA, double *B, int *LDB,
            double *W, double *WORK , int *LWORK,  int *INFO);
}

template <typename T>
class Matrix
{
    public:
    const std::size_t size_x;
    const std::size_t size_y;
    T* data; // row major order

    Matrix(int size_x_, int size_y_);
    Matrix(std::size_t size_x_, std::size_t size_y_);
    

    ~Matrix();                                     // Destructor
    Matrix(const Matrix<T>& M);                    // Copy constructor

    auto operator()(std::size_t x, std::size_t y) -> T&;
    auto operator()(std::size_t x, std::size_t y) const -> T;

    auto operator= (const Matrix<T> &rhs) -> Matrix<T>&;
    auto operator+=(const Matrix<T> &rhs) -> Matrix<T>&;
    auto operator-=(const Matrix<T> &rhs) -> Matrix<T>&;
    auto operator*=(const T rhs)          -> Matrix<T>&;

    friend auto operator<<(std::ostream& os, const Matrix<T> &matrix) -> std::ostream&
    {
        for(std::size_t i = 0; i < matrix.size_x; ++i)
        {
            os << "|";
            for(std::size_t j = 0; j < matrix.size_y - 1; ++j)
            {
                //os << matrix(i,j) << ", ";
                printf("%+1.4e, ", matrix(i,j));
            }
            os << matrix(i, matrix.size_y - 1) << "|\n";
        }
        return os;
    }

    auto transpose() -> Matrix<T>;
};

template <typename T>
class SquareMatrix : public Matrix<T>
{
    public:
    SquareMatrix(int size_)
    : Matrix<T>(size_, size_) {}

    auto get_size() -> std::size_t { return this->size_x; }
};

namespace MatrixTools
{
    auto solveEigenSystem(SquareMatrix<double> &A, SquareMatrix<double> &B) -> std::vector<double>;

    auto inline innerProduct(std::vector<double> r_grid, std::vector<double> bra, std::vector<double> ket) -> double;

    //auto computeMatrixElements(std::vector<double> r_grid, std::vector<std::vector<double>> bra, std::vector<std::vector<double>> ket)
    //-> SquareMatrix<double>*;
}; // namespace MatrixTools

#endif