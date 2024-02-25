#include "matrix.hpp"

template class Matrix<int>;
template class Matrix<double>;

template <typename T>
Matrix<T>::Matrix(int size_x_, int size_y_)
: size_x(size_x_), size_y(size_y_)
{
    int length = size_x_ * size_y_;
    data = new T[length]; // size_t?
    
    for (std::size_t i = 0; i < size_x; ++i)
    {
        for (std::size_t j = 0; j < size_y; ++j)
        {
            (*this)(i,j) = 0.0;
        }
    }
}

template <typename T>
Matrix<T>::Matrix(std::size_t size_x_, std::size_t size_y_)
: size_x(size_x_), size_y(size_y_)
{
    int length = static_cast<int>(size_x_ * size_y_);
    data = new T[length]; // size_t?
    
    // don't like repeating code...
    for (std::size_t i = 0; i < size_x; ++i)
    {
        for (std::size_t j = 0; j < size_y; ++j)
        {
            (*this)(i,j) = 0.0;
        }
    }
}

template <typename T>
inline
auto Matrix<T>::operator()(std::size_t x, std::size_t y) -> T&
{ 
    return data[x * size_y + y];
}

template <typename T>
inline
auto Matrix<T>::operator()(std::size_t x, std::size_t y) const -> T
{
    return data[x * size_y + y];
}

template <typename T>
Matrix<T>::~Matrix()
{
    delete [] data;
}


template <typename T>
auto Matrix<T>::transpose() -> Matrix<T>
{
    Matrix<T> result(size_y, size_x);

    for (std::size_t i = 0; i < size_x; ++i)
    {
        for (std::size_t j = 0; j < size_y; ++j)
        {
            result(j,i) = (*this)(i,j);
        }
    }
    return result;
}

// uses LAPACK function DYSGV
// solves problems of the form Av = eBv.
auto MatrixTools::solveEigenSystem(SquareMatrix<double> *A, SquareMatrix<double> *B) -> double*
{
    std::cout << "Calling LAPACK subroutine DYSGV...";
    int itype = 1;
    char jobz = 'V';
    char uplo = 'U';
    int N = static_cast<int>(A->get_size());
    double *eigenvalues = new double[N];

    // blank array work of length lwork - memory used by LAPACK. We use the recommended value.
    int lwork = 6 * N;
    double *work = new double[lwork];

    // error code. info = 0 if successful.
    int info = 0;

    dsygv_(&itype, &jobz, &uplo, &N, A->data, &N , B->data, &N, eigenvalues, work, &lwork , &info);

    std::cout << " done. (INFO = " << info << ").\n\n";

    delete [] work;

    return eigenvalues;
}

auto MatrixTools::innerProduct(std::vector<double> r_grid, std::vector<double> bra, std::vector<double> ket) -> double
{
    return integrate_trap(r_grid, bra * ket);
}

/*
auto computeMatrixElements(std::vector<double> r_grid, std::vector<std::vector<double>> bra, std::vector<std::vector<double>> ket)
-> SquareMatrix<double>*
{
    // not sure if this is best way to do it...
    size_t N = bra.size();
    SquareMatrix<double> *matrix = new SquareMatrix<double>(N);

    for(int i = 0; i < N; ++i)
    {
        for(int j = 0; j <= i; ++j)
        {
            (*matrix)(i, j) = integrate_trap(r_grid, bra.at(i) * ket.at(j));
        }
    }
    return matrix;
}
*/