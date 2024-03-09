#include "matrix.hpp"

template class Matrix<int>;
template class Matrix<double>;


template <typename T>
Matrix<T>::Matrix()
: size_x(0), size_y(0), m_data(std::vector<T>()) {};


template <typename T>
Matrix<T>::Matrix(const std::size_t size_x_, const std::size_t size_y_)
: size_x(size_x_), size_y(size_y_), m_data(std::vector<T>(size_x_ * size_y_)) {};

template <typename T>
Matrix<T>::Matrix(const std::size_t size_x_, const std::size_t size_y_, std::vector<T> data_)
: size_x(size_x_), size_y(size_y_), m_data(data_) {};

template <typename T>
inline
auto Matrix<T>::operator()(std::size_t x, std::size_t y) -> T&
{ 
    return m_data.at(x * size_y + y);
}

template <typename T>
inline
auto Matrix<T>::operator()(std::size_t x, std::size_t y) const -> T
{
    return m_data.at(x * size_y + y);
}

template <typename T>
auto Matrix<T>::operator+=(const Matrix<T> &rhs) -> Matrix<T>&
{
    assert(this->size_x == rhs.size_x && this->size_y == rhs.size_y);
    for (std::size_t i = 0; i < rhs.size_x; ++i)
    {
        for (std::size_t j = 0; j < rhs.size_y; ++j)
        {
            (*this)(i, j) += rhs(i, j);
        }
    }
    return *this;
}

template <typename T>
auto Matrix<T>::operator-=(const Matrix<T> &rhs) -> Matrix<T>&
{
    assert(this->size_x == rhs.size_x && this->size_y == rhs.size_y);
    for (std::size_t i = 0; i < rhs.size_x; ++i)
    {
        for (std::size_t j = 0; j < rhs.size_x; ++j)
        {
            (*this)(i, j) -= rhs(i, j);
        }
    }
    return *this;
}
template <typename T>
auto Matrix<T>::operator*=(const T rhs) -> Matrix<T>&
{
    for (std::size_t i = 0; i < this->size_x; ++i)
    {
        for (std::size_t j = 0; j < this->size_y; ++j)
        {
            (*this)(i, j) *= rhs;
        }
    }
    return *this;
}

template <typename T>
auto operator+(Matrix<T> lhs, const Matrix<T> &rhs) -> Matrix<T> { return lhs += rhs; }

template <typename T>
auto operator-(Matrix<T> lhs, const Matrix<T> &rhs) -> Matrix<T> { return lhs -= rhs; }

template <typename T>
auto operator* (Matrix<T> lhs, T rhs) -> Matrix<T> { return lhs *= rhs; };

template <typename T>
auto operator* (T lhs, Matrix<T> rhs) -> Matrix<T> { return rhs *= lhs; };

template <typename T>
auto Matrix<T>::get_row(std::size_t row) -> std::vector<T>
{
    std::vector<T> result(size_y);
    for (std::size_t i = 0; i < size_y; ++i)
        result.at(i) = (*this)(row, i);
    return result;
}


// uses LAPACK function DYSGV
// solves problems of the form Av = eBv.
auto MatrixTools::solve_eigen_system(SquareMatrix<double> A, SquareMatrix<double> B) -> std::pair<SquareMatrix<double>, std::vector<double> >
{
    IO::log("Calling LAPACK subroutine DYSGV");
    
    int  itype { 1 };
    char jobz  {'V'};
    char uplo  {'U'};

    int N = static_cast<int>(A.get_size());

    std::vector<double> eigenvalues(A.get_size());

    // blank array work of length lwork - memory used by LAPACK. We use the recommended value of lwork.
    int lwork = 6 * N;
    //double *work = new double[lwork];
    std::vector<double> work(6 * A.get_size());

    // error code. info = 0 if successful.
    int info = 0;

    dsygv_(&itype, &jobz, &uplo, &N, A.data(), &N , B.data(), &N, eigenvalues.data(), work.data(), &lwork , &info);

    // Best way to do this?
    if (info != 0)
    {
        IO::log_params(LogType::error, "Error: dysgv_ failed", {{"INFO", info}});
    }
    else
    {
        IO::done();
    }

    //delete [] work;

    return {A, eigenvalues};
}

/*
auto inline MatrixTools::innerProduct(double dr, std::vector<double> bra, std::vector<double> ket) -> double
{
    return trapz_linear(r_grid, bra * ket);
}


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