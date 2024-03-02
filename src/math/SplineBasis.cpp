#include "SplineBasis.hpp"

#include <cstdio>

#include <limits> // std::epsilon()

SplineBasis::SplineBasis(const std::vector<double> &r_grid_, int k_spline_, int n_spline_)
:   num_spl(n_spline_ - 3),
    k_spline(k_spline_),
    r_grid(r_grid_),
    dr(r_grid.at(1) - r_grid.at(0)),
    bspl(vector_set(num_spl, std::vector<double>(r_grid.size()))),
    bspl_derivative(vector_set(num_spl, std::vector<double>(r_grid.size()))),
    Bmatrix(SquareMatrix<double>(num_spl))
{
    IO::msg::construct<int>("B-Spline basis", {{"k", k_spline}, {"n", num_spl}}, true);

    construct_spline_vectors();
    construct_matrix();

    IO::msg::done(true);
}

auto SplineBasis::construct_spline_vectors() -> void
{
    IO::msg::construct("basis vectors");
    
    BSpline splines(k_spline, num_spl + 3, r_grid.front(), r_grid.back());

    auto h = sqrt(std::numeric_limits<double>::epsilon());

    #pragma omp parallel for
    for (int i = 2; i < num_spl + 2; ++i)
    {
        for (std::size_t j = 0; j < r_grid.size(); ++j)
        {
            bspl.at(i-2).at(j) = splines.b(i, r_grid.at(j));
            
            // Dividing by the difference rather than 2.0 * h is mathematically equivalent
            // but reduces error due to subtraction of small numbers.
            bspl_derivative.at(i-2).at(j) =
            (splines.b(i, r_grid.at(j) + h) - splines.b(i, r_grid.at(j) - h)) / ((r_grid.at(j)+h) - ((r_grid.at(j)-h)));
        }
    }

    IO::msg::done();
}

auto SplineBasis::construct_matrix() -> void
{
    IO::msg::construct("B-matrix");

    #pragma omp parallel for
    for(int i = 0; i < num_spl; ++i)
    {
        for(int j = 0; j <= i; ++j)
        {
            Bmatrix(i,j) = trapz_linear(dr, bspl.at(i) * bspl.at(j));

            // Only need to calculate bottom half.
            Bmatrix(j,i) = Bmatrix(i,j);
        }
    }

    IO::msg::done();
}