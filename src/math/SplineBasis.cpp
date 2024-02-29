#include "SplineBasis.hpp"

#include <cstdio>

SplineBasis::SplineBasis(const std::vector<double> &r_grid_, int k_spline, int n_spline)
: r_grid(r_grid_), num_spl(n_spline - 3)
{
    std::cout << "> Constructing BSpline basis, k = " << k_spline << ", n = " << n_spline - 3 << "... ";
    // construct splines
    BSpline splines(k_spline, n_spline, r_grid.front(), r_grid.back()); // k, n, r0, r_max

    std::size_t N = r_grid.size();

    bspl.resize(num_spl, std::vector<double>(N));
    bspl_derivative.resize(num_spl, std::vector<double>(N));

    // for linear grid, this is constant
    auto dr = (r_grid.at(1) - r_grid.at(0))/2.0;

    #pragma omp parallel for
    for (int i = 2; i < n_spline - 1; ++i)
    {
        for (std::size_t j = 0; j < N; ++j)
        {
            bspl.at(i-2).at(j) = splines.b(i, r_grid.at(j));
            // one less division each time?
            bspl_derivative.at(i-2).at(j) = (splines.b(i, r_grid.at(j) + dr) - splines.b(i, r_grid.at(j) - dr)) / (2.0 * dr);
        }
    }
    std::cout << "done.\n";
}