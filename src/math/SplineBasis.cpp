#include "SplineBasis.hpp"

#include <cstdio>

SplineBasis::SplineBasis(const std::vector<double> &r_grid_, int k_spline, int n_spline)
: r_grid(r_grid_), num_spl(n_spline - 3)
{
    std::cout << "> Constructing BSpline basis, k = " << k_spline << ", n = " << n_spline - 3 << "... ";
    // construct splines
    BSpline splines(k_spline, n_spline, r_grid.front(), r_grid.back()); // k, n, r0, r_max

    auto N = r_grid.size();

    bspl = std::vector<std::vector<double>>(num_spl, std::vector<double>(N));
    bspl_derivative = std::vector<std::vector<double>>(num_spl, std::vector<double>(N));

    // reduce number of divisions.
    dr = r_grid.at(1) - r_grid.at(0); // is this correct?
    auto dr2 = dr / 2.0;

    #pragma omp parallel for
    for (int i = 2; i < n_spline - 1; ++i)
    {
        for (std::size_t j = 0; j < N; ++j)
        {
            bspl.at(i-2).at(j) = splines.b(i, r_grid.at(j));
            // one less division each time?
            bspl_derivative.at(i-2).at(j) = (splines.b(i, r_grid.at(j) + dr2) - splines.b(i, r_grid.at(j) - dr2)) / dr;
        }
    }

    std::cout << "done.\n";
}