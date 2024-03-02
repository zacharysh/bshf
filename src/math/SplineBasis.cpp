#include "SplineBasis.hpp"

#include <cstdio>

#include <limits> // std::epsilon()

SplineBasis::SplineBasis(const std::vector<double> &r_grid_, int k_spline, int n_spline)
: r_grid(r_grid_), num_spl(n_spline - 3)
{
    std::cout << "> Constructing BSpline basis, k = " << k_spline << ", n = " << num_spl << "... ";
    
    // construct splines
    BSpline splines(k_spline, n_spline, r_grid.front(), r_grid.back()); // k, n, r0, r_max

    bspl = std::vector<std::vector<double>>(num_spl, std::vector<double>(r_grid.size()));
    bspl_derivative = std::vector<std::vector<double>>(num_spl, std::vector<double>(r_grid.size()));

    // reduce number of divisions.
    dr = r_grid.at(1) - r_grid.at(0); // is this correct?

    double h = sqrt(std::numeric_limits<double>::epsilon());

    // TODO - Is splines.b(i) absolutely correct???

    #pragma omp parallel for
    for (int i = 2; i < n_spline - 1; ++i)
    {
        for (std::size_t j = 0; j < r_grid.size(); ++j)
        {
            bspl.at(i-2).at(j) = splines.b(i, r_grid.at(j));
                        
            bspl_derivative.at(i-2).at(j) = (splines.b(i, r_grid.at(j) + h) - splines.b(i, r_grid.at(j) - h))
            / ((r_grid.at(j)+h) - ((r_grid.at(j)-h)));
        }
    }

    std::cout << "done.\n";
}