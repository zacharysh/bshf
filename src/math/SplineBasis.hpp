#ifndef SPLINE_BASIS_HPP_
#define SPLINE_BASIS_HPP_

#include "bspline.hpp"

#include <vector>

class SplineBasis
{
    public:
    int num_spl;
    std::vector<double> r_grid;
    std::vector<std::vector<double>> bspl {};
    std::vector<std::vector<double>> bspl_derivative {};

    SplineBasis(std::vector<double> r_grid_, int k_spline, int n_spline);
};

#endif