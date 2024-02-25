#include "init.hpp"

auto constructGridLin(double r_min, double r_max, int n_points) -> std::vector<double>
{
    std::vector<double> r(n_points);

    auto delta_r = (r_max - r_min) / (n_points - 1);

    for(int i = 0; i < n_points; ++i) { r.at(i) = i * delta_r; }

    return r;
}