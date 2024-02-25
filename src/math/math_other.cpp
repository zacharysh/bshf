#include "math_other.hpp"

auto construct_grid_linear(double r_min, double r_max, int n_points) -> std::vector<double>
{
    std::vector<double> r(n_points);

    auto delta_r = (r_max - r_min) / (n_points - 1);

    for(int i = 0; i < n_points; ++i) { r.at(i) = i * delta_r; }

    return r;
}

// r_grid and func need to be vectors of same length

auto integrate_trap(std::vector<double> r_grid, std::vector<double> func) -> double
{
    double result = 0.0;

    // assuming linear grid
    auto dr = r_grid.at(1) - r_grid.at(0);

    for (size_t i = 0; i < r_grid.size() - 1; ++i)
    {
        result += (func.at(i+1) + func.at(i)) * dr / 2.0;
    }

    return result;
}

const std::vector<double> operator*(const std::vector<double> lhs, const std::vector<double> &rhs)
{
    auto result = lhs;
    
    for(size_t i = 0; i < lhs.size(); ++i)
    {
        result.at(i) *= rhs.at(i);
    }
    return result;
}
const std::vector<double> operator+(const std::vector<double> lhs, const std::vector<double> &rhs)
{
    auto result = lhs;
    
    for(size_t i = 0; i <= lhs.size(); ++i)
    {
        result.at(i) += rhs.at(i);
    }
    return result;
}