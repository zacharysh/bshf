#include "math_other.hpp"

auto construct_grid_linear(double r_min, double r_max, int n_points) -> std::vector<double>
{
    std::vector<double> r(n_points);

    auto delta_r = (r_max - r_min) / (n_points - 1);

    for(int i = 0; i < n_points; ++i)
        r.at(i) = i * delta_r + r_min;

    return r;
}

auto trapz_linear(double dr, std::vector<double> func) -> double
{
    double result {0.0};

    for (std::size_t i = 0; i < func.size() - 1; ++i)
    {
        result += func.at(i+1) + func.at(i);
    }

    return result * dr / 2.0;
}

auto simpson_linear(double dr, std::vector<double> func) -> double
{
    auto N = func.size();
    auto result {func.at(0) + func.at(N - 1)};

    for (std::size_t i = 1; i < N - 1; ++i)
    {
        auto w = (i % 2 == 0) ? 2.0 : 4.0; // weight
        result += func.at(i) * w;
    }
    return result * dr / 3.0;
}
// r_grid and func need to be vectors of same length
/*
auto integrate_trap(std::vector<double> r_grid, std::vector<double> func) -> double
{
    double result = 0.0;

    // assuming linear grid
    auto dr = r_grid.at(1) - r_grid.at(0);

    for (std::size_t i = 0; i < r_grid.size() - 1; ++i)
    {
        result += (func.at(i+1) + func.at(i)) * dr / 2.0;
    }

    return result;
}
*/

// TODO NAMESPACE FIX UP!

auto operator*(const std::vector<double> &lhs, const std::vector<double> &rhs) -> std::vector<double> 
{
    auto result = lhs;
    
    for(std::size_t i = 0; i < lhs.size(); ++i)
    {
        result.at(i) *= rhs.at(i);
    }
    return result;
}

auto operator+(const std::vector<double> &lhs, const std::vector<double> &rhs) -> std::vector<double>
{
    // assert they're the same size
    auto result = lhs;
    
    for(std::size_t i = 0; i < lhs.size(); ++i)
    {
        result.at(i) += rhs.at(i);
    }
    return result;
}