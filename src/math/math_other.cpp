#include "math_other.hpp"

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