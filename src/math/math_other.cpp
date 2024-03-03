#include "math_other.hpp"

auto trapz_linear(double dr, std::vector<double> func) -> double
{
    double result {0.0};

    for (std::size_t i = 0; i < func.size() - 1; ++i)
    {
        result += func.at(i+1) + func.at(i);
    }

    return result * dr / 2.0;
}

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