#include "math_other.hpp"

auto trapz_linear(double dr, const std::vector<double> &func) -> double
{
    double result = 0.0;

    for (std::size_t i = 0; i < func.size() - 1; ++i)
    {
        result += func.at(i+1) + func.at(i);
    }

    return 0.5 * result * dr;
}

auto simpson_linear(double dr, const std::vector<double> &func) -> double
{
    auto result = func.front() + func.back();

    for (std::size_t i = 1; i < func.size() - 1; ++i)
    {
        auto w = (i % 2 == 0) ? 2.0 : 4.0;
        result += func.at(i) * w;
    }

    return result * dr / 3.0;
}

// TODO: namespace fix up
auto operator*(const std::vector<double> &lhs, const double &rhs) -> std::vector<double>
{
    std::vector<double> result {lhs};

    std::transform(result.begin(), result.end(), result.begin(),
        [rhs] (auto x) { return (x * rhs); });
    return result;
}

auto operator*(const double &lhs, const std::vector<double> &rhs) -> std::vector<double> 
{
    return rhs * lhs;
}


auto operator+=(std::vector<double> &lhs, const std::vector<double> &rhs) -> std::vector<double>
{
    // assert they're the same size
    // use std::fma
    std::transform(lhs.begin(), lhs.end(), lhs.begin(),
        [n = 0, rhs] (auto x) mutable { return (x + rhs.at(n++)); });

    return lhs;
}

auto operator*=(std::vector<double> &lhs, const std::vector<double> &rhs) -> std::vector<double>
{
    // assert they're the same size
    
    std::transform(lhs.begin(), lhs.end(), lhs.begin(),
        [n = 0, rhs] (auto x) mutable { return (x * rhs.at(n++)); });

    return lhs;
}

auto operator+(const std::vector<double> &lhs, const std::vector<double> &rhs) -> const std::vector<double>
{
    std::vector<double> result = {lhs};
    return result += rhs;
}

auto operator*(const std::vector<double> &lhs, const std::vector<double> &rhs) -> const std::vector<double> 
{
    std::vector<double> result = {lhs};
    return result *= rhs;
}