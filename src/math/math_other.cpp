#include "math_other.hpp"

auto trapz_linear(double dr, const std::vector<double> &func) -> double
{
    double result = 0.0;

    for (std::size_t i = 0; i < func.size() - 1; ++i)
    {
        result += func.at(i+1) + func.at(i);
    }

    return double(0.5 * result * dr);
}

auto simpson_linear(double dr, const std::vector<double> &func) -> double
{
    auto result = func.front() + func.back();
    //auto result = func.back();

    for (std::size_t i = 1; i < func.size() - 1; ++i)
    //for (std::size_t i = 1; i < func.size(); ++i)// we start at func.front() since our integral has to have lower bound at r=0.
    {
        auto w = (i % 2 == 0) ? 2.0 : 4.0;
        result += func.at(i) * w; // we start at func.front() since our integral has to have lower bound at r=0.
    }

    return double(result * dr / 3.0);
}


auto ykab(const int k, const std::vector<double> &Pa,  const std::vector<double> &Pb, const LinearGrid &r_grid) -> std::vector<double>
{
    const auto N = r_grid.grid_size;

    std::vector<double> y(N);

    double a = 0.0;

    for (std::size_t i = 1; i < N; ++i)
    {
        const auto rat = r_grid.at(i - 1) * r_grid.inv_at(i);
        a = (a + Pa.at(i - 1) * Pb.at(i - 1) * r_grid.inv_at(i-1)) * pow(rat, k + 1);

        y.at(i) = a * r_grid.dr;
    }

    double b = Pa.at(N - 1) * Pb.at(N - 1) * r_grid.inv_at(N - 1);

    y.at(N - 1) += b * r_grid.dr;

    for (auto i = N - 1; i >= 1; --i)
    {
        auto rho = Pa.at(i - 1) * Pb.at(i - 1);
        const auto rat = r_grid.at(i - 1) * r_grid.inv_at(i);
        b = b * pow(rat, k) + rho * r_grid.inv_at(i - 1);
        y.at(i - 1) += b * r_grid.dr;
    }

    return std::vector(y);
}

// TODO: namespace fix up
auto operator*(const std::vector<double> &lhs, const double &rhs) -> std::vector<double>
{
    std::vector<double> result {lhs};

    std::transform(result.begin(), result.end(), result.begin(),
        [rhs] (auto x) { return x * rhs; });
    return std::vector<double>(result);
}

auto operator*(const double &lhs, const std::vector<double> &rhs) -> std::vector<double> 
{
    return rhs * lhs;
}


auto operator+=(std::vector<double> &lhs, const std::vector<double> &rhs) -> std::vector<double>
{
    // Should use std::fma...
    std::transform(lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(),
        [] (auto x, auto y) { return x + y; });

    return lhs;
}

auto operator*=(std::vector<double> &lhs, const std::vector<double> &rhs) -> std::vector<double>
{
    std::transform(lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(),
        [] (auto x, auto y) { return x * y; });

    return lhs;
}

auto operator+(const std::vector<double> &lhs, const std::vector<double> &rhs) -> const std::vector<double>
{
    std::vector<double> result = {lhs};
    return std::vector<double>(result += rhs);
}

auto operator*(const std::vector<double> &lhs, const std::vector<double> &rhs) -> const std::vector<double> 
{
    std::vector<double> result = {lhs};
    return std::vector<double>(result *= rhs);
}