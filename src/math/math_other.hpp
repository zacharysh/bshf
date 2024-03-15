#ifndef MATH_OTHER_HPP_
#define MATH_OTHER_HPP_

#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>

#include "../io.hpp"

#include "basis/grid.hpp"

auto trapz_linear(double dr, const std::vector<double> &func) -> double;
auto simpson_linear(double dr, const std::vector<double> &func) -> double;

auto ykab(const int k, const std::vector<double> &Pa, const std::vector<double> &Pb, const LinearGrid &r_grid) -> std::vector<double>;

auto operator*(const std::vector<double> &lhs, const double &rhs)                   -> std::vector<double>;
auto operator*(const double &lhs, const std::vector<double> &rhs)                   -> std::vector<double>;
auto operator+=(std::vector<double> &lhs, const std::vector<double> &rhs)           -> std::vector<double>;
auto operator*=(std::vector<double> &lhs, const std::vector<double> &rhs)           -> std::vector<double>;
auto operator+(const std::vector<double> &lhs, const std::vector<double> &rhs)      -> const std::vector<double>;
auto operator*(const std::vector<double> &lhs, const std::vector<double> &rhs)      -> const std::vector<double>;

#endif