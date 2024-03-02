#ifndef MATH_OTHER_HPP_
#define MATH_OTHER_HPP_

#include <vector>

#include "../IO/io.hpp"

//template <typename T>
auto trapz_linear(double dr, std::vector<double> func) -> double;
auto simpson_linear(double dr, std::vector<double> func) -> double;

auto construct_grid_linear(double r_min, double r_max, int n_points) -> std::vector<double>;


//TODO FIX NAMESPACE
// also. const? or not?
auto operator*(const std::vector<double> &lhs, const std::vector<double> &rhs) -> std::vector<double>;
auto operator+(const std::vector<double> &lhs, const std::vector<double> &rhs) -> std::vector<double>;

#endif