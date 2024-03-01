#ifndef MATH_OTHER_HPP_
#define MATH_OTHER_HPP_

#include <vector>

//template <typename T>
auto trapz_linear(double dr, std::vector<double> func) -> double;
auto simpson_linear(double dr, std::vector<double> func) -> double;

auto construct_grid_linear(double r_min, double r_max, int n_points) -> std::vector<double>;

const std::vector<double> operator*(const std::vector<double> lhs, const std::vector<double> &rhs);
const std::vector<double> operator+(const std::vector<double> lhs, const std::vector<double> &rhs);

#endif