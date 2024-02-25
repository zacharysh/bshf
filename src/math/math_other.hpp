#ifndef MATH_OTHER_HPP_
#define MATH_OTHER_HPP_

#include <vector>

//template <typename T>
auto integrate_trap(std::vector<double> r_grid, std::vector<double> func) -> double;

const std::vector<double> operator*(const std::vector<double> lhs, const std::vector<double> &rhs);
const std::vector<double> operator+(const std::vector<double> lhs, const std::vector<double> &rhs);

#endif