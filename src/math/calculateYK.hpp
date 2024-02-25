#pragma once
#include <cassert>
#include <cmath>
#include <vector>

namespace YK {

//******************************************************************************
/*!
Calculates y^k_ab(r)
Method:
y^0_bb(r) = Int_0^infty P(x)P(x) (1/r_max) dx
          = (1/r) Int_0^r P(x)P(x) dx + Int_r^infty P(x)P(x) (1/x) dx
          = (1/r) * a(r) + b(r)
Calculate a and b first, then merge to form y
*/
inline std::vector<double> ykab(int k, const std::vector<double> &Pa,
                                const std::vector<double> &Pb,
                                const std::vector<double> &r) {
  assert(Pa.size() == Pb.size() && "Pa and Pb must have same size in ykab");
  assert(Pb.size() == r.size() && "Pb and r must have same size in ykab");
  assert(r.size() > 3 && "r must be at least 3 points long!");

  const auto n = r.size();
  const auto dr = r.at(1) - r.at(0);
  // test if grid is uniformly spaced
  {
    const auto dr2 = r.at(2) - r.at(1);
    assert(std::abs(dr2 - dr) < 1.0e-12 && "Radial grid must be uniform");
  }

  std::vector<double> y(r.size());
  y[0] = 0.0;
  double a = 0.0;
  for (auto i = 1ul; i < n; ++i) {
    const auto rat = r[i - 1] / r[i];
    auto rho = Pa[i - 1] * Pb[i - 1];
    a = (a + rho / r[i - 1]) * std::pow(rat, k + 1);
    y[i] = a * dr;
  }

  double b = (Pa[n - 1] * Pb[n - 1] / r[n - 1]);
  y[n - 1] += b * dr;
  for (auto i = n - 1; i >= 1; --i) {
    auto rho = Pa[i - 1] * Pb[i - 1];
    const auto rat = r[i - 1] / r[i];
    b = b * std::pow(rat, k) + rho / r[i - 1];
    y[i - 1] += b * dr;
  }

  return y;
}

} // namespace YK
