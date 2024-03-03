#pragma once
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

namespace Range {

// Returns a range of 'number' values, spread out between 'first' and 'last' on
// a scale according to function f. First and last of range is guarenteed to be
// those given
template <typename Real, typename Function>
std::vector<Real> F_range(Real first, Real last, int number, Function f) {
  assert(number >= 2);
  std::vector<Real> range;
  range.reserve(static_cast<std::size_t>(number));
  range.push_back(first); // guarentee first is first
  for (int i = 1; i < number - 1; ++i) {
    const auto eps = static_cast<Real>(i) / (number - 1);
    const auto value = f(first, eps); // first + (eps * interval);
    range.push_back(value);
  }
  range.push_back(last); // guarentee last is last
  return range;
}

// Returns a range of 'number' values, spread out between 'first' and 'last' on
// a logarithmic scale. First and last of range is guarenteed to be those given
template <typename Real>
std::vector<Real> logarithmic(Real first, Real last, int number) {
  const auto log_ratio = std::log(static_cast<double>(last / first));
  return F_range(first, last, number, [log_ratio](Real first, Real eps) {
    return first *
           static_cast<Real>(std::exp(log_ratio * static_cast<double>(eps)));
  });
}

// Returns a range of 'number' values, spread out between 'first' and 'last' on
// a uniform scale. First and last of range is guarenteed to be those given
template <typename Real>
std::vector<Real> uniform(Real first, Real last, int number) {
  const auto interval = last - first;
  return F_range(first, last, number, [interval](Real first, Real eps) {
    return first + (eps * interval);
  });
}

} // namespace Range

//******************************************************************************
class BSpline {

  //****************************************************************************
private:
  int K;                 // order of B-splines (have degraa K-1)
  int N;                 // Number of B-splines
  double r0;             // First internal knot
  double rmax;           // Last internal knot
  bool loqQ;             // Space knots logarithmically?
  std::vector<double> t; // knot vector

  // nb: A spline of order K is a piecewise polynomial function of degree K-1
  // So, cubic splines have K=4

  //****************************************************************************
public:
  BSpline(int k, int n, double in_r0, double in_rmax, bool logarithmic = false)
      : K(k), N(n), r0(in_r0), rmax(in_rmax), loqQ(logarithmic),
        t(set_knots()) {}

  // ith spline: B_i(r)
  double b(int i, double r) const {
    assert(i >= 0 && i < N && "Spline index range must be i = [0,N)");
    return bki(K, i, r);
  }

  // overload of (); same as b(i,r)
  double operator()(int i, double r) const { return b(i, r); }

  // Prints list of knots to screen in nice format
  void print_knots() const {
    std::cout << "Spline Knots:\n";
    int count = 0;
    for (auto tk : t) {
      ++count;
      printf("%7.2e", tk);
      if (count == (int)t.size()) {
        std::cout << "\n";
      } else {
        std::cout << ", ";
      }
      if (count % 7 == 0 && count != (int)t.size())
        std::cout << "\n";
    }
  }

  //****************************************************************************
private:
  static bool fequal(double a, double b, double eps = 1.0e-12) {
    return std::abs(a - b) < eps;
  }

  //----------------------------------------------------------------------------
  // Sets up the lst of knots
  std::vector<double> set_knots() const {
    std::vector<double> tmp_t;

    // The 1st spline (index 0) is non-zero only below r0
    // The (k+1)th spline (and above) [index k] are non-zero only above r0
    // Only the nth spline (index n-1) is non-zero at rmax

    // n_mid is number of "internal" (middle) knots; does not include zero, does
    // include end-point (but only once)
    const auto n_mid = N - K + 1;
    assert(K > 0 && "Require spline order K > 0");
    assert(n_mid > K && "Require #(internal) N - K + 1 > K; increase N");

    // Start and end of knots extend _slightly_ beyond requested range.
    // Ensures we have good behaviour at exactly r=r0 and r=rmax
    const auto eps = 1.0e-8; // 1.0e-4; // 0.0;
    const auto r0_eff = (1.0 - eps) * r0;
    const auto rmax_eff = (1.0 + eps) * rmax;

    // First knot (at r=0.0) repeated K times
    for (int i = 0; i < K; ++i) {
      tmp_t.push_back(0.0);
    }

    // Mid (internal) knots: spaced logarithmically
    const auto mid_knots = loqQ ? Range::logarithmic(r0_eff, rmax_eff, n_mid)
                                : Range::uniform(r0_eff, rmax_eff, n_mid);
    // .insert fails on smp-teaching??
    // tmp_t.insert(t.end(), mid_knots.begin(), mid_knots.end());
    for (auto tt : mid_knots) {
      tmp_t.push_back(tt);
    }

    // Last knot (at r=rmax) repeated K times (alrady once above)
    for (int i = 0; i < K - 1; ++i) {
      tmp_t.push_back((1.0 + eps) * rmax);
    }

    return tmp_t;
  }

  //----------------------------------------------------------------------------
  // base B-spline
  double b1(int i, double x) const {
    return (t.at(i) <= x && x <= t.at(i + 1)) ? 1.0 : 0.0;
  }

  //----------------------------------------------------------------------------
  // B-spline "weights" (used in recursion formula)
  double wki(int k, int i, double x) const {
    // if (std::abs(t.at(i + k) - t.at(i)) < 1.0e-10)
    if (fequal(t.at(i + k), x))
      return 1.0;
    if (fequal(t.at(i + k), t.at(i)))
      return 0.0;
    return (x - t.at(i)) / (t.at(i + k) - t.at(i));
  }

  //----------------------------------------------------------------------------
  // B-spline recursion formula
  double bki(int k, int i, double x) const {
    if (k == 1)
      return b1(i, x);
    const auto w1 = wki(k - 1, i, x);
    const auto w2 = wki(k - 1, i + 1, x);
    return w1 * bki(k - 1, i, x) + (1.0 - w2) * bki(k - 1, i + 1, x);
  }
};
