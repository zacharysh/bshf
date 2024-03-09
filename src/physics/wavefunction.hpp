#ifndef WAVEFUNCTION_HPP_
#define WAVEFUNCTION_HPP_

#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>

#include "../math/basis/spline_basis.hpp"
#include "../math/math_other.hpp"

#include "../io.hpp"

#include "../physics/physical_constants.hpp"

class Electron
{
    private:

    public:
    int n;
    int l;

    // Assume the electron is actually there.
    bool filled {true};

    std::string state_label;

    double energy;

    bool has_correction {false};
    double energy_correction {0.0};
    
    //std::vector<double> coeffs;
    std::size_t basis_size;

    std::vector<double> P;
    std::vector<double> amplitude;
    
    Electron()
    : n(), l(), state_label(), energy(), basis_size(), P(), amplitude() {};

    Electron(const int n_, const int l_, const double energy_, const std::vector<double> &coeffs_, const SplineBasis &basis);

    // < r ^ k >
    auto calculate_radial_moment(const LinearGrid &r_grid, int k) const -> double;
};

auto calculate_lifetime(const Electron &a, const Electron &b, const LinearGrid &r_grid) -> double;


#endif