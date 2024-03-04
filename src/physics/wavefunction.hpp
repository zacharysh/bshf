#ifndef WAVEFUNCTION_HPP_
#define WAVEFUNCTION_HPP_

#include <vector>
#include <iostream>

#include "../math/basis/spline_basis.hpp"
#include "../math/math_other.hpp"

#include "../IO/io.hpp"

#include "../physics/physical_constants.hpp"

class Electron
{
    private:


    public:
    int n;
    int l;
    int m;

    // Assume the electron is actually there.
    bool filled {true};

    std::string state_label;

    double energy;

    bool has_correction {false};
    double energy_correction {0.0};
    
    std::vector<double> coeffs;
    std::size_t basis_size;

    std::vector<double> P;
    std::vector<double> amplitude;

    double r1_moment;
    double r2_moment;
    
    Electron()
    : n(), l(), m(), state_label(), energy(), coeffs(), basis_size(), P(), amplitude(), r1_moment(), r2_moment() {};

    Electron(const int n_, const int l_, const int m_, const double energy_, const std::vector<double> coeffs_, const SplineBasis &basis);

    // < r ^ k >
    auto calculate_radial_moment(const LinearGrid &r_grid, int k) -> double;
};

auto calculate_lifetime(Electron a, Electron b, const LinearGrid &r_grid) -> double;

auto print_states(std::vector<Electron> states) -> void;
auto print_state_label(const Electron &psi, bool &excited_state_present) -> void;

#endif