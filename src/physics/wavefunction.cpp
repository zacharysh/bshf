#include "wavefunction.hpp"

Electron::Electron(const int n_, const int l_, const double energy_, const std::vector<double> &coeffs_, const SplineBasis &basis_)
: n(n_), l(l_),
state_label(std::to_string(n)),
energy(energy_),
basis_size(coeffs_.size()),
P(std::vector<double>(basis_.grid_size())),
amplitude(std::vector<double>(basis_.grid_size()))
{
    // Must be a better way
    if( l == 0)
        state_label += "s";
    else if( l == 1)
        state_label += "p";
        
    IO::log(LogType::info, "Constructing electron state", state_label);

    // Multiply the basis vectors by the correct weighting.
    // Removes need for multiple for-loops.
    P = std::inner_product(basis_.bspl.begin(), basis_.bspl.end(), coeffs_.begin(), P,
        [] (auto a, auto b) {return a + b;}, [] (auto &b, auto &c) {return b * c;});

    amplitude = P * basis_.r_grid.range_inv;

    //r1_moment = calculate_radial_moment(basis_.r_grid, 1);
    //r2_moment = calculate_radial_moment(basis_.r_grid, 2);
    
    IO::done();
}

auto Electron::calculate_radial_moment(const LinearGrid &r_grid, int k) const -> double
{   
    auto r_k = r_grid.range;

    // better way?
    for (int i = 1; i < k; ++i)
        r_k = r_k * r_grid.range;

    return simpson_linear(r_grid.dr, P * P * r_k);
}

// We assume a is the only empty state and b is the excited state.
auto calculate_lifetime(const Electron &a, const Electron &b, const LinearGrid &r_grid) -> double
{   
    // Use experimental omega.
    constexpr auto omega_3 = 0.06791 * 0.06791 * 0.06791;

    assert((a.l == 0 && b.l == 1 && a.n == 2 && b.n == 2) && "Only works for 2p -> 2s transition currently.");

    auto Rab = simpson_linear(r_grid.dr, a.P * r_grid.range * b.P);

    auto gamma = 2.0 * Rab * Rab * omega_3 / 3.0 * 1.071e10;

    return 1/gamma;
}