#ifndef ATOM_HPP_
#define ATOM_HPP_

#include <vector>
//#include <algorithm> // std::generate
#include <iostream>

#include "../io.hpp"

#include "../math/basis/spline_basis.hpp"
#include "../math/basis/spline_basis.hpp"

#include "potential.hpp"
#include "electron.hpp"


class Atom
{
    public:
    int Z;

    int n_max;
    int l_max;

    SplineBasis basis;
    Potential nuclear_potential;
    Potential interaction_potential;

    SquareMatrix<double> kinetic;

    std::vector<Electron> electrons {};

    Atom(int Z_, int n_max_, int l_max_, Potential::Type nuclear_potential_type_, Potential::Type interaction_type_, const SplineBasis &basis_)
    : 
    Z(Z_), 
    n_max(n_max_),
    l_max(l_max_),
    basis(basis_),
    nuclear_potential(Potential(Z, nuclear_potential_type_, basis.r_grid)),
    interaction_potential(Potential(Z, interaction_type_, basis.r_grid)),
    kinetic(basis.num_spl)
    {
        // Compute the kinetic term once, since we don't need to redo it each time.
        // Only need to fill the bottom half since DSYGV_ anticipates a guaranteed symmetric-definite matrix.
        #pragma omp parallel for
        for(int i = 0; i < basis.num_spl; ++i)
        {
            for(int j = 0; j <= i; ++j)
            {
                kinetic(i, j) = 0.5 * simpson_linear(basis.r_grid.dr, basis.bspl_derivative.at(i) * basis.bspl_derivative.at(j));
            }
        }
    }
    
    auto get_energies() -> std::vector<double>
    {
        std::vector<double> energies(electrons.size());
        for(std::size_t i = 0; i < electrons.size(); ++i)
        {
            energies.at(i) = electrons.at(i).energy;
        }
    return energies;
    }


    auto core() -> Electron& { return electrons.front(); }
    auto valence() -> Electron& { return electrons.back(); }

    auto core() const -> const Electron& { return electrons.front(); }
    auto valence() const -> const Electron& { return electrons.back(); }

    auto get_range() const -> const std::vector<double>& { return basis.r_grid.range; }
    auto get_range_inv() const -> const std::vector<double>& { return basis.r_grid.range_inv; }

    auto print_states() const -> void;
}; // class Atom


#endif