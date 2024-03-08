#ifndef ATOM_HPP_
#define ATOM_HPP_

#include <vector>
//#include <algorithm> // std::generate
#include <iostream>

#include "../IO/IO.hpp"

#include "../math/basis/spline_basis.hpp"

#include "potential.hpp"
#include "wavefunction.hpp"


class Atom
{
    public:
    int Z;

    SplineBasis basis;
    Potential nuclear_potential;
    Potential interaction_potential;

    std::vector<Electron> electrons {};


    Atom() : Z(), basis() {};

    Atom(int Z_, Potential::Type nuclear_potential_type_, Potential::Type interaction_type_, const SplineBasis &basis_)
    : Z(Z_), basis(basis_),
    nuclear_potential(Potential(Z, nuclear_potential_type_, basis.r_grid)),
    interaction_potential(Potential(Z, interaction_type_, basis.r_grid))
    {
        IO::msg::construct<int>("Atomic system", {{"Z", Z_}}, true);
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

    auto get_range() const -> const std::vector<double>& { return basis.r_grid.range; }
    auto get_range_inv() const -> const std::vector<double>& { return basis.r_grid.range_inv; }
}; // class Atom


#endif