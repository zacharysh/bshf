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

    Atom(int Z_, Potential::Type nuclear_potential_type_, Potential::Type interaction_type_, const SplineBasis &basis_)
    : Z(Z_), basis(basis_)
    {
        IO::msg::construct<int>("Atomic system", {{"Z", Z_}}, true);
        nuclear_potential = Potential::construct_potential(Z, nuclear_potential_type_, basis.r_grid);
        interaction_potential = Potential::construct_potential(Z, interaction_type_, basis.r_grid);
    }
}; // class Atom


#endif