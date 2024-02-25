#ifndef ATOM_HPP_
#define ATOM_HPP_

#include <vector>
//#include <algorithm> // std::generate
#include <iostream>

#include "init.hpp"

namespace Atom
{
    auto constructPotential(int Z, std::vector<double> r_grid, int n_points) -> std::vector<double>;

class AtomicSystem
{
    public:
    int Z;

    std::vector<double> potential;

    AtomicSystem(int Z_, std::vector<double> r_grid)
    : Z(Z_) //, potential(constructPotential(Z_, r_grid, n_points))
    {
        std::cout << "Creating Atomic system Z = " << Z_ << ".\n";
        std::cout << "Creating Coulomb potential... ";
        potential = constructPotential(Z_, r_grid, static_cast<int>(r_grid.size()));
        std::cout << "Done.\n";
    }
};

}; // namespace Atom

#endif