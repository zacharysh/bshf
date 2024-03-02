#ifndef ATOM_HPP_
#define ATOM_HPP_

#include <vector>
//#include <algorithm> // std::generate
#include <iostream>

namespace Atom
{
    auto construct_potential(int Z, const std::vector<double> &r_grid) -> std::vector<double>;

class AtomicSystem
{
    public:
    int Z;

    std::vector<double> potential;

    AtomicSystem(int Z_, const std::vector<double> &r_grid)
    : Z(Z_) //, potential(constructPotential(Z_, r_grid, n_points))
    {
        std::cout << "> Creating Atomic system Z = " << Z_ << ":\n";
        std::cout << "  > Creating Coulomb potential... ";
        potential = construct_potential(Z_, r_grid);
        std::cout << "done.\n";
    }
};

}; // namespace Atom

#endif