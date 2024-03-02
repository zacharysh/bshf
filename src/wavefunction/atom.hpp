#ifndef ATOM_HPP_
#define ATOM_HPP_

#include <vector>
//#include <algorithm> // std::generate
#include <iostream>

#include "../IO/io.hpp"

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
        
        IO::msg::construct<int>("Atomic system", {{"Z", Z_}}, true);
        potential = construct_potential(Z_, r_grid);
    }
};

}; // namespace Atom

#endif