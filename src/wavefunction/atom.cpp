#include "atom.hpp"


// Coulombic potential
auto Atom::constructPotential(int Z, std::vector<double> r_grid, int n_points) -> std::vector<double>
{
    std::vector<double> potential {};//(n_points);
    potential.reserve(n_points);

    //std::generate(potential.front(), potential.back(), 
    //    [Z, r_grid](double r) -> double {return - Z / r;} );
    
    for (auto r : r_grid)
        potential.push_back(- Z / r);

    return potential;
}
