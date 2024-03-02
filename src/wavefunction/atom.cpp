#include "atom.hpp"


// Coulombic potential
auto Atom::construct_potential(int Z, const std::vector<double> &r_grid) -> std::vector<double>
{
    IO::msg::construct("coulombic potential");
    std::vector<double> potential(r_grid.size());

    //std::generate(potential.front(), potential.back(), 
    //    [Z, r_grid](double r) -> double {return - Z / r;} );
    
    for (std::size_t i = 0; i < r_grid.size(); ++i) 
        potential.at(i) = -Z / r_grid.at(i);

    IO::msg::done();

    return potential;
}
