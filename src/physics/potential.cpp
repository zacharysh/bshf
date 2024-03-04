#include "potential.hpp"

// This must be fixed ASAP!
Potential::Potential(int Z, Type type_, const LinearGrid &r_grid_)
: m_type(type_), m_values()
{
    switch(type_)
    {
        case Type::Coulomb: { m_values = {Coulomb(Z, r_grid_)}; break; }
        case Type::Greens:  { m_values = {Greens(Z, r_grid_)};  break; }
        case Type::Unknown:
        default:            { m_values = {}; break; }
    }
}

auto Potential::Coulomb(int Z, const LinearGrid &r_grid) -> std::vector<double>
{
    IO::msg::construct<std::string>("potential", {{std::string(), "Coulomb"}});

    std::vector<double> potential(r_grid.grid_size);

    //std::generate(potential.front(), potential.back(), 
    //    [Z, r_grid](double r) -> double {return - Z / r;} );
    
    for (std::size_t i = 0; i < r_grid.grid_size; ++i) 
        potential.at(i) = -Z / r_grid.at(i);

    IO::msg::done();

    return potential;
}

auto Potential::Greens(int Z, const LinearGrid &r_grid) -> std::vector<double>
{
    IO::msg::construct<std::string>("potential", {{std::string(), "Green's"}});

    // Valid for Lithium
    const auto h = 1;
    const auto d = 0.2;

    std::vector<double> potential(r_grid.grid_size);

    for (std::size_t i = 0; i < r_grid.grid_size; ++i)
    {
        auto h_e = h * exp(r_grid.at(i) / d);
        potential.at(i) = (Z - 1) / r_grid.at(i) * (h_e - h) / (1 + h_e - h);
    }

    IO::msg::done();

    return potential;
}