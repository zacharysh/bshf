#include "potential.hpp"

Potential::Potential(int Z, Type type_, const LinearGrid &r_grid_)
: type(type_), values()
{
    switch(type_)
    {
        case Type::Coulomb: { values = {Coulomb(Z, r_grid_)};                    break; }
        case Type::Greens:  { values = {Greens(Z, r_grid_)};                     break; }
        case Type::Unknown: { values = {std::vector<double>(r_grid_.grid_size)}; break; }
        default:            { values = {};                                       break; }
    }
}

auto Potential::Coulomb(int Z, const LinearGrid &r_grid) -> std::vector<double>
{
    IO::log(LogType::info, "Constructing potential", "Coulomb");

    std::vector<double> potential(r_grid.grid_size);

    std::transform(r_grid.range_inv.begin(), r_grid.range_inv.end(), potential.begin(),
        [Z](auto r_inv) {return -Z * r_inv;} );

    IO::done();

    return std::vector<double>(potential);
}

auto Potential::Greens(int Z, const LinearGrid &r_grid) -> std::vector<double>
{
    IO::log(LogType::info, "Constructing potential", "Green's");

    // Valid for Lithium.
    constexpr auto h = 1;
    constexpr auto d = 0.2;

    std::vector<double> potential(r_grid.grid_size);
    
    std::transform(r_grid.range.begin(), r_grid.range.end(), potential.begin(),
        [Z, h, d](auto r) {
            auto h_e = h * exp(r / d);
            return (Z - 1) * (h_e - h) / (1 + h_e - h) / r;
        });

    IO::done();

    return std::vector<double>(potential);
}