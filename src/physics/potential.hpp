#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include <vector>
#include <cmath> //exp
#include <algorithm> // std::transform

#include "../io.hpp"

#include "../math/basis/grid.hpp"

class Potential
{
    public:
    enum class Type
    {
        Coulomb,
        Greens,
        HF_Direct,
        HartreeFock,
        Unknown
    };

    public:
    Type type;
    std::vector<double> values;
    
    Potential() : type(Type::Unknown), values(std::vector<double>()) {};
    Potential(int Z, Type type_, const LinearGrid &r_grid);

    static auto Coulomb(int Z, const LinearGrid &r_grid) -> std::vector<double>;
    static auto Greens (int Z, const LinearGrid &r_grid) -> std::vector<double>;
};



#endif