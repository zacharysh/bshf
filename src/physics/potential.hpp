#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include <vector>
#include <cmath> //exp

#include "../IO/IO.hpp"

#include "../math/basis/grid.hpp"

class Potential
{   
    public:
    enum class Type
    {
        Coulomb,
        Greens,
        HartreeFock,
        Unknown
    };

    private:
    Type m_type;
    std::vector<double> m_values;
    
    public:
    Potential() : m_type(Type::Unknown), m_values(std::vector<double>()) {};
    Potential(Type type_, std::vector<double> values_) : m_type(type_), m_values(values_) {};
    
    static auto Coulomb(int Z, const LinearGrid &r_grid) -> Potential;
    static auto Greens (int Z, const LinearGrid &r_grid) -> Potential;

    static auto construct_potential(int Z, Potential::Type type_, const LinearGrid &r_grid) -> Potential;

    auto get_values() const -> const std::vector<double> { return m_values; };
};



#endif