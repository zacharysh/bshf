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
        HF_Direct,
        HF_Exchange,
        HartreeFock,
        Unknown
    };

    //private:

    public:
    Type m_type;
    std::vector<double> m_values;
    
    Potential() : m_type(Type::Unknown), m_values(std::vector<double>()) {};
    Potential(int Z, Type type_, const LinearGrid &r_grid);


    //auto operator=(const Potential &) -> Potential& = default;
    //auto operator=(Potential &&) -> Potential& = default;
    
    static auto Coulomb(int Z, const LinearGrid &r_grid) -> std::vector<double>;
    static auto Greens (int Z, const LinearGrid &r_grid) -> std::vector<double>;

    auto get_values() const -> const std::vector<double> { return m_values; };
    auto get_type() -> Type { return m_type; }; // const ?


    auto set_type(Type new_type_) -> void { (*this).m_type = new_type_;};
    auto set_values(std::vector<double> new_values_) -> void { (*this).m_values = new_values_; };
};



#endif