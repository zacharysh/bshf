#ifndef PHYSICAL_CONSTANTS_HPP_
#define PHYSICAL_CONSTANTS_HPP_

namespace PhysicalConstants
{
    
constexpr auto Eh_to_eV(const double E) -> double
{
    return E * 27.2114079527; // citation?
}

}; //namespace PhysicalConstants

#endif