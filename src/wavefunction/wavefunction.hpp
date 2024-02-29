#ifndef WAVEFUNCTION_HPP_
#define WAVEFUNCTION_HPP_

#include <vector>
#include <iostream>

class Wavefunction
{
    public:
    int n;
    int l;
    int m;

    std::vector<double> energy {};
    std::vector<double> expansion_coeffs {};

    Wavefunction(int n_, int l_, int m_);
};

#endif